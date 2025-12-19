# Sequential Indicator Simulation (SIS) for Categorical Vegetation Classes (per `strManagementUnit`)

This workflow generates, **per management unit**, (1) **probability rasters** for each vegetation class and (2) a **single categorical raster** representing the most likely class (with an optional “Uncertain” label and MMU cleanup). It also produces a small **diagnostics bundle** to document model selection and key parameters.

Target categorical field: `strAveHeight_cm_PctCover` with classes:

- Short Open
- Tall Open
- Mid Mod
- Short Dense
- Tall Dense

Coordinates are assumed to be **UTM Zone 10, NAD83, meters** (EPSG:26910) in fields `UTM_X`, `UTM_Y`.

---

## 1. Conceptual overview

### Why SIS?
With categorical outcomes and no covariates, **indicator geostatistics** is a standard approach for producing **class probability surfaces** and **stochastic realizations** that reflect spatial uncertainty. **Sequential Indicator Simulation (SIS)** produces multiple conditional realizations that (a) reproduce local class proportions, (b) reflect spatial continuity, and (c) allow uncertainty quantification via ensembles (probabilities, confidence).

- SIS provides probability maps as the **fraction of realizations** taking each class at each cell.
- It provides a categorical map by selecting the class with **maximum probability** at each cell (optionally applying uncertainty and MMU filtering).

**Key references (indicator geostatistics & SIS):**
- Journel, A. G. (1983). *Nonparametric estimation of spatial distributions.* Mathematical Geology.  
- Journel, A. G., & Huijbregts, C. J. (1978). *Mining Geostatistics.* Academic Press.  
- Goovaerts, P. (1997). *Geostatistics for Natural Resources Evaluation.* Oxford University Press.  
- Deutsch, C. V., & Journel, A. G. (1998). *GSLIB: Geostatistical Software Library and User’s Guide.* Oxford University Press.

---

## 2. Unit-based processing

All steps are run independently **per `strManagementUnit`** to avoid blending potentially different spatial structures and class prevalences across study areas.

Each unit is processed as:

1. Read points for the unit
2. Print unit summary (class counts + nearest-neighbor stats + derived declustering and block sizes)
3. Define grid and spatial domain (envelope)
4. Compute declustering weights (to reduce clustered sampling bias)
5. Select a shared variogram family + parameters (per unit) via spatial block CV
6. Run SIS to produce realizations → probability rasters → categorical raster
7. Apply uncertainty labeling and MMU cleanup
8. Save diagnostics (variograms, CV scores, fit surface, realization panel)

---

## 3. Spatial domain (“envelope”) mask

Rather than using a bounding box, the workflow uses a **convex hull envelope** of sample points per unit, optionally buffered by `ENVELOPE_BUFFER_M`. Outputs are masked outside this domain.

**Rationale:** Prevents extrapolation into areas that were never sampled and keeps raster extents tight to the sampling footprint (a common best practice for sample-footprint mapping).

Implementation notes:
- Preferred: `shapely` convex hull + buffer
- Fallback: `scipy.spatial.ConvexHull`

---

## 4. Declustering weights

Field sampling often produces **clustered points** (e.g., observers map presences densely once found). Declustering reduces the undue influence of dense clusters on estimated proportions and variograms.

This workflow uses simple **cell declustering**:
- Define a declustering grid cell size = `DECLUSTER_MULT × median NN distance`, clamped to `[DECLUSTER_CELL_MIN_M, DECLUSTER_CELL_MAX_M]`
- Weight each point as `w_i = 1 / n_cell(i)`
- Normalize weights to mean 1.0

Declustered weights are used to estimate:
- Weighted class prevalence `p_k`
- Total indicator sill per class: `sill_k = p_k (1 - p_k)`
- Empirical variogram for the binary field used in variogram selection (see below)

**References (declustering / sampling bias):**
- Isaaks, E. H., & Srivastava, R. M. (1989). *An Introduction to Applied Geostatistics.* Oxford University Press.  
- Deutsch, C. V. (2002). *Geostatistical reservoir modeling.* Oxford University Press. (declustering concepts widely used)

---

## 5. Variogram family selection with spatial block CV

### 5.1 Why a shared (range, nugget) per unit?
Fitting a full variogram for every class can be unstable when some classes are rare. Instead, the workflow fits **one shared range and nugget per unit**, and allows **class-specific sills** driven by declustered prevalence.

This improves stability while still letting each class have different variance magnitude.

### 5.2 What is actually cross-validated?
To select the variogram family (exponential / spherical / gaussian) and shared parameters, the workflow uses a **binary field**:

- Let the **dominant class** be the class with highest decluster-weighted prevalence in the unit.
- Define `J = 1` if a point is **NOT** dominant, else `J = 0`.

This provides a robust binary signal of “dominant vs rest” with enough support to fit.

Then:
- Compute a decluster-weighted empirical variogram of `J`
- Fit exp/sph/gau models via grid search over (range, nugget)
- Score each model family using **spatial block cross-validation** and **log loss** on held-out blocks using **local ordinary kriging** predictions for `J`.

Block size = `BLOCK_SIZE_MULT × median NN distance`, clamped to `[BLOCK_SIZE_MIN_M, BLOCK_SIZE_MAX_M]`.

**Rationale:**  
Spatial CV reduces optimistic bias from spatial autocorrelation; block CV is a common practical strategy for spatial model validation.

**References (spatial cross-validation):**
- Roberts, D. R., et al. (2017). *Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure.* Ecography.  
- Brenning, A. (2012). *Spatial cross-validation and bootstrap for the assessment of prediction rules in remote sensing.* (methodological precedent for spatial CV)

---

## 6. Class-specific sills with shared nugget + range

For each class `k` (1..5):
- Estimate decluster-weighted prevalence `p_k`
- Total sill: `sill_k = p_k (1 - p_k)`
- Shared nugget from the selected unit model is applied to all classes:
  - Partial sill: `partial_k = max(sill_k - nugget, 0)`

Optional: A sill floor `SILL_FLOOR` can prevent extremely small partial sills for very rare classes.

**Rationale:**  
Indicator variance is naturally `p(1-p)`; using this links sill magnitude to prevalence and avoids overfitting sills from sparse data.

**References (indicator sill):**
- Goovaerts (1997) – indicator coding and sill relationships  
- Deutsch & Journel (1998) – indicator kriging and SIS practice

---

## 7. Sequential Indicator Simulation (SIS)

### 7.1 Conditioning data
- Start with observed points inserted into the raster grid (fixed hard data).
- During simulation, each newly simulated cell becomes additional conditioning data.

### 7.2 Neighborhood search
For each unsimulated cell inside the envelope:
- Find neighbors within `search_radius = SEARCH_RADIUS_FACTOR × range`
- Use up to `MAX_NEIGHBORS` and require at least `MIN_NEIGHBORS` for local kriging; otherwise revert to marginal `p_k`.

### 7.3 Conditional probabilities
At each cell, compute probabilities for each class:

1. Compute ordinary kriging weights once per cell from neighbor geometry using a reference covariance scale.
2. For each class, compute local indicator mean using neighbor class indicators.
3. Apply a shrinkage blend toward the marginal `p_k`:
   - `lambda_k = partial_k / (partial_k + nugget)`
   - `P_k = lambda_k * m_k + (1 - lambda_k) * p_k`
4. Normalize probabilities across classes.
5. Sample the class according to `P_k`.

This produces one conditional realization. Repeat for `N_REALIZATIONS`.

**Why shrinkage?**  
It stabilizes probabilities for rare classes and avoids overly confident local estimates when partial sill is tiny.

---

## 8. Final raster products

### 8.1 Probability rasters (per class)
For each class k:
- `prob_<ClassName>.tif` = (# realizations assigned to class k) / N_REALIZATIONS

### 8.2 Confidence raster
- `p_max_confidence.tif` = max_k(prob_k)

### 8.3 Categorical raster
- `categorical_winner.tif` = argmax_k(prob_k), encoded as:
  - 0 = Uncertain / NoData (masked or below confidence threshold)
  - 1..5 = class codes (see metadata)

### 8.4 Uncertain label (optional)
Cells with `p_max < UNCERTAIN_TOLERANCE` are set to 0 (“Uncertain”).

### 8.5 MMU cleanup (optional)
A minimum mapping unit is enforced on the categorical raster:
- Patches smaller than `MMU_CELLS` are removed and filled using local neighborhood majority filtering.
- Requires `scipy.ndimage`.

---

## 9. Diagnostics bundle (per unit)

The workflow writes the following diagnostics in each unit folder:

1. **Empirical variogram + fitted curves**  
   `variogram_J_fit.png`  
   - Declust-weighted empirical variogram of J
   - Fitted curves for exponential, spherical, gaussian (all-data fit)

2. **Spatial block CV scores**  
   `cv_scores.csv`  
   `cv_scores.png`  
   - Per-fold log loss by variogram family

4. **Variogram fit surface for selected family**  
   `variogram_fit_surface_<best_model>.png`  
   - Heatmap of weighted SSE over (range, nugget) grid
   - Best-fit marked

5. **Realizations panel**  
   `realizations_panel.png`  
   - A small selection of realizations shown as categorical maps

Additional:
- `metadata.json` with parameters, classes, chosen variogram, etc.

---

## 10. Inputs

### Required CSV fields
- `strManagementUnit` (unit identifier)
- `strAveHeight_cm_PctCover` (class label; must be one of the five supported classes)
- `UTM_X`, `UTM_Y` (meters, EPSG:26910)

### Parameters (key)
- `GRID_CELL_SIZE_M`: raster resolution
- `ENVELOPE_BUFFER_M`: buffer around convex hull envelope
- `DECLUSTER_MULT` + clamps: declustering cell size
- `BLOCK_SIZE_MULT` + clamps: spatial block CV fold geometry
- `SEARCH_RADIUS_FACTOR`: neighborhood radius as a fraction of selected range
- `N_REALIZATIONS`: number of SIS realizations
- `UNCERTAIN_TOLERANCE`: minimum confidence for categorical assignment
- `MMU_CELLS`: minimum mapping unit size in cells

---

## 11. Outputs (directory layout)

`OUTPUT_DIR/<unit_name>/`
- `prob_Short_Open.tif`
- `prob_Tall_Open.tif`
- `prob_Mid_Mod.tif`
- `prob_Short_Dense.tif`
- `prob_Tall_Dense.tif`
- `p_max_confidence.tif`
- `categorical_winner.tif`
- `metadata.json`
- `variogram_J_fit.png`
- `cv_scores.csv`
- `cv_scores.png`
- `variogram_fit_surface_<exponential|spherical|gaussian>.png`
- `realizations_panel.png`

`OUTPUT_DIR/`
- `sis_summary_by_unit.csv`

---

## 12. Requirements

### Python
- Python 3.9+ recommended

### Required packages
- `numpy`
- `pandas`
- `rasterio`
- `matplotlib`

### Optional (strongly recommended)
- `scipy` (for KD-tree NN distances and MMU cleanup)
- `shapely` (for robust convex hull + buffering)
- `scikit-learn` (fallback for nearest neighbor queries if scipy is missing)

---

## 13. Practical notes & cautions

- **Rare classes:** If a class is extremely rare in a unit, its probability surface may be dominated by marginal prevalence; this is expected without covariates.
- **Stationarity:** The workflow assumes within-unit stationarity for spatial continuity. Strong internal gradients may require splitting the unit or using covariates.
- **Computation:** SIS can be expensive at fine cell sizes and many realizations. Start with coarser grids and fewer realizations, then refine.
- **Uncertainty threshold:** `UNCERTAIN_TOLERANCE` controls how conservative the categorical map is. Higher values produce more “Uncertain” cells.

---

## 14. How to run

In a Jupyter notebook:

1) Set `INPUT_CSV` and `OUTPUT_DIR` in the script header  
2) Adjust parameters as needed  
3) Run:

```python
summary_df = run_sis_workflow(INPUT_CSV, OUTPUT_DIR, seed=RANDOM_SEED)
summary_df
