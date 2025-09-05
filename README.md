# astrolabe <img src="man/figures/logo.png" align="left" width="120" />

**astrolabe** is an R package for inferring **non-linear causal relations** among variables.  
It couples random forests with an information-theoretic objective: the entropy of residuals. Higher scores imply “cleaner” residuals on the chosen predictors — a proxy for causality.

---

## 🔬 Method in a nutshell

For each candidate relationship, **astrolabe** fits a tuned random forest, computes residuals, and evaluates a **Kozachenko–Leonenko** k-NN entropy on the joint space of predictors and residual. It then drills down by pruning low-importance predictors, and validates findings via bootstrap and permutation. The full pipeline can also add a **pairwise direction matrix** for 2-variable edges and merge both sources into a single DAG. Given a model $m: \ \ Y\sim X_1 + X_2+X_3+\dots+X_p$, with $p$ the number of predictors, the core metric on which the method lies is: 

$$
H(m) = e^{-S_{\text{knn}}(m)} ,
$$

where the $S_{knn}(m)$ is defined as follows:

$$
S_{\text{knn}}(m) \approx \psi(n) - \psi(k) + \log(c_d) + \frac{d}{n} \sum_{i=1}^n \log \varepsilon_i ,
$$

---

## 📦 Installation

```r
# install.packages("remotes") # if needed
remotes::install_github("astrolabe9698/astrolabe")
```

---

## 🧪 Quick start

```r
library(astrolabe)

# load data
df <- your_data

# End-to-end pipeline (multi-outcome scan + robust validation + pairwise)
out <- complete_function(
  df,
  which = "all",       # "robust" or "binary" also allowed
  n_boot = 100,
  n_perm = 50,
  ntree  = 500,
  importance_method = "neg_exp",
  plot = TRUE
)

```

---

## 📘 Mini tutorial — core functions

### 1) Pairwise direction matrix: `cor_forest_matrix_robust_perm()`

**What it does.**  
For every unordered pair (X, Y), it decides the preferred direction on real data, then estimates robustness with bootstrap and an **empirical permutation test**. Returns:

- `real`: chosen directions on real data (binary adjacency)  
- `freq`: bootstrap hit counts  
- `significant`: directions surviving permutation at `alpha`  
- `pval`: empirical p-values  

**Example.**
```r
mats <- cor_forest_matrix_robust_perm(
  df,
  B = 50, P = 100,
  ntree = 300,
  alpha = 0.05,
  importance_method = "neg_exp",
  verbose = TRUE
)
which(mats$significant == 1, arr.ind = TRUE)
```

---

### 2) Robust validation: `robust_scan_all_outcomes()`

**What it does.**  
Validates the **multi-outcome scan results** with:

1. **Bootstrap on real data** → how often predictor sets reappear  
2. **Permutation test** (with embedded bootstrap) → empirical p-value  

**Example.**
```r
scan_res <- scan_all_outcomes_complete(
  df, ntree = 500,
  importance_method = "neg_exp",
  verbose = TRUE
)

robust <- robust_scan_all_outcomes(
  df, scan_results = scan_res,
  n_boot = 300, n_perm = 50,
  alpha = 0.05,
  ntree = 500,
  importance_method = "neg_exp"
)

robust[[1]]$p_empirical
robust[[1]]$significant
```

---

### 3) Complete pipeline: `complete_function()`

**What it does.**  
Runs everything in sequence:

1. Preprocessing (numeric/factor coercion, rare-level handling)  
2. Multi-outcome scan with drill-down pruning  
3. Robust validation (bootstrap + permutation)  
4. Optional pairwise direction matrix  
5. Merge into a single DAG (`pairwise`, `complex`, or `both`)  

**Example.**
```r
out <- complete_function(
  df,
  which = "all",
  n_boot = 100, n_perm = 50,
  alpha = 0.05,
  ntree = 500,
  importance_method = "neg_exp",
  plot = TRUE, layout = "kk"
)

out$scan_res   # raw scan results per outcome
out$robust     # bootstrap+permutation validation
out$binary     # pairwise direction matrices
out$res_all    # merged edge list
out$graph      # ggplot DAG
```

---

## 🧭 Choosing `importance_method`

- `"fixed"` — static thresholds for numeric vs categorical predictors  
- `"neg_exp"` — adaptive thresholds decaying with dimensionality (recommended)  
- `"net_clust"` — DBSCAN + small neural net classifier on the “high” cluster(s)  

---

## ⚙️ Practical tips

- Always `set.seed()` for reproducibility  
- Use `n_cores` to parallelize (defaults to `detectCores()-1`)  
- Plot layouts: `"kk"`, `"fr"`, `"sugiyama"`  
- Mixed data types are handled automatically (numeric/factors)
- Use `always_predictios` to set as predictor a variable that can never be an outcome (for example time in a longitudinal dataset)

---

## 📚 Citation & License

See the package manual for full references and function details.  
MIT License.
