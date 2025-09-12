# astrolabe <img src="inst/figures/logo.png" align="left" width="120" />

**astrolabe** is an R package for inferring **non-linear causal relations** among variables.  
It couples random forests with an information-theoretic objective: the entropy of residuals. Higher scores imply “cleaner” residuals on the chosen predictors — a proxy for causality.

---

## 🔬 Method in a nutshell

For each candidate relationship, **astrolabe** fits a tuned random forest, computes residuals, and evaluates a **Kozachenko–Leonenko** k-NN entropy on the joint space of predictors and residual. It then drills down by pruning low-importance predictors, and validates findings via bootstrap and permutation. The full pipeline can also add a **pairwise direction matrix** for 2-variable edges and merge both sources into a single DAG. Given a model $m: \ \ Y\sim X_1 + X_2+X_3+\dots+X_p$, with $p$ the number of predictors, the core metric on which the method lies is: 

$$
H(m) = e^{-S_{knn}(m)} ,
$$

where the $S_{knn}(m)$ is defined as follows:

$$
S_{knn}(m) \approx \psi(n) - \psi(k) + \log(c_d) + \frac{d}{n} \sum_{i=1}^n \log \varepsilon_i ,
$$

The $H$-score allows to compare different models with one another giving the possibility to choose the best one.

---

## 📦 Installation

```r
# install.packages("remotes") # if needed
remotes::install_github("astrolabe9698/astrolabe")
```

---

## 📘 Mini tutorial — core functions

### 1. Complete pipeline: `complete_function()`

**What it does.**  
Runs everything in sequence:

1. Preprocessing (numeric/factor coercion, rare-level handling)  
2. Multi-outcome scan with drill-down pruning  
3. Robust validation (bootstrap + permutation)  
4. Optional pairwise direction matrix  
5. Merge into a single DAG (`pairwise`, `complex`, or `both`)  

**Example.**
```r
library(astrolabe)

# load data
df <- your_data

# End-to-end pipeline (multi-outcome scan + robust validation + pairwise)
out <- complete_function(df,
                         which = "all",       # "robust" or "binary" also allowed
                         n_boot = 100,
                         n_perm = 50,
                         ntree  = 500,
                         importance_method = "fixed",   # "neg_exp" or "net_clust" also allowed
                         plot = TRUE
)

```


### 2. PCA tools: `pca_scan_and_augment()`

In case you have many variables the `pca_scan_and_augment()` allows you to perform a PCA before running **astrolabe** to select the most important variables in the model. It then uses those variables to to causal inference on the variables that explain more.

```r
library(astrolabe)

# load data
df <- your_data

pca_results <- pca_scan_and_augment(df,
                                   ncomp = 5,
                                   top_k = 5,
                                   alpha = 0.05,
                                   verbose = TRUE,
                                   plot = TRUE,
                                   fixed_variables = NULL) 
```

In *fixed_variables* you can choose variables that you want to put in the model regardless the PCA results.

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

## 🏅 Extracting results

```r
out$res_all
```

---

## 📚 Citation & License

See the package manual for full references and function details.  
MIT License.
