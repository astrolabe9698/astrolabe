

**astrolabe** is an R package that implements a novel methodology for inferring **non-linear causal relationships** among variables.  
The approach combines machine learning with information-theoretic principles, using entropy of residuals to evaluate causal directions.

## 🔬 Methodology

The method relies on the entropy of residuals from models fitted with **random forests**.  
We define a score for each model:

\[
H(m) = e^{-S_{\text{knn}}(m)} ,
\]

where \( S_{\text{knn}}(m) \) is a **Kozachenko–Leonenko entropy estimator** applied to the joint space of predictors and residuals.  
This formulation is robust for continuous residuals and avoids sensitivity to binning choices, as entropy is estimated from the geometry of nearest-neighbour distances.

The procedure consists of:
1. Fitting random forest models for candidate causal directions.  
2. Computing the residual entropy score \( H \).  
3. Performing **bootstrap–permutation testing** to assess significance.  
4. Selecting relevant predictors via **variable importance measures**.  

This provides a framework to detect non-linear causality in high-dimensional datasets where traditional linear approaches may fail.

## 📦 Installation

### Development version (from GitHub)
```r
# install.packages("remotes") # if not already installed
remotes::install_github("your-username/astrolabe")
