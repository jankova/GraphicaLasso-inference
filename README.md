# Inference in high-dimensional graphical models with the Graphical Lasso
The package provides confidence intervals for edges in high-dimensional undirected Gaussian graphical models and p-values for testing their significance. For the high-dimensional case, confidence intervals and tests are based on a de-biased version of the graphical lasso [3] or the nodewise (neighbourhood) lasso [1]. For the low-dimensional case, confidence intervals can be computed using the precision matrix estimator.

## Installation
```
install.packages("devtools")

library(devtools)

install_github("jankova/GGMinference_R_package/GGMinference")

library(GGMinference)
```

<!---## Graphical models

## Methods--->


## Examples
The following code demonstrates the usage on a dataset generated from a Gaussian graphical model with a tri-diagonal precision matrix.

```
# Inference for edge weights using data generated from a Gaussian graphical model
library(MASS) 
set.seed(1)

p <- 100
n <- 150
rho <- 0.3

Theta <- diag(p) + cbind(rho*diag(p)[,-1], rep(0,p)) + t(cbind(rho*diag(p)[,-1], rep(0,p)))
X <- mvrnorm(n, rep(0,p), Sigma = solve(Theta))

# p-values using the de-biased graphical lasso
glasso.inference <- glasso.pvals(X, standardize = FALSE, alpha = 0.05, rho0 = sqrt(log(p)/n), pmethod = "BH", visual = FALSE)
glasso.inference$p.values.adjusted

# p-values using the de-biased nodewise lasso
nodelasso.inference <- nodelasso.pvals(X, standardize = FALSE, alpha = 0.05, visual = FALSE) 
nodelasso.inference$p.values.adjusted

# p-values using the precision matrix estimator
precmat.inference <- precmat.pvals(X, standardize = FALSE, alpha = 0.05, pmethod = "BH", visual = FALSE)
precmat.inference$p.values.adjusted

```

## References

[1] Janková, J. and van de Geer, S. (2018) Inference in high-dimensional undirected graphical models, Handbook of Graphical Models, CRC Press.

[2] Janková, J. and van de Geer, S. (2017) Honest confidence regions and optimality in high-dimensional precision matrix estimation, Test 26 (1), 143-162.

[3] Janková, J. and van de Geer, S. (2015) Confidence intervals for high-dimensional inverse covariance estimation, Electronic Journal of Statistics 9 (1), 1205-1229. 

