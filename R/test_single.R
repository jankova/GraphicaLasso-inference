# functions:
# estimate of sqrt prec mat - use rajens code??
# F1: CIS - choose nw/glasso,
#  - all pars / bonf correction optional
# testing:
#-F2: single -> nw -okay-two options, glasso - not clear
#-F3: groups -> MB
# F4: varsel -> padjust
#' Confidence intervals and tests for significance of edges in high-dimensional sparse
#' undirected graphical models
#'
#' Can test the significance of individual predictors or (potentially large) groups of predictors
#' in a low- and high-dimensional undirected graphical models.
#' Outputs confidence interval(s) or a p-value(s) for the corresponding test.
#'
#' @param X Input matrix with \code{nobs} rows, each an observation vector.
#' @param y Response vector.
#' @param alpha Level of the confidence intervals.
#' @param method nodewise (default) or glasso.
#' @param y Response vector.
#' @param y Response vector.
#'
#' @param B If \code{test} is \code{group}, B is the number of bootstrap samples. Note that the p-value for grouptest will always be
#'   at least 1/(B+1).
#' @param usecv If \code{test} is \code{nonlin}, B is the number of bootstrap samples. Note that the p-value for grouptest will always be
#'   at least 1/(B+1).
#' @details
#'   When \code{test} is \code{nonlin}, the function tests nonlinearity in the
#'   condition mean. ?Should we explain the method in detail?
#'   When \code{test} is \code{group}, the function tests significance of a set of variables
#'   whos indices are specified by G.
#' @return The output is a matrix .
#' @references Janková, J. and van de Geer, S. (2015,2017,2018)
#' \emph{Goodness-of-fit testing in high-dimensional generalized linear models}
#' \url{https://arxiv.org/abs/1908.03606}
#' @seealso \code{\link{sqrt_lasso}}
#' @examples
#' # Testing for nonlinearity
#' set.seed(1)
#' X <- scale(matrix(runif(300*50), 300, 50))
#' z <- X[, 1] + X[, 1]^4 + rnorm(nrow(X))
#' pr <- 1/(1 + exp(-z))
#' y <- rbinom(n, 1, pr)
#' (out <- GRPtest(X, y, test = "nonlin", family = "binomial", usecv = TRUE, lam = 0.1))
#'
#' # Testing significance of a group
#' z <- X[, 1:5] %*% rep(1, 5) + rnorm(nrow(X))
#' pr <- 1/(1 + exp(-z))
#' y <- rbinom(n, 1, pr)
#' (out <- GRPtest(X, y, test = "group", family = "binomial", G = 5:10, usecv = TRUE, lam = 0.1, B = 1000))
#'
#' @export
#' @useDynLib
#' @importFrom Rcpp sourceCpp
#' @import stats
#' @import MASS
#' @import glasso
#' @import RPtests
#'

test_single <- function(X, alpha, method = "nodewise", rho0 = sqrt(log(p)/n), theta0 = matrix(0, p, p)){

  u <- qnorm(1 - alpha/2)
  n <- nrow(X)
  p <- ncol(X)
  S <- var(X)

  if(method == "glasso"){
    theta.hat <- glasso(S, penalize.diagonal = FALSE, rho = rho0)$wi
  }else if(method == "nodewise"){
    theta.hat <- nodewise(X)
  }
  theta.debiased <- theta.hat + t(theta.hat) - theta.hat %*% S %*% t(theta.hat)

  # estimate of standard deviation
  sigma.hat <- sqrt(theta.hat^2 + diag(theta.hat) %*% t(diag(theta.hat)))

  # confidence intervals
  CI.lower.end <- theta.debiased - u * sigma.hat / sqrt(n)
  CI.upper.end <- theta.debiased + u * sigma.hat / sqrt(n)

  test.stat <- sqrt(n) * (theta.debiased - theta0) / sigma.hat
  p.values <- 1 - pnorm( abs(test.stat) )
  p.values <- p.adjust(p.values)
  which(p.values < 0.05, arr.ind = TRUE)
  selected <- which(p.values < 0.05)
  heatmap(theta.hat)

  list(CI.lower.end = CI.lower.end, CI.upper.end = CI.upper.end, debiased.estimator = theta.debiased)

}
