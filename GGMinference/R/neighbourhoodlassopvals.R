#' P-values and confidence intervals based on the nodewise (neighbourhood) lasso.
#'
#' @param X Input matrix with \code{n} rows and \code{p} columns, each row an observation vector.
#' @param standardize Should design matrix be standardized to unit column standard deviation.
#' @param alpha Significance level for testing and confidence intervals.
#' @param pmethod Method to use for multiple testing correction ("holm", "hochberg", "hommel", 'bonferroni', 'BH',
#' "BY","fdr" or "none"). The default is "BH".
#' @param theta0 (Optional) A \code{pxp} matrix for custom testing of hypothesis theta_ij = theta0_ij.
#' @param visual Logical indicating whether a heatmap of adjusted p-values should be plotted. Default is False.
#' @description
#' P-values and confidence intervals for edge weights in high-dimensional
#' undirected Gaussian graphical models.
#' The method is based on the nodewise (neighbourhood) lasso as the initial estimator. For details on the methodology, see
#' Janková, J. and van de Geer, S. (2017)
#' \emph{Honest confidence regions and optimality in high-dimensional precision matrix estimation}.
#' @return Returns a list containing p-values, adjusted p-values, confidence intervals
#' and significant variables as specified below.
#' \item{pvals}{\code{pxp} matrix containing individual p-values for each parameter.}
#' \item{pvals.adjusted}{\code{pxp} matrix containing p-values for each parameter, corrected with a multiple-testing adjustment specified by \code{pmethod}.}
#' \item{conf.ints}{A list with two \code{pxp} matrices: the matrices hold upper and lower bounds of individual confidence intervals respectively (without
#' multiple testing adjustment).}
#' \item{selected}{A matrix with 2 columns where rows are pairs of coordinates of edges that are selected as significant at level \code{alpha} (with multiple testing adjustment as specified by
#' \code{pmethod}).}
#' @references Janková, J. and van de Geer, S. (2018)
#' \emph{Inference in high-dimensional undirected graphical models}
#' \url{https://arxiv.org/abs/1801.08512}
#' @examples
#' # Inference for edge weights using data generated from a Gaussian graphical model
#'
#' library(MASS)
#' set.seed(1)
#'
#' p <- 100
#' n <- 150
#' rho <-  0.3
#' Theta <- diag(p) + cbind(rho*diag(p)[,-1], rep(0,p)) + t(cbind(rho*diag(p)[,-1], rep(0,p)))
#'
#' X <- mvrnorm(n, rep(0,p), Sigma = solve(Theta))
#'
#' nlasso.inference <- nodewiselasso.pvals(X, standardize = FALSE, alpha = 0.05, visual = FALSE)
#' nlasso.inference$p.values.adjusted
#'
#'
#' @export
#' @import stats
#' @import MASS
#' @import RPtests
#' @import glmnet
#' @import glasso
#'


nodewiselasso.pvals <- function(X,
                            standardize = FALSE,
                            alpha = 0.05,
                            theta0 = matrix(0, ncol(X), ncol(X)),
                            pmethod = "BH",
                            visual = FALSE){

  run.checks(X, standardize, alpha, theta0, pmethod)

  if(standardize == TRUE)
    X <- scale(X)

  inference(X, alpha, method = "neighbourhood", rho0 = NULL, theta0, pmethod, visual)

}
