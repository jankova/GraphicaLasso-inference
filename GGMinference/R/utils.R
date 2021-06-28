#

nodewise <- function(x, lambdaj){
  p <- ncol(x)
  n <- nrow(x)
  C <- matrix(rep(0,p*p),nrow=p)
  for(j in 1:p){
    xj <- x[,j]
    xmj <- x[,-j]
    gamaj <- RPtests::sqrt_lasso(xmj, xj)
    tau2j <- as.numeric(t(x[,j] - x[,-j]%*%gamaj) %*% (x[,j] - x[,-j]%*%gamaj)/n)
    tau.tilde2j <- tau2j + lambdaj*sqrt(tau2j)*sum(abs(gamaj))
    Cj <- gamaj/tau.tilde2j;
    if(j==1){ first <- numeric(0);}else first <- Cj[1:(j-1)];
    if(j==p){ second <- numeric(0);} else second <- Cj[j:(p-1)];
    C[j,] <- c(-first,1/tau2j,-second);
  }
  C
}


run.checks <- function(X, standardize, alpha, theta0, pmethod){
  if(!is.logical(standardize))
    stop("standardize should be a logical (TRUE or FALSE)")

  if (!is.matrix(X)){
    stop("X should be a matrix with at least one column.")
  }


  if(!is.numeric(alpha)){
    stop("alpha should be numeric")
  }

  if(alpha < 0 || alpha > 1 )
    stop("alpha is significance level and should be between 0 and 1")

  if(!is.character(pmethod))
    stop("pmethod should be a string")

  if(!pmethod %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                     "fdr", "none")){
    stop("pmethod should be one of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH',
         'BY','fdr', 'none'")
  }

  if (!is.matrix(theta0)){
    stop("theta0 should be a matrix with at least one column.")
  }
}

inference <- function(X, alpha, method, rho0, theta0, pmethod, visual){

  np <- dim(X)
  if (is.null(np) | (np[2] < 1L)){
    stop("X should be a matrix with at least one column.")
  }
  n <- as.integer(np[1])
  p <- as.integer(np[2])

  S <- var(X)

  u <- qnorm(1 - alpha/2)

  if(method == "glasso"){
    theta.hat <- glasso(S, penalize.diagonal = FALSE, rho = rho0)$wi
  }else if(method == "neighbourhood"){
    theta.hat <- nodewise(X, sqrt(log(p)/n))
  }else if(method == "precmat"){
    theta.hat <- solve(S)
  }

  theta.debiased <- theta.hat + t(theta.hat) - theta.hat %*% S %*% t(theta.hat)

  # estimate of standard deviation
  sigma.hat <- sqrt(theta.hat^2 + diag(theta.hat) %*% t(diag(theta.hat)))

  # confidence intervals
  CI.lower.end <- theta.debiased - u * sigma.hat / sqrt(n)
  CI.upper.end <- theta.debiased + u * sigma.hat / sqrt(n)

  test.stat <- sqrt(n) * (theta.debiased - theta0) / sigma.hat
  p.values <- 1 - pnorm( abs(test.stat) )
  p.values.adjusted <- matrix(p.adjust(p.values, method = pmethod), ncol = p)
  selected <- which(p.values.adjusted < 0.05, arr.ind = TRUE)
  selected <- selected[apply(selected, 1, function(x) x[1] != x[2]), ]

  if(visual == TRUE){
    heatmap(p.values.adjusted, revC = TRUE,
            Rowv = NA, Colv = NA, keep.dendro = FALSE)
  }

  conf.ints <- list(CI.lower.end, CI.upper.end)

  list(conf.ints = conf.ints, p.values.adjusted = p.values.adjusted,
       p.values = p.values, selected = selected)
}
