########################################
# Gaussian process
########################################

# Element-wise matrix autocorrelation function
cov.m <- function(rows, k, d) {
  params <- k$params
  k <- k$k
  # Covariance of the elements from matrix d
  k(d[rows[1],], d[rows[2],], params)
}

# Matrix autocorrelation
K.m <- function(X, k, symmetric = F, diagonal=1) {
  # Get row indexes in the right order to fill the lower triangle
  # of the covariance matrix
  rows <- combn(1:nrow(X), 2)
  # get the covariances of the rows
  covs <- as.matrix(apply(rows, 2, cov.m, k = k, d = X))
  
  Ks <- list()
  
  # If the parameter is a vector (e.g. the lengthscale parameter in an ARD kernel)
  if(ncol(covs) > 1) {
    for(i in 1:nrow(covs)) {
      K <- diag(diagonal, nrow(X))
      K[lower.tri(K)] <- covs[i,]
      K <- t(K)
      
      if (symmetric)
        K[lower.tri(K)] <- t(K)[lower.tri(K)]
      
      Ks[[i]] <- K
    }
    Ks
  } else { #if the parameter is a scalar
    K <- diag(diagonal, nrow(X))
    K[lower.tri(K)] <- covs
    K <- t(K)
    
    if (symmetric)
      K[lower.tri(K)] <- t(K)[lower.tri(K)]
    K
  }
}

# Element-wise Matrix-Matrix covariance
# rows - set of indexes
cov.mm <- function(rows, k, A, B) {
  params <- k$params
  k <- k$k
  
  k(A[rows[1],],
    B[rows[2],], params)
}

# Matrix-Matrix covariance
K.mm <- function(X, Y, k) {
  rows <- expand.grid(1:nrow(X), 1:nrow(Y))
  covs <- apply(
    X = rows,
    MARGIN = 1,
    FUN = cov.mm,
    A = X,
    B = Y,
    k = k
  )
  matrix(covs, nrow = nrow(X), ncol = nrow(Y))
}

# trace of a matrix
trace <- function(M) {
  sum(diag(M))
}

# Partial derivative of the marginal likelihood
gp.marginal.pd <- function(gp, k.pd) {
  # Compute the matrix of element-wise partial derivatives of the kernel function
  d_K <- K.m(gp$X, 
             k = list(k=k.pd, params=gp$kernel$params), 
             symmetric = T, 
             diagonal = 0)
  
  pds <- c()
  if(is.list(d_K)) {
    for (i in 1:length(d_K)) {
      pd_i <- (1 / 2) * trace((gp$alpha %*% t(gp$alpha) - chol2inv(gp$L)) %*% d_K[[i]])
      pds <- c(pds, pd_i)
    }
    pds
  } else {
    (1 / 2) * trace((gp$alpha %*% t(gp$alpha) - chol2inv(gp$L)) %*% d_K)
  }
}

# Both the marginal likelihood and gradient functions have the form: function(param_values, gp)
# This is because the optim function requires this format. param_values is a vector of kernel that is converted
# back to kernel parameters by the gp.kernel.refit function

# Marginal Likelihood
gp.marginal <- function(param_values, gp) {
  # TODO: check if params are different
  gp <- gp.kernel.reset(param_values, gp)
  gp <- gp.refit(gp)
  
  # marginal likelihood
  D <- gp$det_Ky
  n <- nrow(gp$X)
  ml <- -1 / 2 * t(gp$y) %*% gp$alpha - D - n/2 * log(2 * pi)

  if(gp$verbose){
    print('--------optimizer params------------')
    print(param_values)
    print('--------marginal likelihood---------')
    print(ml)
  }
  
  ml
}

gp.gradient <- function(param_values, gp) {
  gp <- gp.kernel.reset(param_values, gp)
  gp <- gp.refit(gp)
  
  gradient <- c()
  # Append each partial derivative (pd) to the gradient
  for (param in names(gp$kernel$params)) {
    pd <- gp$kernel$pds[[param]]
    g <- gp.marginal.pd(gp, pd)
    gradient <- c(gradient, g)
  }
  
  if(gp$verbose){
    print('-----------gradient-----------')
    print(gradient)
    print(length(gradient))
    print('------------------------------')
  }
  
  gradient
}

gp.optimize <- function(gp, lower, upper) {
  init_params = unlist(gp$kernel$params)
  
  optimized <- optim(
    par = init_params,
    fn = gp.marginal,
    gr = gp.gradient,
    method = 'L-BFGS-B',
    lower=lower, # upper bound of parameters for optimizer
    upper=upper, # lower bound
    control = list(fnscale = -1), #maximize marginal likelihood, optim minimizes by default
    gp=gp
  )
  
  best_params <- optimized$par
  
  if(gp$verbose) {
    print('-------optimized---------')
    print(optimized)
  }
  
  gp <- gp.kernel.reset(best_params, gp)
  
  gp <- gp.refit(gp)
  
  gp
}

# Resets the kernel parameters before re-fitting the model
gp.kernel.reset <- function(params, gp) {
  offset <- 0
  for (i in 1:length(gp$kernel$params)) {
    # Handle the special case of vector kernel parameters.
    # Usually these vector parameters are the length scales in ARD kernels
    if(length(gp$kernel$params[[i]]) > 1) {
      offset <- offset + length(gp$kernel$params[[i]]) - 1
      gp$kernel$params[[i]] <- params[i:(i+offset)]
    } else {
      gp$kernel$params[[i]] <- params[(i+offset)]
    }
  }
  gp
}

# Fit a gaussian process using the supplied kernel
gp.fit <- function(X, y, kernel, sigma_n, verbose=F) {
  # Covariance matrix
  K <- K.m(X, kernel) + diag(sigma_n, nrow(X))
  # Cholesky factors of covariance matrix
  L <- chol(K)
  
  
  gp <- list()
  gp$X <- X
  gp$y <- y
  gp$kernel <- kernel
  gp$alpha <- chol2inv(L) %*% y
  gp$L <- L
  gp$det_Ky <- det(K.m(gp$y, kernel))
  
  gp$sigma_n <- sigma_n
  
  gp$verbose = verbose
  
  gp
}

# Simplified 
gp.refit <- function(gp) {
  K <- K.m(gp$X, gp$kernel) + diag(gp$sigma_n, nrow(gp$X))
  L <- chol(K)
  
  gp$L <- L
  gp$alpha <- chol2inv(L) %*% gp$y
  gp$det_Ky <- det(K.m(gp$y, gp$kernel))
  gp
}

gp.predict <- function(gp, x_test) {
  covs <- K.mm(gp$X, x_test, gp$kernel)
  f_mean <- t(covs) %*% gp$alpha
  vars <- K.m(x_test, gp$kernel) - t(covs) %*% chol2inv(gp$L) %*% covs
  list(means = f_mean, vars = vars)
}

#Kernel functions
# each function should be of the form function(x, y, params) 
# params is a list of parameters for the function

# Radial Basis Function (RBF) kernel
rbf <- function(x, y, params) {
  sigma <- params$sigma
  exp(-1/2*sum(sigma^-2*(x - y)^2))
}

# List of partial derivative functions for each parameter
rbf.pds <- list(
  sigma = function(x, y, params){
    sigma <- params$sigma
    exp(-1/2*sum(sigma^-2*(x - y)^2)) * (sigma^(-3)*(x-y)^2)
  }
)

# List holding the covariance function, partial derivatives and parameters
rbf.k <- function(sigma) {
  list(
    k = rbf,
    pds = rbf.pds,
    params = list(sigma=sigma)
  )
}

# RBF kernel with Automatic Relevence Determination 
# The lengthscale parameter (l) is a vector 
rbf.ard <- function(x, y, params) {
  l <- unlist(params$l)
  sigma <- params$sigma
  sigma^2 * exp(-1/2*sum(l^-2*(x - y)^2))
}

rbf.ard.pds <-  list(
  l = function(x, y, params) {
    l <- unlist(params$l)
    sigma <- params$sigma
    sigma^2 * exp(-1/2*sum(l^-2*(x - y)^2)) * (l^-3*(x-y)^2)
  },
  sigma = function(x, y, params) {
    l <- unlist(params$l)
    sigma <- params$sigma
    2*sigma * exp(-1/2*sum(l^-2*(x - y)^2))
  }
)

rbf.ard.k <- function(sigma, l) {
  list(
    k = rbf.ard,
    pds = rbf.ard.pds,
    params = list(sigma=sigma, l=l)
  )
}

# Example of fitting a complex 1 dimensional function
gp.example.1d <- function() {
  missing <- c(seq(.1, 5, .1), seq(8.5, 10, .1))
  full <- seq(.1, 10, .1)
  
  sparse <- seq(.1, 20, 3)
  # 1-d example
  X <- matrix(full)
  y <- matrix(log(abs(cos(X))) ^ 2 + cos(X) ^ 2 + sin(log(X)) + sin(X) + log(abs(cos(X) - cos(X) * log(abs(cosh(X))))))
  # y <- matrix(log(cos(X)^4) ^ 2 + cos(X) ^ 3 + runif(nrow(X), -.1, .1))
  
  gp <- gp.fit(X, y, sigma_n = 0.2, k = rbf.k(sigma=.05))
  gp <- gp.optimize(gp,
                    lower=c(-Inf),
                    upper=c(Inf))
  print(gp$kernel$params)
  
  x_test <- matrix(seq(.1, max(X), .1))
  p <- gp.predict(gp, x_test)
  
  plot(X, y)
  lines(x_test, p$means, pch = 'x', col = 'red')
  points(x_test,
         p$means + diag(p$vars),
         pch = '-',
         col = 'blue')
  points(x_test,
         p$means - diag(p$vars),
         pch = '-',
         col = 'green')
}

# 2-d example
# Requires the rgl library
# This takes a long time to run
# Runtime can be reduced by changing the full variable
gp.example.2d <- function() {
  library('rgl')
  full <- seq(.1, 10, 1)+runif(10,-.01,.01)
  X <- data.matrix(expand.grid(full, full), rownames = F)
  y <- sin(X[, 1]*pi^2) + cos(-X[,2])*log(X[,1]^2)
  y <- matrix(y + runif(nrow(X),-.5, .5))
  persp3d(full, full, y, color = "green", alpha = .5)
  
  gp <- gp.fit(X, y, 
               sigma_n = 0.2, 
               k = rbf.ard.k(sigma = .5, l = c(1,1)))
  
  gp <- gp.optimize(gp,
                    lower=c(0,-Inf,-Inf),
                    upper=c(1,Inf,Inf))
  
  print('------final params-----------')
  print(gp$kernel$params)
  
  x_test <- X
  
  p <- gp.predict(gp, x_test)
  
  points3d(X[, 1], X[, 2], p$means, col = 'red')
}

gp.example.1d()
# gp.example.2d()

# TODO: simplifying code for handling scalar vs vector partial derivatives
# TODO: matern 3/2, matern 5/2, periodic kernel, etc.
# TODO: 
