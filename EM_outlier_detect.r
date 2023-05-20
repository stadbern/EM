# Setting up test sets
X <- iris[,1:4]
outliers <- matrix(runif(n=8, min=1, max=7), nrow=2)
colnames(outliers) <- colnames(X)
X_out <- rbind(X,outliers)

## Code for Expectation-Maximization of multivariate Gaussian mixture models
## with uniform component




# code for k-means algorithm to set means and covariance matrices for EM alogrithme

km.maison <- function(X, k){
  p <- ncol(X)  # nombres de parametres
  n <- nrow(X)  # nombres d'observations
  
  Delta <- 1; iter <- 0; itermax <- 30 #Delta = mesure for convergence
  while(Delta > 1e-4 && iter <= itermax){
    # initiation, random centroid selection
    if(iter == 0){
      centroid <- X[sample(nrow(X), k),]
      centroid_mem <- centroid
    }
    
    # Euclidian distance between centroid and data points
    d <- sapply(1:k, function(c) sapply(1:n, 
                                        function(i) sum((centroid[c,] - X[i,])^2) ))
    # Centroid association
    cluster <- apply(d, 1, which.min)
    
    # New centroid settings
    centroid <- t(sapply(1:k, function(c) 
      apply(X[cluster == c,], 2, mean)))
    
    Delta <- sum((centroid - centroid_mem)^2)
    iter <- iter + 1; centroid_mem <- centroid
  }
  return(list(centroid = centroid, cluster = cluster))
}



# helper functions

mvnorm.cov.inv <- function(Sigma) {
  # Eigendecomposition of covariance matrix
  E <- eigen(Sigma)
  Lambda.inv <- diag(E$values^-1)   # diagonal matrix
  Q <- E$vectors
  return(Q %*% Lambda.inv %*% t(Q))
}
#multivariate Gaussian pdf (probability density function)
mvn.pdf.i <- function(xi, mu, Sigma)
  1/sqrt( (2*pi)^length(xi) * det(Sigma) ) * 
  exp(-(1/2) * t(xi - mu) %*% mvnorm.cov.inv(Sigma) 
      %*% (xi - mu) )

mvn.pdf <- function(X, mu, Sigma)
  apply(X, 1, function(xi) mvn.pdf.i(as.numeric(xi), as.numeric(mu), Sigma))

#multivariate Uniform pdf for outlier detection
mvu.pdf <- function(X) {
  mi <- apply(X, 2, min)
  ma <- apply(X, 2, max)
  uni_v <- 1/(ma-mi) 
  prod(uni_v)
}


#' EM multivarate Gaussians
gmm.maison.normal <- function(X, k){
  p <- ncol(X)  # number of parameters
  n <- nrow(X)  # number of observations
  
  Delta <- 1; iter <- 1; itermax <- 30; lh <- rep(NULL,itermax); bic <- rep(NULL,itermax)
  w.s <- array(dim=c(itermax,k)); cen.s <- array(dim=c(k,p,itermax))
  cov.s <- array(dim=c(p,p,k,itermax))
  while(Delta > 1e-4 && iter <= itermax){
    # initiation
    if(iter == 1){
      km.init <- km.maison(X, k)
      # centroids of k-means run
      mu <- km.init$centroid; mu_mem <- mu
      # pk start uniform
      w <- rep(1/k,k)
      cov <- array(dim = c(p, p, k))
      for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
        1/n * sum((X[km.init$cluster == c, i] - as.numeric(mu[c, i])) *
                    (X[km.init$cluster == c, j] - as.numeric(mu[c, j])))
    }
    
    # E-step
    mvn.c <- sapply(1:k, function(c) mvn.pdf(X, mu[c,], cov[,, c]))
    #print(mvn.c)
    r_ic <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
    if(sum(is.nan(r_ic)>=1)){
      break
    }
    #log-likelihood
    lh[iter] <- sum(log(rowSums(t(w*t(mvn.c)))))
    
    #BIC
    v <- k*(p*(p+1)/2 + p + 1) - 1
    bic[iter] <- sum(log(rowSums(t(w*t(mvn.c))))) - (v/2 * (log(n)))
    
    #Updates for output
    cen.s[,,iter] <- mu
    w.s[iter,] <- w
    cov.s[,,,iter] <- cov
    
    # M-step
    n_c <- colSums(r_ic)
    w <- n_c/sum(n_c)
    mu <- t(sapply(1:k, function(c) 1/n_c[c] * colSums(r_ic[, c] * X)))
    for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <-
      1/n_c[c] * sum(r_ic[, c] * (X[, i] - mu[c, i]) * r_ic[, c] *
                       (X[, j] - mu[c, j]))
    if(sum(is.nan(mu)>=1)){
      break
    }
    Delta <- sum((mu - mu_mem)^2)
    
    iter <- iter + 1; mu_mem <- mu
  }
  return(list(softcluster.norm = r_ic, BIC.norm = bic, likelihood.norm = lh,
              cluster.norm = apply(r_ic, 1, which.max), mu.s.norm = na.omit(cen.s),
              pk.s.norm = na.omit(w.s), cov.s.norm = na.omit(cov.s)))
}



#' EM multivariate Gaussian and uniform distribution

gmm.maison.uni <- function(X, k){
  p <- ncol(X)  
  n <- nrow(X)  
  uni_d <- mvu.pdf(X)
  
  Delta <- 1; iter <- 1; itermax <- 30; lh <- rep(NULL,itermax); bic <- rep(NULL,itermax)
  w.s <- array(dim=c(itermax,k+1)); cen.s <- array(dim=c(k,p,itermax))
  cov.s <- array(dim=c(p,p,k,itermax))
  while(Delta > 1e-4 && iter <= itermax){
    # initiation
    if(iter == 1){
      km.init <- km.maison(X, k)
      # centroids of k-means run
      mu <- km.init$centroid; mu_mem <- mu
      # pks start uniform
      w <- rep(1/(k+1),k+1)
      cov <- array(dim = c(p, p, k))
      for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
        1/n * sum((X[km.init$cluster == c, i] - as.numeric(mu[c, i])) *
                    (X[km.init$cluster == c, j] - as.numeric(mu[c, j])))
    }
    
    # E-step
    mvn.c <- sapply(1:k, function(c) mvn.pdf(X, mu[c,], cov[,, c]))
    #print(mvn.c)
    mvn.c <- cbind(mvn.c, rep(uni_d,length(X[,1])))
    
    r_ic <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
    #print(r_ic)
    if(sum(is.nan(r_ic)>=1)){
      break
    }
    #log-likelihood
    lh[iter] <- sum(log(rowSums(t(w*t(mvn.c)))))
    
    #BIC
    v <- k*(p*(p+1)/2 + p + 1) - 1 + p # +p à la fin uniform
    bic[iter] <- sum(log(rowSums(t(w*t(mvn.c))))) - (v/2 * (log(n)))
    
    #Updates for output
    cen.s[,,iter] <- mu
    w.s[iter,] <- w
    cov.s[,,,iter] <- cov
    # M-step
    n_c <- colSums(r_ic)
    w <- n_c/sum(n_c)
    #print(w)
    mu <- t(sapply(1:k, function(c) 1/n_c[c] * colSums(r_ic[, c] * X)))
    for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <-
      1/n_c[c] * sum(r_ic[, c] * (X[, i] - mu[c, i]) * r_ic[, c] *
                       (X[, j] - mu[c, j]))
    if(sum(is.nan(mu)>=1)){
      break
    }
    Delta <- sum((mu - mu_mem)^2)
    mu_mem <- mu
    
    iter <- iter + 1 
  }
  return(list(softcluster.uni = r_ic, BIC.uni = bic, likelihood.uni = lh,
              cluster.uni = apply(r_ic, 1, which.max), mu.s.uni = na.omit(cen.s),
              pk.s.uni = na.omit(w.s), cov.s.uni = na.omit(cov.s)))
}



# run functions
get.norm.uni <- function(X, k_clusters, n_runs=3){
  x.norm <- lapply(1:n_runs, function(i){
    gmm.maison.normal(X,k_clusters)
    
  }) 
  x.uni <- lapply(1:n_runs, function(i){
    
    gmm.maison.uni(X,k_clusters)
    
  }) 
  
  return (list(x.norm = x.norm, x.uni = x.uni))
  
}

# result functios
get.all <- function(X, k_clusters, n_runs=3){
  cat("------------------------------------------------------\n")
  cat("------------------------------------------------------\n")
  cat("Model selection for k Gaussians and possible uniform\n")
  cat(" distribution for outlier detection\n")
  a <- get.norm.uni(X,k_clusters, n_runs)
  n.u <- rep(0,2)
  BICs.n <- unlist(lapply(1:n_runs, function(i){tail(a$x.norm[[i]]$BIC.norm,1)}))
  BICs.u <- unlist(lapply(1:n_runs, function(i){tail(a$x.uni[[i]]$BIC.uni,1)}))
  n.u[1] <- BICs.n[which.max(BICs.n < 0)]
  n.u[2] <- BICs.u[which.max(BICs.u < 0)]
  select.m <- which.max(n.u)
  
  if(select.m == 1){
    m <- which.max(BICs.n < 0)
    
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("       Model with uniform distribution selected\n")
    cat("               Detection of outliers: No\n")
    cat("------------------------------------------------------\n")
    cat("         BIC value with uniform component: ",n.u[2],"\n")
    cat("         BIC value without uniform component: ",n.u[1],"\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("                    Model parameters\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("\n")
    pks.norm <- tail(a$x.norm[[m]]$pk.s.norm,1)
    cat("              pks for the different components\n")
    print(pks.norm)
    cat("------------------------------------------------------\n")
    cat("\n")
    mus.norm <- a$x.norm[[m]]$mu.s.norm[,,length(a$x.norm[[m]]$likelihood.norm)]
    colnames(mus.norm) <- colnames(X)
    cat("       µs for the different Gaussian components\n")
    cat("------------------------------------------------------\n")
    print(mus.norm)
    cat("------------------------------------------------------\n")
    cat("\n")
    covs.norm <- a$x.norm[[m]]$cov.s.norm[,,,length(a$x.norm[[m]]$likelihood.norm)]
    #colnames(mus.uni) <- colnames(X)
    cat("         Sigmas for the different Gaussian components\n")
    cat("------------------------------------------------------\n")
    print(covs.norm)
    cat("\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("                    Data\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("Probabilities for different components (1-n)\n")
    cat("       (last column = uniform component)\n")
    cat("------------------------------------------------------\n")
    print(a$x.norm[[m]]$softcluster.norm)
    cat("\n")
    cat("------------------------------------------------------\n")
    cat(" Dataframe and component allocation (last column)\n")
    cat("------------------------------------------------------\n")
    cat("\n")
    c <- cbind(X,a$x.norm[[m]]$cluster.norm)
    colnames(c)[length(colnames(c))] <- "Cluster"
    print(c)
    cat("\n")
    cat("------------------------------------------------------\n")
    cat("\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("             Plots ===================>>>>>>\n")
    cat("------------------------------------------------------\n")
    cat("Plots: 1: Scatterplots showing pairs of variables\n")
    cat("       2: Evolution of likelihood over iterations\n")
    cat("------------------------------------------------------\n")
    pairs(X_out, lower.panel = NULL, col = a$x.norm[[m]]$cluster.norm,
          main = "Pairs of variables/Colour = cluster")
    plot(x = seq(1,length(a$x.norm[[m]]$likelihood.norm)), y = a$x.norm[[m]]$likelihood.norm, type="l",
         main = "Evolution of likelihood", xlab = "n iteration", ylab = "log-likelihood")
    return(list(pks.norm=pks.norm, mus.norm=mus.norm, covs.norm=covs.norm))
    
  }
  if(select.m == 2){
    m <- which.max(BICs.u < 0)
    
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("        Model with uniform distribution selected")
    cat("               Detection of outliers: Yes\n")
    cat("------------------------------------------------------\n")
    cat("         BIC value with uniform component: ",n.u[2],"\n")
    cat("         BIC value without uniform component: ",n.u[1],"\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("                    Model parameters\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    pks.uni <- tail(a$x.uni[[m]]$pk.s.uni,1)
    cat("              pks for the different components\n")
    print(pks.uni)
    cat("\n")
    cat("------------------------------------------------------\n")
    mus.uni <- a$x.uni[[m]]$mu.s.uni[,,length(a$x.uni[[m]]$likelihood.uni)]
    colnames(mus.uni) <- colnames(X)
    cat("       µs for the different Gaussian components\n")
    cat("------------------------------------------------------\n")
    print(mus.uni)
    cat("------------------------------------------------------\n")
    cat("\n")
    covs.uni <- a$x.uni[[m]]$cov.s.uni[,,,length(a$x.uni[[m]]$likelihood.uni)]
    #colnames(mus.uni) <- colnames(X)
    cat("         Sigmas for the different Gaussian components\n")
    cat("------------------------------------------------------\n")
    print(covs.uni)
    cat("\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("                    Data\n")
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("Probabilities for different components (1-n)\n")
    cat("       (last column = uniform component)\n")
    cat("------------------------------------------------------\n")
    print(a$x.uni[[m]]$softcluster.uni)
    cat("\n")
    cat("------------------------------------------------------\n")
    cat(" Dataframe and component allocation (last column)\n")
    cat("------------------------------------------------------\n")
    cat("\n")
    c <- cbind(X,a$x.uni[[m]]$cluster.uni)
    colnames(c)[length(colnames(c))] <- "Cluster"
    print(c)
    cat("\n")
    
    cat("------------------------------------------------------\n")
    cat("------------------------------------------------------\n")
    cat("             Plots ===================>>>>>>\n")
    cat("------------------------------------------------------\n")
    cat("Plots: 1: Scatterplots showing pairs of variables\n")
    cat("       2: Evolution of likelihood over iterations\n")
    cat("------------------------------------------------------\n")
    pairs(X_out, lower.panel = NULL, col = a$x.uni[[m]]$cluster.uni,
          main = "Pairs of variables/Colour = cluster")
    plot(x = seq(1,length(a$x.uni[[m]]$likelihood.uni)), y = a$x.uni[[m]]$likelihood.uni, type="l",
         main = "Evolution of likelihood", xlab = "n iteration", ylab = "log-likelihood")
    return(list(pks.uni=pks.uni, mus.uni=mus.uni, covs.uni=covs.uni))
  }
  
}

# Enter dataframe, number of ks and number of runs
r <- get.all(X_out, 3, n_runs=3)




