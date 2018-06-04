library("jsonlite")
library("ggplot2")
library("readr")
setwd("~/Documents/Project/YelpAnalysis/")

business <- "dataset/business.json"
review <-read_lines(business, n_max = 200000)
business.df <- fromJSON(paste("[", paste(review, collapse = ","), "]"))
business.df

compare_m_gibbs <- function(y, ind, maxiter = 5000)
{
  
  ### weakly informative priors
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters, gamma
  eta0 <-1/2 ; t0 <- 50 ## tau_b hyperparameters, gamma
  mu0<-50 ## mean, norm
  gamma0 <- 1/25 
  ###
  
  ### starting values
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  an <- a0 + sum(n_m)/2
  ###
  
  ### setup MCMC
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  ###
  
  ### MCMC algorithm
  for(s in 1:maxiter) 
  {
    # sample new values of the thetas
    for(j in 1:m) 
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    bn <- b0 + ss/2
    tau_w <- rgamma(1, an, bn)
    
    #sample a new value of mu
    gammam <- m * tau_b + gamma0
    mum <- (mean(theta) * m * tau_b + mu0 * gamma0) / gammam
    mu <- rnorm(1, mum, 1/ sqrt(gammam)) 
    
    # sample a new value of tau_b
    etam <- eta0 + m/2
    tm <- t0 + sum((theta-mu)^2)/2
    tau_b <- rgamma(1, etam, tm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

business.df.filtered <- business.df[,c("stars", "neighborhood")]
business.df.filtered <- business.df.filtered[!(is.na(business.df.filtered$neighborhood) | business.df.filtered$neighborhood==""), ]
business.df.filtered

tapply(business.df.filtered$stars, business.df.filtered$neighborhood, mean)

fit <- compare_m_gibbs(business_df$stars, business_df$neighborhood)

