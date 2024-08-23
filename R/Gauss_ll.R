#Coded by Ricardo Carrizo V
#Version at 22-08-2024
library(TruncatedNormal)
library(mnormt)

Gauss_ll <- function( x , m , C , 
                      Reg = T , eps = 1e-6 , 
                      Continuity.corr = F){
  #Function to compute the multi-dimensional Gaussian log-likelihood of a data vector.
  
  #Input:
  #
  #     x <- data vector
  #     m <- mean vector of the Gaussian vector
  #     C <- covariance matrix of the Gaussian vector (must be non-negative definite)
  #     Reg <- Boolean. Set to TRUE if a regularization of the covariance is desired.
  #     eps <- If Reg = TRUE, positive value to be added to the variances of the components.
  #     Continuity.corr <-  Boolean. If TRUE, a continuity-corrected version of the log-likelihood is computed.
  #
  
  #Output:
  #     ll <- If Continuity.corr = F, the Gaussian log-likelihood of the data x.
  #                         The mnormt Package is used.
  #           
  #           If Continuity.corr = T, the log-probability of the Gaussian vector to belong to [x-0.5 , x + 0.5] (useful when x is integer and the distribution is approximatively Gaussian).
  #                         The TruncatedNormal Package is used to compute the continuity-corrected Gaussian probability.
  
  if(Reg){
    if(length(x) > 1){
      diag(C) <- diag(C) + eps
    }else{
      C <- C + eps
    }
    if(Continuity.corr){
      return(log(TruncatedNormal::pmvnorm( mu = m , 
                                           sigma = C ,
                                           lb = x - 0.5 , 
                                           ub = x + 0.5 ,
                                           B = 1e6) )  )
    }
    return( mnormt::dmnorm( x , mean = m , varcov = C , log = T )  )
  }
  
  if(Continuity.corr){
    return(log(TruncatedNormal::pmvnorm( mu = m , 
                                         sigma = C ,
                                         lb = x - 0.5 , 
                                         ub = x + 0.5 ,
                                         B = 1e6) )  )
  }
  
  eig.C <- eigen(C)
  r <- sum( eig.C$values > 0 )
  if(r == 0){
    if( isTRUE(all.equal(x,m)) ){
      return(Inf)
    } 
    return(-Inf)
  }
  y <- t(eig.C$vectors)%*%(x-m)
  if(r == length(x)){
    return( sum(dnorm(  y , sd = sqrt( eig.C$values ) , log = T ) ) )
  }
  if( isTRUE( all.equal(  0 ,  y[ (r+1):length(x) ]   ) ) ){
    return( sum(dnorm( y[1:r] , sd = sqrt(eig.C$values[1:r]), log = T ) ))
  }
  return(-Inf)
}


# #Example with functions from this repository (must be defined before):
# N <- 1000
# t0 <- 0
# tL <- 0.5
# tH <- 2
# dt <- 1/60
# d <- 2
# Dc.B.L = rbind(  c(0,0) , c(1,0) , c(-1,1) , c(-1 , -1) , c(0,1))
# dx <- 1
# theta.X.free <- c(3,-1 , 1)
# alpha <- 0.1
# 
# Q <- Sim_Capture( N = N ,
#                    t0 = t0 , tL = tL , tH = tH , dt = dt ,
#                    d = d ,
#                    Dc.B.L = Dc.B.L , dx = dx ,
#                    alpha = alpha  ,
#                    theta.X.free = theta.X.free  )
# 
# Q <- as.vector(Q)
# 
# mom.Q <- Moments_Capture( N = N ,
#                           t0 = t0 , tL = tL , tH = tH , dt = dt ,
#                           d = d ,
#                           Dc.B.L = Dc.B.L , dx = dx ,
#                           alpha = alpha  ,
#                           theta.X.free = theta.X.free  )
# 
# Gauss_ll( x = Q , 
#           m = mom.Q$mean.vec , 
#           C = mom.Q$cov.mat , Reg = F )
# 
# Gauss_ll( x = Q , 
#           m = mom.Q$mean.vec , 
#           C = mom.Q$cov.mat , Reg = T )
# 
# Gauss_ll( x = Q , 
#           m = mom.Q$mean.vec , 
#           C = mom.Q$cov.mat , Reg = F , Continuity.corr = T )
