#Coded by Ricardo Carrizo V
#Version at 22-08-2024
library(mnormt)

Moments_Snapshot <- function( N ,
                              t0  , t  , 
                              d  , 
                              Regions.B.L , dx  , 
                              p , 
                              Traj = "Brownian+Advection", theta.X){
  #Function to obtain the mean vector and covariance matrix of a Snapshot model random matrix (wrapped into a vector following the R conventions).
  #
  #     N <- Number of individuals
  #     t0 <- initial time
  #     t <-  vector of times of snapshot takes (without initial time), assumed to 
  #           be ordered increasingly (t[1] > t0).
  #     d <- spatial dimension
  #     p <- probability of detection in snapshot
  #     
  #     SNAPSHOT REGIONS SPECIFICATION: 
  #     In this public function only an application to same-length snapshot regions is provided.
  #     In addition, the snapshot squares are supposed to be the same for every snapshot time.
  #
  #       Regions.B.L <-  Matrix of dimensions ( n.s , d ), where n.s is the number of snapshot regions.
  #                       Each row indicates the coordinates of the bottom-left of a square snapshot region.      
  #       dx <- Double with the length of the side of the square snapshot regions.
  #
  #     TRAJECTORY SPECIFICATION
  #
  #     Traj <- String indicating the type of trajectory of the individuals."Brownian" and "Brownian+Advection" are the currently available trajectories. 
  #                 If Traj = "Brownian" the individuals move according to a Brownian motion (iid components centered at the origin).
  #             
  #                 If Traj ="Brownian+Advection" the individuals move according to a Brownian motion (iid components centered at the origin) plus a translation according to a velocity vector.
  #
  #     theta.X <-  Vector of parameters for the trajectory model. 
  #                   For "Brownian", a double indicating the standard deviation. 
  #                   For "Brownian+Advection", a vector of size d+1:
  #                           theta.X[1]        <- standard deviation of the Brownian motion 
  #                           theta.x[2:(d+1)]  <- advection velocity vector.
  #
  #Output:
  #
  #   l <- A list of two objects:
  #               l$mean.vec <- the mean vector of a wrapped-into-a-vector random matrix from a Snapshot model.                            
  #               l$cov.mat <- the covariance matrix of a wrapped-into-a-vector random matrix from a Snapshot model.                            
  
  n.t <- length(t)
  n.s <- nrow(Regions.B.L )
  
  P.X.t.site <- matrix( 1 , nrow = n.t , ncol = n.s )
  P.X.ts.sites <- array( 1 , dim = c(n.t,n.t,n.s,n.s) )
  
  if(Traj == "Brownian"){
    Traj <- "Brownian+Advection"
    theta.X <- c( theta.X , rep(0,d) ) #On se casse pas la tête...
  }
  
  if(Traj == "Brownian+Advection"){
    for(  k in 1:n.t ){
      for( l in 1:d ){
        P.X.t.site[k,] <- P.X.t.site[k,]*(
          pnorm( Regions.B.L[,l] + dx ,
                 mean = theta.X[1+l]*(t[k]-t0) , 
                 sd = theta.X[1]*sqrt(t[k]-t0)  ) -
            pnorm(  Regions.B.L[,l],
                    mean = theta.X[1+l]*(t[k]-t0) , 
                    sd = theta.X[1]*sqrt(t[k]-t0)  )
        )
      }
    }
    
    #Covariance matrix of a one-dimensional Brownian motion
    Cov.Brown <- matrix( 0 , n.t , n.t)
    for(k in 1:n.t){
      Cov.Brown[k , k:n.t] <- (theta.X[1]^2)*(t[k] - t0)
      Cov.Brown[k:n.t , k] <- (theta.X[1]^2)*(t[k] - t0)
    }
    
    for(k in 1:n.t){
      if(n.s==1){
        P.X.ts.sites[k,k, , ] <- as.matrix( P.X.t.site[k, ] , nrow = 1 , ncol = 1 )
      }else{
        P.X.ts.sites[k,k, , ] <- diag(  P.X.t.site[k, ] ) 
      }
      if(k > 1){
        for(k2 in 1:(k-1)){
          for(l in 1:n.s){
            for(l2 in 1:n.s){
              P.X.ts.sites[k,k2,l,l2] <- sadmvn(
                lower = c( Regions.B.L[l,] , Regions.B.L[l2,] ) , 
                upper = c( Regions.B.L[l,]+dx , Regions.B.L[l2,]+dx ),
                mean = c( (t[k]-t0)*theta.X[1+(1:d)] , (t[k2]-t0)*theta.X[1+(1:d)]),
                varcov =rbind( cbind( diag(Cov.Brown[k,k],nrow=d) , diag(Cov.Brown[k2,k2],nrow=d) ),
                               cbind( diag(Cov.Brown[k2,k2],nrow=d) , diag(Cov.Brown[k2,k2],nrow=d) )),
                abseps = 1e-8)
            }
          }
        }
      }
    }
  }
  
  C <- matrix( 0 , n.t*n.s , n.t*n.s  )
  for(k in 1:n.t){
    C[n.s*(k-1)+1:n.s , n.s*(k-1)+1:n.s] <- outer(P.X.t.site[k,] , (-N*p^2)*P.X.t.site[k,] )
    if(n.s == 1){
      C[k,k] <- C[k,k] + (p*N)*P.X.t.site[k,]
    }else{
      diag(C[n.s*(k-1)+1:n.s , n.s*(k-1)+1:n.s]) <- diag(C[n.s*(k-1)+1:n.s , n.s*(k-1)+1:n.s]) +
        (p*N)*P.X.t.site[k,]
    }
    if(k > 1){
      for(k1 in 1:(k-1)){
        C[ n.s*(k-1)+1:n.s , n.s*(k1-1)+1:n.s ] <- (N*p^2)*( P.X.ts.sites[k,k1, , ] - 
                                                               outer(P.X.t.site[k,] , P.X.t.site[k1,]  ) )
      }
    }
  }
  C[upper.tri(C)] = t(C)[upper.tri(C)]
  
  l <- list( (N*p)*as.vector(t(P.X.t.site)) , C )
  names(l) <- c("mean.vec","cov.mat")
  return(l)
}


