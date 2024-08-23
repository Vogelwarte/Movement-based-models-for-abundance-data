#Coded by Ricardo Carrizo V
#Version at 22-08-2024

library(mnormt)

Sim_Snapshot <- function( N ,
                          t0 , t  , 
                          d, 
                          Regions.B.L, dx, 
                          p, 
                          Traj = "Brownian+Advection", theta.X){
  #Function to simulate one instance of a Snapshot model counting matrix
  
  #Inputs:
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
  #   Q: A matrix of dimensions (n.s , n.t), with n.t being the number of time steps, containing the counts of detected individuals in each respective regions and times.
  #
  
  n.t <- length(t)
  n.s <- nrow(Regions.B.L )
  
  if(Traj == "Brownian"){
    Traj <- "Brownian+Advection"
    theta.X <- c( theta.X , rep(0,d) ) #on se casse pas la tête...
  }
  
  
  if(Traj == "Brownian+Advection"){
    
    #Covariance matrix of a one-dimensional Brownian motion
    Cov.Brown <- matrix( 0 , n.t , n.t)
    for(k in 1:n.t){
      Cov.Brown[k , k:n.t] <- (theta.X[1]^2)*(t[k] - t0)
      Cov.Brown[k:n.t , k] <- (theta.X[1]^2)*(t[k] - t0)
    }
    Chol.Brown <- chol(Cov.Brown)
    
    #Trajectories simulation
    X <- array( 0 , dim = c( N , d ,  n.t  ) )
    for(l in 1:d){
      X[ , l , ] <- rmnorm(N , 
                           mean = theta.X[1+l]*(t - t0),
                           sqrt = Chol.Brown )
    }
    
    #Counting in each space-time square
    Q <- matrix(0 , nrow = n.s , ncol = n.t )
    for(k in 1:n.t){
      Belonging.mat <- matrix(T , nrow = N , ncol = n.s)
      for(l in 1:d){
        Belonging.mat <- Belonging.mat&outer(X[ ,l,k],
                                             Regions.B.L[,l],
                                             function(x,y){x > y})
        Belonging.mat <- Belonging.mat&outer(X[,l,k],
                                             Regions.B.L[,l]+dx,
                                             function(x,y){x <= y})
      }
      Q[,k] <- colSums(Belonging.mat)
    }
    
    #Snapshot detection error
    Q <- rbinom( n.s*n.t , as.vector(Q) , p  )
    
    return( matrix(Q , nrow = n.s , ncol = n.t  ) )
  }
}


# #Example:
# N <- 1000
# t0 <- 0
# t <- c(0.5 , 2.5 , 4 , 10)
# d <- 2
# Regions.B.L = rbind(  c(0,0) , c(1,0) , c(-1,1) , c(-1 , -1) , c(0,1))
# dx <- 1
# theta.X <- c(2,-1 , 1)
# p <- 0.8
# 
#
# Q <- Sim_Snapshot( N = N , 
#                    t0 = t0 , t = t , 
#                    d = d , 
#                    Regions.B.L = Regions.B.L , 
#                    dx = dx , 
#                    p = p  , 
#                    theta.X = theta.X  )








