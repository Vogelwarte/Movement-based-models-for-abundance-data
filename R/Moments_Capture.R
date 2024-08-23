#Coded by Ricardo Carrizo V
#Version at 22-08-2024
library(mnormt)

Moments_Capture <- function(N , 
                            t0 , tL , tH , dt,
                            d , 
                            Dc.B.L , dx , 
                            alpha , 
                            free.traj = "Brownian+Advection" , 
                            theta.X.free ){
  
  #Function to compute the mean and the covariance of a Capture model random matrix evaluated at the liberation and horizon time.
  #The mean is a vector and the covariance is a matrix, corresponding to the moments of a Capture model random matrix wrapped into a vector following the R conventions (see Sim_Capture function for the structure of a Capture random matrix).
  
  
  #Inputs:
  #
  #     N <- Number of individuals
  #     t0 <- initial time
  #     tL <- Liberation time (tL > t0)
  #     tH <- Horizon time (tH > tL)
  #     dt <- Discretization time grid (must be such that the quantity of time-steps after and before liberation are integer)
  #     d <- spatial dimension
  #     alpha <- rate of growth of probability of capture
  #     
  #     CAPTURE DOMAIN SPECIFICATION: 
  #     In this public function only an application to capture domains constructed as a union of same-length squares is provided.
  #     IMPORTANT: in this method an approximation method for bi-variate Gaussian probabilities is used, which works better when the capture domain is partitioned in small squares. Therefore, it is suggested to give a fine enough partition, even if the whole capture domain is, for instance, square.
  #
  #       Dc.B.L <- Matrix of dimensions ( n.s , d ), where n.s is the number of squares partitioning the capture domain.
  #                 Each row indicates the coordinates of the bottom-left of a square in the capture domain partition.      
  #       dx <- length of the side of the squares in the capture domain partition.
  #
  #     FREE-TRAJECTORY SPECIFICATION
  #
  #     Traj <- String indicating the type of free-trajectory of the individuals."Brownian" and "Brownian+Advection" are the currently available trajectories. 
  #                 If Traj = "Brownian" the individuals move according to a Brownian motion (iid components centered at the origin).
  #             
  #                 If Traj ="Brownian+Advection" the individuals move according to a Brownian motion (iid components centered at the origin) plus a translation according to a velocity vector.
  #
  #     theta.X.free <-  Vector of parameters for the free-trajectory model. 
  #                       For "Brownian", a double indicating the standard deviation. 
  #                       For "Brownian+Advection", a vector of size d+1:
  #                           theta.X[1]        <- standard deviation of the Brownian motion 
  #                           theta.x[2:(d+1)]  <- advection velocity vector.
  #
  
  
  #Output:
  #
  #   l <- A list of two objects:
  #           l$mean.vec  <- mean vector of the wrapped Capture count random matrix.
  #           l$cov.mat   <- covariance matrix of the wrapped Capture count random matrix.
  #
  #(Here we do not consider the positions of individuals, captured or not, in other times than tL and tH).
  
  
  eps <- 1e-161  #Reference epsilon for small stuff (helps avoiding underflow problems)
  n.s <- nrow(Dc.B.L) #number of squares in the partition of Dc
  
  if(free.traj == "Brownian"){
    free.traj <- "Brownian+Advection"
    theta.X.free <- c( theta.X.free , rep(0 , d)   )  #On se casse pas la tête
  }
  #Naming parameters otherwise for simplification (since only Brownian + advection is essentialy available here)
  sigma <- theta.X.free[1]
  v <- theta.X.free[1+1:d]
  
  #Construction of the time grid
  #Verifying that the time grid contains consistent quantity of time points.
  n.L <- (tL - t0)/dt #Number of time-steps before liberation (tL considered, t0 not)
  if(n.L%%1 != 0 ){
    print("ERROR: in constructed discrete grid, number of times before liberation is not integer.")
    return(NULL)
  }
  n.AL <- (tH-tL)/dt #Number of time-steps after liberation (tH considered, tL not)
  if(n.AL%%1 != 0){
    print("ERROR: in constructed discrete grid, number of times after liberation is not integer.")
    return(NULL)
  }
  t <- seq( from = t0 + dt , to = tH , by = dt )
  n.t <- length(t)
  
  
  #Initialization of the returning object
  l.names <- c("mean.vec","cov.mat")
  l <- list( rep(0 , 2*n.t) , matrix(0 , 2*n.t , 2*n.t)  )
  names(l) <- l.names
  
  
  #Matrix of probabilities of the free-trajectory being in the capture sub-domains
  # at the different time-points.
  # PP.1.S[k,l] = P( X.free( t[k] ) \in Regions[l] )
  # k = 1 , ... , n.t
  # l = 1 , ... , n.s
  PP.1.S <- matrix( 1 , nrow = n.t , ncol = n.s)
  for(j in 1:d){
    PP.1.S <- PP.1.S*(
      matrix( pnorm( rep(Dc.B.L[ , j] + dx, 1 , each = n.t) ,
                     mean = v[j]*(t-t0) ,
                     sd = sigma*sqrt(t-t0) ) , 
              n.t , n.s) -
        matrix( pnorm( rep(Dc.B.L[ , j] , 1 , each = n.t) ,
                       mean = v[j]*(t-t0) ,
                       sd = sigma*sqrt(t-t0) ) , 
                n.t , n.s)
    )
  }
  if(any(PP.1.S < eps)){
    PP.1.S[PP.1.S < eps] <- eps #Positivity correction on PP.1.S
  }
  
  
  #Vector of probabilities of being in the capture domain at the time-points
  #PP.1[k] = P( X.free(t[k]) \in D_c )
  PP.1 <- rowSums( PP.1.S )
  if(any(PP.1 > 1 - eps )){
    #Unitarity correction: 
    #Theoretically (Brownian case with bounded capture domain), all the probabilities 
    #in PP.1 should be strictly smaller than 1.
    #For those which are not, we change the corresponding row of PP.1.S by re-weighting so the sum is 1-eps.
    for(k in which(PP.1 > 1 - eps) ){
      PP.1.S[k, ] <- PP.1.S[k, ]*(  (1-eps)/PP.1[k] )  
      PP.1[k] <- sum(PP.1.S[k,])
    }
  }
  
  #time distances between time-tag points
  t.dists <- dt*1:(n.t-1) 
  
  #Now we compute an auxiliary array A.aux, which shall help us to compute approximated bi-variate Gaussian probabilities without computing explicitly PP.2.S (see pseudo-code).
  #Auxiliary array A.aux, intended to mean:
  # A.aux[k,l2,l1] = P( X.free(t + k*dt) \in A_l1   |   X.free(t) = z_l2   ) 
  # k = 1 , ... , n.t - 1
  # l2 = 1 , ... , n.s
  # l1 = 1 , ... , n.s
  #Here z_l2 means the center of A_l2.
  #
  #In general, this should also depend upon "t", but X.free is supposed
  #Brownian. In such case, Aux[k,l2,l1] does not depend upon t
  #(Markov-homogeneous process)
  A.aux <- array( 1 , dim=c( n.t-1 , n.s , n.s ) )
  z <- Dc.B.L + dx/2
  for( j in 1:d ){
    A.aux <- A.aux*sapply( Dc.B.L[,j] , function(x){
      pnorm( x + dx ,
             mean = t( z[,j] + matrix(  v[j]*t.dists ,
                                        n.s ,
                                        n.t-1 ,
                                        byrow = T)),
             sd = sigma*sqrt(t.dists) ) -
        pnorm( x ,
               mean = t( z[,j] + matrix(  v[j]*t.dists ,
                                          n.s ,
                                          n.t-1 ,
                                          byrow = T)),
               sd = sigma*sqrt(t.dists) )
    } , simplify = "array"  )
  }
  
  
  #Auxiliary array PP.cond.P.1
  #PP.cond.P.1 is intended to mean:
  # PP.cond.P.1[k,l] = P( X.free(t+k*dt) \in D_c | X.free(t) = z_l  ),
  # k = 1,...,n.t-1,
  # l = 1,...,n.s
  #It does not depend upon t in our Markov-homogeneous case
  PP.cond.P.1 <- rowSums( A.aux , dims = 2  )
  PP.cond.P.1[eps >1-PP.cond.P.1] <- 1-eps #Unitarity correction
  
  
  #PP.2 matrix. It is lower-triangular.
  #PP.2[k1,k2] = P( X.free(t[k1]) in Dc & X.free(t[k2]) in Dc )
  #k1 = 1 ... n.t
  #k2 = 1 ... k1
  PP.2 <- matrix( 0 , n.t , n.t)
  diag(PP.2) <- PP.1
  for(k2 in 1:(n.t-1) ){
    PP.2[ (k2 + 1):n.t , k2 ] <-  PP.cond.P.1[ 1:(n.t-k2)  , ]%*%PP.1.S[k2 , ]
  }
  #Comparison with PP.1 correction.
  #It should hold PP.2[k1,k2] <= min( PP.1[k1] , PP.1[k2] )
  out.min <- outer( PP.1 , PP.1 , pmin)
  bool.not.valid <- PP.2 > out.min
  if( any(bool.not.valid)){
    PP.2[bool.not.valid] <- out.min[bool.not.valid]
  }
  
  
  #Conditional-to-past probabilities of presence lower-triangular matrix.
  #CPP[k1,k2] = P( X.free(t[k1]) in Dc |  X.free(t[k2]) in Dc  )
  #k1 = 1 ... n.t
  #k2 = 1 ... k1
  CPP <- t(t(PP.2)/PP.1)
  
  
  #The density f.Tc
  #Solving the Volterra equation for f.Tc
  M <- (alpha*dt)*CPP 
  diag(M) <- diag(M) + rep(1,n.t)
  f.Tc <- forwardsolve( M , alpha*PP.1 )
  
  
  #phi.TC function 
  phi.Tc <- f.Tc/PP.1
  
  
  #COMPUTATION OF THE MOMENTS
  
  
  #Initialization
  cov.mat <- matrix( 0 , 2*n.s , 2*n.s )
  mean.vec <- rep( 0 , 2*n.s)
  
  #Non-diagonal blocks of the covariance matrix
  phi.coef <- ( N*alpha*(dt^2) )*phi.Tc #This may make things faster
  for(k2 in 1:n.L){
    cov.mat[1:n.s , n.s + 1:n.s] <- cov.mat[1:n.s , n.s + 1:n.s] +
      colSums(A.aux[1:n.AL ,  , ]*phi.coef[k2 + 1:n.AL])*PP.1.S[k2, ]
  }
  cov.mat[n.s+1:n.s , 1:n.s] <- t(cov.mat[1:n.s , n.s+1:n.s])
  
  
  #Mean vector computation
  mean.vec[1:n.s] <- (N*dt)*colSums(PP.1.S[1:n.L, ]*phi.Tc[1:n.L])
  mean.vec[n.s + 1:n.s] <- (N*dt)*colSums( 
    PP.1.S[(n.L+1):n.t , ]*phi.Tc[(n.L+1):n.t] ) + 
    colSums(cov.mat[  1:n.s , n.s + 1:n.s ])
  
  
  #The diagonal blocks of the (non-centered) second order matrix are diagonal.
  diag(cov.mat) <- mean.vec
  
  
  #Until now, cov.mat is the non-centered second order matrix.
  #The following transforms it into the properly speaking covariance matrix.
  extra <- mean.vec/sqrt(N)
  cov.mat <- cov.mat - outer( extra , extra  )
  
  
  l$mean.vec <- mean.vec
  l$cov.mat <- cov.mat 
  return(l)
}



# #Example:
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
# 
# mom <- Moments_Capture( N = N ,
#                    t0 = t0 , tL = tL , tH = tH , dt = dt ,
#                    d = d ,
#                    Dc.B.L = Dc.B.L , dx = dx ,
#                    alpha = alpha  ,
#                    theta.X.free = theta.X.free  )
# 
# mom$mean.vec
# mom$cov.mat


