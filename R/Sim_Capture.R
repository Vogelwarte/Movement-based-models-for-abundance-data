#Coded by Ricardo Carrizo V
#Version at 22-08-2024

Sim_Capture <- function( N , 
                         t0 , tL , tH , dt,
                         d = 2 , 
                         Dc.B.L = NULL , dx = NULL , 
                         alpha , 
                         free.traj = "Brownian+Advection" , 
                         theta.X.free){
  #Function to simulate one instance of a Capture model counting matrix over a liberation time and an horizon time.
  
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
  #   Q <- matrix with the counts of captured individuals in each square of the partition of the capture domain at the liberation time and at the horizon time.
  #
  #(Here we do not simulate the positions of individuals, captured or not, in other times than tL and tH).
  
  
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
  
  #Initialization of the returning matrix
  Q <- matrix(0 , n.s , 2 )
  
  
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
  
  
  ##SIMULATION OF THE CAPTURE TIMES
  
  #Probabilities associated to the capture time.
  P.Tc.1 <- c( dt*f.Tc , 1 - sum(dt*f.Tc) ) #NOTE: no unitarity verification so far, watch out for complicated values of the inputs
  #Simulation of the indexes of the capture times. 
  Tc.1.ind <- sample.int( n.t + 1 , size = N , replace = T ,
                          prob = P.Tc.1 ) #NOTE: if the index is n.t+1, it means it is not captured within the horizon time.
  
  
  #Indices of individuals captured at least once
  ind.captured <- which(Tc.1.ind <= n.t)
  
  
  #If no-one was captured, there is nothing else to compute.
  if(length(ind.captured) == 0){
    return(Q)
  }
  
  
  #indices of individuals which were captured before tL
  ind.capt.bef.tL <- which(Tc.1.ind <= n.L)
  #indices of individuals which were captured for the first time after tL
  ind.capt.1.aft.tL <- setdiff(ind.captured , ind.capt.bef.tL)
  
  
  #Initialization of vectors of indexes of position of capture
  #capt.pos.1[i] = 0, means the individual i was not captured
  #capt.pos.1[i] in {1 , ... , n.s} AND capt.pos.2[i] = 0 means i was captured once but not twice
  capt.pos.1 <- rep( 0 , N )
  capt.pos.2 <- rep( 0 , N )
  
  
  #Capture Position for those who where captured for the first time after tL
  if(length(ind.capt.1.aft.tL) > 0){
    #For these individuals there is no need to simulate the second capture time.
    capt.pos.1[ind.capt.1.aft.tL] <- sapply( ind.capt.1.aft.tL ,
                                             function(i){
                                               sample.int(n.s , 
                                                          1 , 
                                                          replace = T , 
                                                          prob = PP.1.S[
                                                            Tc.1.ind[i] , ]
                                               )
                                             },
                                             simplify = "vector")
    
    #Addition to the Q matrix at the second time
    aux <- tabulate( capt.pos.1[ind.capt.1.aft.tL] )
    Q[ 1:length(aux) , 2] <- Q[ 1:length(aux) , 2] + aux
  }
  
  
  #If no one was capture before tL
  if(length(ind.capt.bef.tL)==0){
    #Everybody was captured after tL, and their positions have already been computed and aggregated to the Q matrix. There is nothing else to do.
    return(Q)
  }
  
  
  #Individuals who where captured before tL so they could get captured again.
  
  #Auxiliary function phi.Tc
  phi.Tc <- f.Tc/PP.1
  
  
  #Matrix of conditional probabilities for lag to second capture time
  #(we do not compute the density f.Tc.2, but F.Tc.2 = dt*f.Tc.2 directly)
  F.Tc.2 <- matrix( 0 , nrow = n.AL, ncol = n.L )
  for(k2 in 1:n.L){
    F.Tc.2[ , k2] <- ( (dt*alpha)/phi.Tc[k2])*( 
      CPP[ k2 + 1:n.AL , k2]*phi.Tc[k2 + 1:n.AL] )
  }
  F.Tc.2 <- rbind( F.Tc.2 , 1 - colSums(F.Tc.2) )
  
  
  #Simulation of the second capture times
  Tc.2.ind <- rep(Inf , N) #Initialization
  #probabilities for the lag time for each individual
  P.Tc.2 <- as.matrix( F.Tc.2[ ,Tc.1.ind[ind.capt.bef.tL]] , 
                       nrow = n.AL+1 , ncol = length(ind.capt.bef.tL) ) 
  #Simulation
  Tc.2.ind[ind.capt.bef.tL] <- n.L + sapply(  1:length(ind.capt.bef.tL) , 
                                              function(i){ sample.int(n.AL+1 , 1 , replace = T , 
                                                                      prob = P.Tc.2[,i]) } ,  
                                              simplify = "vector" )
  
  
  #indices of individuals captured twice
  ind.cap.twice <- which(Tc.2.ind <= n.t)
  
  
  #Capture positions for individuals captured twice:
  if(length(ind.cap.twice) > 0 ){
    #index for the two capture positions within a matrix (index from 1 to n.s^2)
    mat.ind.pos <- sapply(  ind.cap.twice , 
                            function(i){
                              sample.int( n.s^2 , 
                                          1 ,
                                          replace = T ,
                                          prob = as.vector(
                                            t(A.aux[ Tc.2.ind[i] - n.L ,  ,  ]*PP.1.S[
                                              Tc.1.ind[i] ,  ])
                                          )  )
                            } , 
                            simplify = "vector")
    
    #Capture position at the first time of capture
    capt.pos.1[ind.cap.twice] <- (mat.ind.pos-1)%/%n.s + 1
    #Capture position at the second time of capture
    capt.pos.2[ind.cap.twice] <- (mat.ind.pos-1)%%n.s + 1
    
    #Aggregation to the matrix Q
    aux <- tabulate( capt.pos.1[ind.cap.twice] )
    Q[ 1:length(aux) , 1] <- Q[ 1:length(aux) , 1] + aux
    aux <- tabulate( capt.pos.2[ind.cap.twice] )
    Q[ 1:length(aux) , 2] <- Q[ 1:length(aux) , 2] + aux
  }
  
  
  #Capture position for individuals captured only once, before tL
  
  #Indexes of individuals captured only once, before tL
  ind.capt.1.bef.tL <- setdiff( ind.capt.bef.tL , ind.cap.twice )
  
  if(length(ind.capt.1.bef.tL)>0){
    #Auxiliary array G.aux
    #G.aux[k1,k2,l] = P(  X.free(t[k2]) in A_l  &  X.free(t[k1]) in Dc  )*phi.Tc[k1]
    #k1 = 1 ... n.t
    #k2 = 1 ... min(k1,n.L)    (only k2 <= n.L is necessary, if k2 > k1 the value is not used)
    #l = 1 ... n.s
    G.aux <- array( 0 , dim = c( n.t , n.L , n.s )  )
    for(k2 in 1:min(n.L,n.t-1) ){
      G.aux[k2,k2,  ] <- phi.Tc[k2]*PP.1.S[k2, ]
      G.aux[(k2+1):n.t,k2, ] <- t( t(PP.cond.P.1[1:(n.t-k2), ])*PP.1.S[k2, ] )*phi.Tc[(k2+1):n.t]
    }
    
    #Simulation of the positions of capture
    capt.pos.1[ind.capt.1.bef.tL] <- sapply( ind.capt.1.bef.tL , 
                                             function(i){
                                               sample.int( n.s , 
                                                           size = 1 , 
                                                           replace = T , 
                                                           prob = PP.1.S[Tc.1.ind[i], ]-(alpha*dt/phi.Tc[Tc.1.ind[i]])*
                                                             colSums(G.aux[Tc.1.ind[i] + 1:n.AL , Tc.1.ind[i], ] )   
                                               )
                                             },
                                             simplify = "vector")
    #Aggregation to the matrix Q
    aux <- tabulate( capt.pos.1[ind.capt.1.bef.tL] )
    Q[1:length(aux) , 1] <- Q[1:length(aux) , 1] + aux
  }
  
  return(Q)
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
# Q <- Sim_Capture( N = N ,
#                    t0 = t0 , tL = tL , tH = tH , dt = dt ,
#                    d = d ,
#                    Dc.B.L = Dc.B.L , dx = dx ,
#                    alpha = alpha  ,
#                    theta.X.free = theta.X.free  )




