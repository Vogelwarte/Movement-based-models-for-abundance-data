#Coded by Ricardo Carrizo V
#Version at 22-08-2024

Sim_EcoDiff <- function(N , 
                        t0 , t , d ,
                        p , sigma , v , 
                        Regions.B.L , dx){
  #Function to simulate one instance of an EcoDiff model counting matrix
  #The underlying EDP is supposed to be
  #
  #       du/dt + grad(u)*v - (sigma^2)/2 Lap(u) = 0,    over [t0 , Inf) x R^d,
  #       u( , t0) = delta_0
  #
  #where grad is the spatial gradient, Lap is the spatial Laplacian, sigma > 0 is a diffusivity parameter and v is an advection vector.
  #NOTE: Only the case with independent Poisson marginals is here considered.
  
  
  #Inputs:
  #
  #     N <- Number of individuals
  #     t0 <- initial time
  #     t <-  vector of times of count, assumed to 
  #           be ordered increasingly (t[1] > t0).
  #     d <- spatial dimension
  #     p <- probability of detection
  #     sigma <- diffusivity controlling parameter
  #     v <- advection vector
  #
  #     
  #     COUNT REGIONS SPECIFICATION: 
  #     In this public function only an application to same-length square regions is provided.
  #     In addition, the square regions are supposed to be the same for every time.
  #
  #       Regions.B.L <-  Matrix of dimensions ( n.s , d ), where n.s is the number of count regions.
  #                       Each row indicates the coordinates of the bottom-left of a square region.      
  #       dx <- Double with the length of the side of the square count regions.
  #
  #
  #Output:
  #
  #   Q: A matrix of dimensions (n.s , n.t), with n.t being the number of time steps, containing the counts of detected individuals in each respective regions and times.
  #
  n.t <- length(t)
  n.s <- nrow(Regions.B.L)
  U <- matrix( 1 , nrow = n.s , ncol = n.t)

  for(k in 1:n.t){
    for(j in 1:d){
      U[ ,k] = U[ ,k]*( pnorm(  Regions.B.L[ ,j] + dx , 
                                mean = (t[k]-t0)*v[j] , 
                                sd = sqrt(t[k]-t0)*sigma ) -
                          pnorm(  Regions.B.L[ ,j] , 
                                  mean = (t[k]-t0)*v[j] , 
                                  sd = sqrt(t[k]-t0)*sigma ) )
    }
  }
  
  return( matrix( rpois(n.s*n.t , (N*p)*as.vector(U) ) , 
                  nrow = n.s , 
                  ncol = n.t ) )
}