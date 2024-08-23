#### ---- Simulation studies for MGLE on the Capture model ---- ####

#This code uses parallel computing.
library(parallel)

#Select here the number of cores to use 
detectCores()
n.cores <- 3
cl <- makeCluster(n.cores)

clusterEvalQ(cl , {
  library(numDeriv)
  
  source( "Gauss_ll.R" )
  source( "Sim_Capture.R" )
  source( "Moments_Capture.R" )
  source( "RoverFlyData.R") #The data is imported to obtain the space-time trap configuration.
  
  #Set here the number of simulations per core
  N.sim <- 2
  
  #Change here to propose the quantity of individuals
  N <- 1e4
  
  #Change here to propose the theoretical setting for the parameters
  theta.theo <- c(5,-1, 1, 0.1) 
  
  #Some not-so-different-from-the-theoretical initial value
  theta.0 <- c( theta.theo[1]*1.2 , 
                0 , 
                0 , 
                min( 1 , theta.theo[4]*1.5) )
  
  #Optimization limits for the parameters, chosen artisanally
  inf.lim <- c(0.01 , -10 , -10 , 0.001)
  sup.lim <- c(2*theta.theo[1] , 10 , 10 , 1)
}
)

#Set seeds (different for each core)
seeds <- 13^(1:n.cores) - 0:(n.cores-1) 


start <- Sys.time()
L <- parLapply(cl , seeds , 
               function(s){
                 #This is the function to be executed en each core
                 
                 
                 #Initialization of the matrix containing the counting vectors for each simulation.
                 Q.sim <- matrix( 0 , N.sim , 2*n.s )
                 
                 #Initialization of the matrix containing the MGLEs
                 MGLE <- matrix( 0 , N.sim , 4)
                 
                 #Initialization of the vector containing the obtained maximal Gauss log-likelihoods
                 LL.MGLE <- rep( 0 , N.sim )
                 
                 #Initialization of the array containing the Hessians at the MGLE
                 MGLE.Hess <- array( 0 , c(N.sim , 4 , 4) )
                 
                 #Setting the seed
                 set.seed(s)
                 for(i in 1:N.sim){
                   Q.sim[ i , ] <- as.vector(
                     Sim_Capture(N = N ,d = d , 
                                 theta.X.free = theta.theo[1:3] ,
                                 alpha = theta.theo[4] ,
                                 t0 = t0 , tL = tL , tH = tH , dt = 1/60,
                                 Dc.B.L = Regions.B.L  , dx = dx ) 
                   )
                   
                   mgle <- optim( theta.0 ,
                                  function(theta){
                                    mom.capt <- Moments_Capture(N = N, d = d , 
                                                                theta.X.free = theta[1:3] ,
                                                                alpha = theta[4] ,
                                                                t0 = t0 , tL = tL , tH = tH , dt = 1/60,
                                                                Dc.B.L = Regions.B.L , dx = dx)
                                    ll <- Gauss_ll( Q.sim[i , ] ,
                                                    m = mom.capt$mean.vec,
                                                    C = mom.capt$cov.mat)
                                    return(-ll)},
                                  method = "L-BFGS-B" ,
                                  lower = inf.lim ,
                                  upper = sup.lim ,
                                  control = list( parscale = c(10 , 10 , 10 , 1),
                                                  ndeps = c(1e-5, 1e-5 , 1e-5 , 1e-5 )) )
                   MGLE[i , ] <- mgle$par
                   LL.MGLE[i] <- -mgle$value
                   MGLE.Hess[i , , ] <- hessian( function(theta){
                     mom.capt <- Moments_Capture(N = N, d = d , 
                                                 theta.X.free = theta[1:3] ,
                                                 alpha = theta[4] ,
                                                 t0 = t0 , tL = tL , tH = tH , dt = 1/60,
                                                 Dc.B.L = Regions.B.L , dx = dx)
                     ll <- Gauss_ll( Q.sim[i , ] ,
                                     m = mom.capt$mean.vec,
                                     C = mom.capt$cov.mat)
                     return(-ll)} , 
                     x = mgle$par   )
                 }
                 l <- list(Q.sim, MGLE , LL.MGLE , MGLE.Hess , N , theta.theo , N.sim)
                 names(l) <- c("Q.sim" , "MGLE" , "LL.MGLE" , "MGLE.Hess" , "N" , "theta.theo" , "N.sim.per.core")
                 gc()
                 return(l)
               } )
stopCluster(cl)
gc()
print(Sys.time()-start)



##Merging the results.

N.sim.per.core <- L[[1]]$N.sim.per.core
n.s <- 227
N <- L[[1]]$N
theta.theo <- L[[1]]$theta.theo


Q.sim <- matrix(0 , N.sim.per.core*n.cores , 2*n.s  )
MGLE <- matrix( 0 , N.sim.per.core*n.cores , 4 )
LL.MGLE <- rep( 0 , N.sim.per.core*n.cores)
MGLE.Hess <- array( 0 , c(N.sim.per.core*n.cores , 4 , 4 ) )
for(i in 1:n.cores){
  Q.sim[ (i-1)*N.sim.per.core + 1:N.sim.per.core , ] <- L[[i]]$Q.sim
  MGLE[ (i-1)*N.sim.per.core + 1:N.sim.per.core , ] <- L[[i]]$MGLE
  LL.MGLE[(i-1)*N.sim.per.core + 1:N.sim.per.core] <- L[[i]]$LL.MGLE
  MGLE.Hess[ (i-1)*N.sim.per.core + 1:N.sim.per.core , , ] <- L[[i]]$MGLE.Hess 
}


#Writing down the results

#Saving the simulations
saveRDS( Q.sim , paste(c("Q.sim, Capture, theta=(", 
                         theta.theo[1], "," ,
                         theta.theo[2], "," ,
                         theta.theo[3], "," ,
                         theta.theo[4], "), N=", N ,".rds"), collapse = "" ))
#Saving the MGL punctual estimates
saveRDS( MGLE , paste(c("MGLE, Capture, theta=(", 
                        theta.theo[1], "," ,
                        theta.theo[2], "," ,
                        theta.theo[3], "," ,
                        theta.theo[4], "), N=", N ,".rds"), collapse = "" ))
#Saving the Gauss log-likelihood maximal values
saveRDS( LL.MGLE , paste(c("LL.MGLE, Capture, theta=(", 
                           theta.theo[1], "," ,
                           theta.theo[2], "," ,
                           theta.theo[3], "," ,
                           theta.theo[4], "), N=", N ,".rds"), collapse = "" ))
#Saving the Hessians at MGLE
saveRDS( MGLE.Hess , paste(c("MGLE.Hess, Capture, theta=(", 
                             theta.theo[1], "," ,
                             theta.theo[2], "," ,
                             theta.theo[3], "," ,
                             theta.theo[4], "), N=", N ,".rds"), collapse = "" ))


