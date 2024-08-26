#### ---- Simulation studies for MLE on the EcoDiff model ---- ####

#This code uses parallel computing.
library(parallel)

#Select here the number of cores to use 
detectCores()
n.cores <- 3
cl <- makeCluster(n.cores)

clusterEvalQ(cl , {
  library(numDeriv)
  
  source( "Sim_EcoDiff.R" )
  source( "EcoDiff_ll.R" )
  source( "RoverFlyData.R") #The data is imported to obtain the space-time trap configuration.
  
  #Set here the number of simulations per core
  N.sim <- 350
  
  #Change here to propose the quantity of individuals
  N <- 1e5
  
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
seeds <- 13*( 2*(1:n.cores) - 1 ) - 2*(0:(n.cores-1)) 


start <- Sys.time()
L <- parLapply(cl , seeds , 
               function(s){
                 #This is the function to be executed en each core
                 
                 
                 #Initialization of the matrix containing the counting vectors for each simulation.
                 Q.sim <- matrix( 0 , N.sim , 2*n.s )
                 
                 #Initialization of the matrix containing the MLEs
                 MLE <- matrix( 0 , N.sim , 4)
                 
                 #Initialization of the vector containing the obtained maximal log-likelihoods
                 LL.MLE <- rep( 0 , N.sim )
                 
                 #Initialization of the array containing the Hessians at the MLE
                 MLE.Hess <- array( 0 , c(N.sim , 4 , 4) )
                 
                 #Setting the seed
                 set.seed(s)
                 for(i in 1:N.sim){
                   q <- Sim_EcoDiff(N = N ,d = d , 
                                    sigma = theta.theo[1] , 
                                    v = theta.theo[2:3],
                                    p = theta.theo[4] ,
                                    t0 = t0 , t = t,
                                    Regions.B.L = Regions.B.L  , dx = dx ) 
                   Q.sim[i , ] <- as.vector(q)
                   
                   mle <- optim( theta.0 ,
                                 function(theta){
                                   ll <- EcoDiff_ll(q ,
                                                    N = N , p = theta[4] , 
                                                    sigma = theta[1] , v = theta[2:3] , 
                                                    t0 = t0 , t = t , d = d , 
                                                    Regions.B.L = Regions.B.L , dx = dx)
                                   return(-ll)},
                                 method = "L-BFGS-B" ,
                                 lower = inf.lim ,
                                 upper = sup.lim ,
                                 control = list( parscale = c(10 , 10 , 10 , 1),
                                                 ndeps = c(1e-5, 1e-5 , 1e-5 , 1e-5 )) )
                   MLE[i , ] <- mle$par
                   LL.MLE[i] <- -mle$value
                   MLE.Hess[i , , ] <- hessian( function(theta){
                     ll <- EcoDiff_ll(q ,
                                      N = N , p = theta[4] , 
                                      sigma = theta[1] , v = theta[2:3] , 
                                      t0 = t0 , t = t , d = d , 
                                      Regions.B.L = Regions.B.L , dx = dx)
                     return(-ll)} , 
                     x = mle$par   )
                 }
                 l <- list(Q.sim, MLE , LL.MLE , MLE.Hess , N , theta.theo , N.sim)
                 names(l) <- c("Q.sim" , "MLE" , "LL.MLE" , "MLE.Hess" , "N" , "theta.theo" , "N.sim.per.core")
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
MLE <- matrix( 0 , N.sim.per.core*n.cores , 4 )
LL.MLE <- rep( 0 , N.sim.per.core*n.cores)
MLE.Hess <- array( 0 , c(N.sim.per.core*n.cores , 4 , 4 ) )
for(i in 1:n.cores){
  Q.sim[ (i-1)*N.sim.per.core + 1:N.sim.per.core , ] <- L[[i]]$Q.sim
  MLE[ (i-1)*N.sim.per.core + 1:N.sim.per.core , ] <- L[[i]]$MLE
  LL.MLE[(i-1)*N.sim.per.core + 1:N.sim.per.core] <- L[[i]]$LL.MLE
  MLE.Hess[ (i-1)*N.sim.per.core + 1:N.sim.per.core , , ] <- L[[i]]$MLE.Hess 
}


#Writing down the results

#Saving the simulations (WATCH OUT, THIS WILL OVERWRITTE THE SIMULATION RESULTS WHEN USING OTHER FITTING METHODS)
#If the seeds are the same, there should be no difference, since the same Sim_EcoDiff function is used to simulate.
saveRDS( Q.sim , paste(c("Q.sim, EcoDiff, theta=(", 
                         theta.theo[1], "," ,
                         theta.theo[2], "," ,
                         theta.theo[3], "," ,
                         theta.theo[4], "), N=", N ,".rds"), collapse = "" ))
#Saving the MLE punctual estimates
saveRDS( MLE , paste(c("MLE, EcoDiff, theta=(", 
                       theta.theo[1], "," ,
                       theta.theo[2], "," ,
                       theta.theo[3], "," ,
                       theta.theo[4], "), N=", N ,".rds"), collapse = "" ))
#Saving the log-likelihood maximal values
saveRDS( LL.MLE , paste(c("LL.MLE, EcoDiff, theta=(", 
                          theta.theo[1], "," ,
                          theta.theo[2], "," ,
                          theta.theo[3], "," ,
                          theta.theo[4], "), N=", N ,".rds"), collapse = "" ))
#Saving the Hessians at MLE
saveRDS( MLE.Hess , paste(c("MLE.Hess, EcoDiff, theta=(", 
                            theta.theo[1], "," ,
                            theta.theo[2], "," ,
                            theta.theo[3], "," ,
                            theta.theo[4], "), N=", N ,".rds"), collapse = "" ))

