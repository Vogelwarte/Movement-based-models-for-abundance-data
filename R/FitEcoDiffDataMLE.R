#### ---- MLE for the EcoDiff model over the Rover fly data ---- ####

library(numDeriv)

source( "EcoDiff_ll.R"   )
source( "RoverFlyData.R") 


#Initial value for the parameters
theta0 <- c(5 , -1 , 1 , 0.1)

#Limits of the parameter domain (chosen through simulations)
inf.lim <- c(0.01 , -10 , -10 , 0.001)
sup.lim <- c(10 , 10 , 10 , 1)

#Here we do not transform the data matrix into a vector


start <- Sys.time()
#Maximum Gaussian Likelihood Estimation using optim
MGLE <- optim( theta0 ,
               function(theta){
                 ll <- EcoDiff_ll(Q.obs , 
                                  N = N , p = theta[4] , 
                                  sigma = theta[1] , v = theta[2:3] , 
                                  t0 = t0 , t = t , d = d , 
                                  Regions.B.L = Regions.B.L , dx = dx)
                 return(-ll)},
               method = "L-BFGS-B" ,
               lower = inf.lim ,
               upper = sup.lim ,
               control = list( parscale = c(10 , 10 , 10 , 1),
                               ndeps = c(1e-5, 1e-5 , 1e-5 , 1e-5 ) ,
                               trace = 6) ,
               hessian = T  )
Hess.MGLE <- hessian( function(theta){
  ll <- EcoDiff_ll(Q.obs , 
                   N = N , p = theta[4] , 
                   sigma = theta[1] , v = theta[2:3] , 
                   t0 = t0 , t = t , d = d , 
                   Regions.B.L = Regions.B.L , dx = dx)
  return(-ll)} , x = MGLE$par  )
print(Sys.time()-start)




print("Estimated parameters:")
print(MGLE$par)

inv.Hess.MGLE <- solve( Hess.MGLE )
print("Estimated covariance matrix:")
print(inv.Hess.MGLE)

print("Estimated eigenvalues of the covariance matrix: ")
print(eigen(inv.Hess.MGLE)$values )

est.corr.mat <- inv.Hess.MGLE/sqrt( outer(  diag(inv.Hess.MGLE) ,
                                            diag(inv.Hess.MGLE))  )
print("Estimated correlation matrix:")
print( est.corr.mat  )


inf.coverage <- MGLE$par - qnorm(0.975)*sqrt( diag( inv.Hess.MGLE ) ) 
sup.coverage <- MGLE$par + qnorm(0.975)*sqrt( diag( inv.Hess.MGLE ) ) 
print("Estimated coverage intervals: ")
print(cbind(inf.coverage , sup.coverage   ))

print("Maximum obtained likelihood:")
print(-MGLE$value)

