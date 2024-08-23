#### ---- MGLE for the Capture model over the Rover fly data ---- ####

library(numDeriv)

source( "Gauss_ll.R"  )
source( "Moments_Capture.R"   )
source( "RoverFlyData.R") 


#Initial value for the parameters
theta0 <- c(5 , -1 , 1 , 0.1)

#Limits of the parameter domain (chosen through simulations)
inf.lim <- c(0.01 , -10 , -10 , 0.001)
sup.lim <- c(10 , 10 , 10 , 1)

#Transform the data matrix into a vector for using the Gauss_ll function.
Q.obs <- as.vector(Q.obs)


start <- Sys.time()
#Maximum Gaussian Likelihood Estimation using optim
MGLE <- optim( theta0 ,
               function(theta){
                 mom.capt <- Moments_Capture(N = N, d = d , 
                                             theta.X.free = theta[1:3] ,
                                             alpha = theta[4] ,
                                             t0 = t0 , tL = tL , tH = tH , dt = 1/60,
                                             Dc.B.L = Regions.B.L , dx = dx)
                 ll <- Gauss_ll( Q.obs ,
                                 m = mom.capt$mean.vec,
                                 C = mom.capt$cov.mat)
                 #print(paste( c( "For theta =" , theta , "the approximated Gauss_ll is " , ll) , collapse = " "  ))
                 return(-ll)},
               method = "L-BFGS-B" ,
               lower = inf.lim ,
               upper = sup.lim ,
               control = list( parscale = c(10 , 10 , 10 , 1),
                               ndeps = c(1e-5, 1e-5 , 1e-5 , 1e-5 ) ,
                               trace = 6) ,
               hessian = T  )
print(Sys.time()-start)


start <- Sys.time()
Hess.MGLE <- hessian( function(theta){
  mom.capt <- Moments_Capture(N = N, d = d , 
                              theta.X.free = theta[1:3] ,
                              alpha = theta[4] ,
                              t0 = t0 , tL = tL , tH = tH , dt = 1/60,
                              Dc.B.L = Regions.B.L , dx = dx)
  ll <- Gauss_ll( Q.obs ,
                  m = mom.capt$mean.vec,
                  C = mom.capt$cov.mat)
  #print(paste( c( "For theta =" , theta , "the approximated Gauss_ll is " , ll) , collapse = " "  ))
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


inf.coverage <- MGLE$par - qnorm(0.975)*sqrt( diag( inv.Hess.MGLE) ) 
sup.coverage <- MGLE$par + qnorm(0.975)*sqrt( diag( inv.Hess.MGLE) ) 
print("Estimated coverage intervals: ")
print(cbind(inf.coverage , sup.coverage   ))


print("Maximum obtained (continuous) likelihood:")
print(-MGLE$value)

mom.Capture.MGLE <- Moments_Capture(N = N, d = d , 
                                    theta.X.free = MGLE$par[1:3] ,
                                    alpha = MGLE$par[4] ,
                                    t0 = t0 , tL = tL , tH = tH , dt = 1/60,
                                    Dc.B.L = Regions.B.L , dx = dx)
start <- Sys.time()
cc.Gauss.ll <- Gauss_ll(Q.obs ,  
                        m = mom.Capture.MGLE$mean.vec ,
                        C = mom.Capture.MGLE$cov.mat ,
                        Continuity.corr = T)
print(Sys.time()-start)

print("Continuity corrected Gaussian log-likelihood at MGLE:")
print(cc.Gauss.ll[1])
print("Continuity corrected Gaussian likelihood at MGLE:")
print(exp(cc.Gauss.ll) )

