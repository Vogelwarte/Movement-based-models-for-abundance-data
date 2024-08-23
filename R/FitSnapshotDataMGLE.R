#### ---- MGLE for the Snapshot model over the Rover fly data ---- ####

library(numDeriv)

source( "Gauss_ll.R"  )
source( "Moments_Snapshot.R"   )
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
                 mom.snap <- Moments_Snapshot(N=N, d = d , p=theta[4],t0=t0,t=t,
                                              Regions.B.L = Regions.B.L, dx = dx,
                                              theta.X = theta[1:3])
                 ll <- Gauss_ll( Q.obs ,
                                 m = mom.snap$mean.vec,
                                 C = mom.snap$cov.mat)
                 print(paste( c( "For theta =" , theta , "the approximated Gauss_ll is " , ll) , collapse = " "  ))
                 return(-ll)},
               method = "L-BFGS-B" ,
               lower = inf.lim ,
               upper = sup.lim ,
               control = list( parscale = c(10 , 10 , 10 , 1), #The scale of the parameter p is one order of maginute smaller than the others.
                               ndeps = c(1e-5, 1e-5 , 1e-5 , 1e-5 ) ,
                               trace = 6) ,
               hessian = F  )

#Computing the Hessian at the MGLE
Hess.MGLE <- hessian( function(theta){
  mom.snap <- Moments_Snapshot(N=N,d = d , p=theta[4],t0=t0,t=t,
                               Regions.B.L = Regions.B.L, dx = dx,
                               theta.X = theta[1:3])
  ll <- Gauss_ll( Q.obs,
                  m = mom.snap$mean.vec,
                  C = mom.snap$cov.mat)
  print(paste( c( "For theta =" , theta , "the approximated Gauss_ll is " , ll) , collapse = " "  ))
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

mom.Snapshot.MGLE <- Moments_Snapshot(N=N, d= d , p=MGLE$par[4],t0=t0,t=t,
                                      Regions.B.L = Regions.B.L, dx = dx,
                                      theta.X = MGLE$par[1:3])
start <- Sys.time()
cc.Gauss.ll <- Gauss_ll(Q.obs ,  
                        m = mom.Snapshot.MGLE$mean.vec ,
                        C = mom.Snapshot.MGLE$cov.mat ,
                        Continuity.corr = T)
print(Sys.time()-start)

print("Continuity corrected Gaussian log-likelihood at MGLE:")
print(cc.Gauss.ll[1])
print("Continuity corrected Gaussian likelihood at MGLE:")
print(exp(cc.Gauss.ll) )

