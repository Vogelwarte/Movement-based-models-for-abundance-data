####---- Script to analyse the results of simulation studies ----
# rm(list = ls())
# gc(full = T)


#Setting the desired model and fitting method to analyse
model <- "Capture"
fit_method <- "MGLE"


#Proposed theoretical value case
theta.theo <- c(5,-1,1,0.1)  
#Proposed quantity of individuals
N <- 1e4


#Loading the simulations
Q.sim <- readRDS(paste(c("Q.sim, " , model , ", theta=(", 
                         theta.theo[1], "," ,
                         theta.theo[2], "," ,
                         theta.theo[3], "," ,
                         theta.theo[4], "), N=", N ,".rds"), collapse = "" ))

#Loading the punctual estimates
est <- readRDS( paste(c( fit_method , ", " , model , ", theta=(", 
                         theta.theo[1], "," ,
                         theta.theo[2], "," ,
                         theta.theo[3], "," ,
                         theta.theo[4], "), N=", N ,".rds"), collapse = "" ))

#Loading the log-likelihood maximal values
LL <- readRDS( paste(c("LL." , fit_method , ", " , model , ", theta=(", 
                       theta.theo[1], "," ,
                       theta.theo[2], "," ,
                       theta.theo[3], "," ,
                       theta.theo[4], "), N=", N ,".rds"), collapse = "" ))

#Loading the Hessians at punctual estimates
Hess <- readRDS( paste(c( fit_method , ".Hess, " , model , ", theta=(", 
                          theta.theo[1], "," ,
                          theta.theo[2], "," ,
                          theta.theo[3], "," ,
                          theta.theo[4], "), N=", N ,".rds"), collapse = "" ))



#Number of simulations in the imported data
N.sim <- nrow(est) #Number of simulations



#We first filter the "erratic" estimations, that is, the cases where the Hessian is not positive definite.
n.erratic <- 0
ind.erratic <- numeric(0)
for(i in 1:N.sim){
  #Important: as positive.definiteness criterion, we require all the eigenvalues to be bigger than 1e-12. 
  if( !all( eigen(Hess[i, , ])$values > 1e-12 ) ){
    print(paste("Warning: not positive variances at simulation " , i))
    print("Eigenvalues of the Hessian:")
    print(eigen(Hess[i, , ])$values)
    n.erratic <- n.erratic + 1
    ind.erratic <- c(ind.erratic , i)
  }
}

#Only non-erratic estimations are considered.
if( n.erratic > 0 ){
  Q.sim <- Q.sim[-ind.erratic , ]
  est <- est[ -ind.erratic , ]
  LL <- LL[ -ind.erratic ]
  Hess <- Hess[ -ind.erratic , , ]
}


N.sim.considered <- 1000 #Number of simulations we want to consider
if(N.sim - n.erratic < N.sim.considered){
  print("WARNING: You are trying to analize more estimations than available and valid.")
  print(paste("Number of desired estimations: " , N.sim.considered))
  print(paste("Number of available estimations: " , N.sim))
  print(paste("Number of erratic estimations: " , n.erratic))
}
#Real final number of analysed simulations
N.sim.considered <- min(N.sim.considered , N.sim - n.erratic) 
print(paste("Final number of effectively analysed simulations: " , N.sim.considered))

#Retaining only the considered ones.
Q.sim <- Q.sim[1:N.sim.considered , ]
est <- est[ 1:N.sim.considered , ]
LL <- LL[ 1:N.sim.considered ]
Hess <- Hess[ 1:N.sim.considered , , ]


##-- ANALYSIS ITSELF --##

#Average of estimations
colMeans(est)

#SD of estimations
apply( est , 2 , sd   )

#covariance matrix of estimations
cov(est)
#correlation matrix of estimations
cor(est)


#Coverage regions in each simulation-estimation

#Initialization
Cover.rect.reg <- array( 0 , dim = c(  N.sim.considered  , 4 , 2  )   ) #For marginal coverage   
Cover.rect.presence <- matrix( F ,   N.sim.considered , 4  ) #For marginal coverage
Cover.ellipse.presence <- rep( F ,  N.sim.considered  ) #For joint ellipsoidal coverage

for(i in 1:N.sim.considered){
  inv.Hess <- solve(Hess[i, , ])
  #Coverage marginal Gauss 95% confidence intervals
  Cover.rect.reg[i ,  , 1] <- est[ i ,  ] - qnorm(0.975)*sqrt( diag(inv.Hess) )
  Cover.rect.reg[i ,  , 2] <- est[ i ,  ] + qnorm(0.975)*sqrt( diag(inv.Hess) )
  #Presence or absence of the theoretical value
  Cover.rect.presence[i , ] <- (theta.theo >= Cover.rect.reg[ i ,  , 1])&(
    theta.theo <= Cover.rect.reg[ i ,  , 2] )
  
  #Presence of absence in the 95% joint Gauss confidence ellipsoid 
  Cover.ellipse.presence[ i ] <- sum((Hess[ i , , ]%*%(est[i, ]-theta.theo))*(
    est[i, ]-theta.theo) )<= qchisq( 0.95 , df = 4 )
  if( !all(  diag(inv.Hess) > 0 ) ){
    print(paste("Warning: not positive variances at simulation " , i))
    Cover.rect.presence[ i , ] <- F
    Cover.ellipse.presence[ i ] <- F
  }
}

#Results:
Cover.rect.presence
Cover.ellipse.presence

#Proportion of presence of the theoretical value at the proposed confidence region.
#(for each parameter separately)
prop.presence.rect <- colMeans( Cover.rect.presence  )
prop.presence.rect


#Proportion of presence of all the theoretical parameters to the proposed confidence rectangle.
overall.prop.presence.rect <- mean( apply(Cover.rect.presence , 1 , all) )
overall.prop.presence.rect


#Proportion of presence of the whole theoretical parameters  in the confidence ellipsoid.
prop.presence.ellipse <- mean(Cover.ellipse.presence)
prop.presence.ellipse

