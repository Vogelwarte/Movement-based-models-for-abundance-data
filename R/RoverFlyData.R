###Rover fly data

#The data was obtained from the collaborators of Edelsparre et al. (2021).

#For sake of data protection, here we provide only the data concerning the Rover strain flies at times t = 0.5 and t = 1.5.
#This is all that is needed to reproduce the results in Carrizo Vergara et al. (2024).

#SPATIAL SETTING
d <- 2
n.s <- 227
Regions.B.L <- matrix( 0 , n.s , d )
Regions.B.L[ , 1] <- c(-2.5,2.5,2.5,-2.5,-7.5,-7.5,-7.5,-2.5,2.5,7.5,-2.5,-12.5,-2.5,12.5,12.5,2.5,-2.5
                       ,-7.5,-17.5,-17.5,-17.5,-7.5,-2.5,2.5,12.5,17.5,17.5,17.5,17.5,12.5,2.5,-2.5,-7.5,-17.5
                       ,-22.5,-22.5,-22.5,-22.5,-22.5,-22.5,-22.5,-17.5,-7.5,-2.5,2.5,12.5,17.5,17.5,17.5,22.5,22.5
                       ,22.5,22.5,22.5,17.5,12.5,2.5,-2.5,-7.5,-17.5,-22.5,-27.5,-27.5,-27.5,-27.5,-27.5,-27.5,-27.5
                       ,-27.5,-27.5,-22.5,-17.5,-7.5,-2.5,2.5,12.5,17.5,22.5,22.5,22.5,22.5,27.5,27.5,27.5,27.5
                       ,17.5,12.5,2.5,-2.5,-7.5,-17.5,-22.5,-32.5,-32.5,-32.5,-32.5,-32.5,-32.5,-32.5,-22.5,-17.5,-7.5
                       ,-2.5,2.5,12.5,17.5,27.5,27.5,27.5,32.5,32.5,32.5,32.5,17.5,12.5,2.5,-2.5,-7.5,-17.5
                       ,-22.5,-37.5,-37.5,-37.5,-37.5,-37.5,-37.5,-37.5,-22.5,-17.5,-7.5,-2.5,2.5,12.5,17.5,32.5,32.5
                       ,32.5,37.5,37.5,37.5,37.5,37.5,32.5,17.5,12.5,-2.5,-17.5,-22.5,-37.5,-42.5,-42.5,-42.5,-42.5
                       ,-42.5,-42.5,-42.5,-42.5,-42.5,-37.5,-22.5,-17.5,-2.5,12.5,17.5,32.5,37.5,37.5,37.5,37.5,37.5
                       ,32.5,17.5,-2.5,-22.5,-37.5,-42.5,-47.5,-47.5,-47.5,-47.5,-47.5,-47.5,-47.5,-47.5,37.5,32.5,17.5
                       ,-2.5,-22.5,-37.5,-42.5,-52.5,-52.5,-52.5,-52.5,-52.5,-52.5,-52.5,37.5,32.5,17.5,-2.5,-22.5,-37.5
                       ,-42.5,-57.5,-57.5,-57.5,-57.5,-57.5,-57.5,-57.5,37.5,32.5,17.5,-22.5,-37.5,-42.5,-57.5,-62.5,-62.5
                       ,-62.5,-62.5,-62.5,-62.5,-62.5,-62.5)
Regions.B.L[ , 2] <- c(-2.5,-2.5,2.5,2.5,2.5,-2.5,-7.5,-7.5,-7.5,-2.5,7.5,-2.5,-12.5,-2.5,2.5,12.5,12.5
                       ,12.5,2.5,-2.5,-7.5,-17.5,-17.5,-17.5,-7.5,-2.5,2.5,12.5,17.5,17.5,17.5,17.5,17.5,17.5
                       ,17.5,12.5,2.5,-2.5,-7.5,-17.5,-22.5,-22.5,-22.5,-22.5,-22.5,-22.5,-22.5,-17.5,-7.5,-2.5,2.5
                       ,12.5,17.5,22.5,22.5,22.5,22.5,22.5,22.5,22.5,22.5,22.5,17.5,12.5,2.5,-2.5,-7.5,-17.5
                       ,-22.5,-27.5,-27.5,-27.5,-27.5,-27.5,-27.5,-27.5,-27.5,-27.5,-22.5,-17.5,-7.5,-2.5,2.5,12.5,17.5
                       ,27.5,27.5,27.5,27.5,27.5,27.5,27.5,17.5,12.5,2.5,-2.5,-7.5,-17.5,-22.5,-32.5,-32.5,-32.5
                       ,-32.5,-32.5,-32.5,-32.5,-22.5,-17.5,-7.5,-2.5,2.5,12.5,17.5,32.5,32.5,32.5,32.5,32.5,32.5
                       ,32.5,17.5,12.5,2.5,-2.5,-7.5,-17.5,-22.5,-37.5,-37.5,-37.5,-37.5,-37.5,-37.5,-37.5,-22.5,-17.5
                       ,-7.5,-2.5,12.5,17.5,32.5,37.5,37.5,37.5,37.5,37.5,37.5,37.5,37.5,37.5,32.5,17.5,12.5
                       ,-2.5,-17.5,-22.5,-37.5,-42.5,-42.5,-42.5,-42.5,-42.5,-42.5,-42.5,-42.5,-42.5,-37.5,-22.5,-17.5,42.5
                       ,42.5,42.5,42.5,42.5,42.5,42.5,42.5,37.5,32.5,17.5,-2.5,-22.5,-37.5,-42.5,47.5,47.5,47.5
                       ,47.5,47.5,47.5,47.5,37.5,32.5,17.5,-2.5,-22.5,-37.5,-42.5,52.5,52.5,52.5,52.5,52.5,52.5
                       ,52.5,37.5,32.5,17.5,-2.5,-22.5,-37.5,-42.5,57.5,57.5,57.5,57.5,57.5,57.5,57.5,57.5,52.5
                       ,37.5,32.5,17.5,-22.5,-37.5,-42.5)
dx <- 5

#TIME SETTING
t0 <- 0
tL <- 0.5
tH <- 1.5
t <- c(tL , tH)
n.t <- length(t)


#Count Data
N <- 5644 #Number of individuals

Q.obs <- matrix( 0 , n.s , n.t )

Q.obs[1,1] <- 55
Q.obs[4,1] <- 4
Q.obs[6,1] <- 1

Q.obs[1,2] <- 87
Q.obs[2,2] <- 2
Q.obs[6,2] <- 3
Q.obs[8,2] <- 1
Q.obs[10,2] <- 1
Q.obs[13,2] <- 1
Q.obs[18,2] <- 3
Q.obs[19,2] <- 1
Q.obs[20,2] <- 1
