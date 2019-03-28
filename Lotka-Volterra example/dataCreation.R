
# First we create some data #
#############################

# Equilibrium points are alpha/beta for predators and gamma/delta for preys
populationTop<-120
alpha <- 0.1 # birth rate for prey
beta <- 2/populationTop*alpha # 0.005 # predator eating prey
gamma <- 0.1 # death of predator
delta <- 2/populationTop*gamma # 0.005 # birth rate predator
seasonality <- function(t){
  #1.5 + cos(2*pi/100*t)/2
  1
}

# Dominating rate
domRate <- 1 * (populationTop * (alpha+gamma) + (beta + delta) * populationTop^2)

# Jump times
times <- cumsum(rexp(100000,domRate))
times <- times[times<100]
tail(times,1)

# Storage
process <- matrix(0,nrow = 1+length(times),3); colnames(process) <- c("Time", "prey", "predator")
process[,1] <- c(0,times)
process[1,2:3] <- c(round(populationTop*0.75),round(populationTop*0.25))
totalJumps = 0

# Populate
for (i in 2:nrow(process)){
  alphaProb <- seasonality(process[i,1]) * alpha * process[i-1,2] * (process[i-1,2]<populationTop)/domRate
  betaProb <- seasonality(process[i,1]) * beta * process[i-1,2] * process[i-1,3] / domRate
  deltaProb <- seasonality(process[i,1]) * delta * process[i-1,2] * process[i-1,3] * (process[i-1,3]<populationTop)/ domRate
  gammaProb <- seasonality(process[i,1]) * gamma*process[i-1,3]/domRate
  threshold <- runif(1)
  
  if( threshold < alphaProb){
    process[i,2] <- process[i-1,2] + 1
    process[i,3] <- process[i-1,3]
    totalJumps <- totalJumps+1
  } else if( threshold < alphaProb + betaProb){
    process[i,2] <- process[i-1,2] - 1
    process[i,3] <- process[i-1,3]
    totalJumps <- totalJumps+1
  } else if( threshold < alphaProb + betaProb + deltaProb){
    process[i,3] <- process[i-1,3] + 1
    process[i,2] <- process[i-1,2]
    totalJumps <- totalJumps+1
  } else if( threshold < alphaProb + betaProb + deltaProb + gammaProb){
    process[i,3] <- process[i-1,3] - 1
    process[i,2] <- process[i-1,2]
    totalJumps <- totalJumps+1
  } else{
    process[i,2] <- process[i-1,2]
    process[i,3] <- process[i-1,3]
  }
}
plot(rep(process[,1],each=2)[-1],rep(process[,2],each=2)[-nrow(process)*2],col="black",type="l",ylim=c(0,populationTop)) # Prey
lines(rep(process[,1],each=2)[-1],rep(process[,3],each=2)[-nrow(process)*2],col="red",type="l") # Predator
print(totalJumps)


# Record observations
obsTimes <- seq(1,100,5)
observations <- cbind(obsTimes, t(sapply(obsTimes, function(i) process[tail(which(i>process[,1]),1),2:3])))
observations[,2] <- rnorm(length(obsTimes), mean = observations[,2], sd =populationTop/25)
observations[,3] <- rnorm(length(obsTimes), mean = observations[,3], sd =populationTop/25)
print(colMeans(observations[,2:3]))
points(observations[,1],observations[,2],pch=20,col="blue")
points(observations[,1],observations[,3],pch=20,col="green")

write.table(observations,file=paste("data_",populationTop,".txt",sep = ""),row.names = F,sep = " ",col.names = F)

