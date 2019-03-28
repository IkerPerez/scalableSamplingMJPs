
# First we create some data #
#############################

# Parameters
birth <- 30 #30 #10 #5 # 2 #1 # 0.5
death <- 2 * 0.01
populationTop <- 3000 #3000 # 1000 #500# 200 # 100 #50
seasonality <- function(t){
  1.5 + cos(2*pi/100*t)/2
}

# Dominating rate
domRate <- 2* death * populationTop + birth
  
# Jump times
times <- cumsum(rexp(100000,domRate))
times <- times[times<100]
tail(times,1)

# Storage
process <- matrix(0,nrow = 1+length(times),2); colnames(process) <- c("Time", "population")
process[,1] <- c(0,times)
process[1,2] <- floor(populationTop)
totalJumps = 0

# Populate
for (i in 2:nrow(process)){
  birthProb <- birth*(process[i-1,2]<populationTop)/domRate
  deathProb <- death * seasonality(process[i,1]) * process[i-1,2] /domRate
  threshold <- runif(1)
  if( threshold < birthProb){
    process[i,2] <- process[i-1,2] + 1
    totalJumps <- totalJumps+1
  } else if( threshold < birthProb + deathProb){
    process[i,2] <- process[i-1,2] - 1
    totalJumps <- totalJumps+1
  } else{
    process[i,2] <- process[i-1,2]
  }
}
plot(rep(process[,1],each=2)[-1],rep(process[,2],each=2)[-nrow(process)*2],col="black",type="l",ylim=c(min(process[,2]),populationTop))
print(totalJumps)

# Record observations
obsTimes <- seq(1,100,2)
observations <- cbind(obsTimes, sapply(obsTimes, function(i) process[tail(which(i>process[,1]),1),2]))
observations[,2] <- rnorm(nrow(observations),mean = observations[,2],sd = populationTop/50)
print(mean(observations[,2]))
points(observations[,1],observations[,2],pch=20,col="red")


# ODE solution for auxiliary variable approach
##############################################
library(deSolve)

bd.model <- function (t, x, params) {
  
  # states
  X <- x[1]

  # parameters
  death <- params[1]
  
  # dynamics
  dX <- birth*(X<populationTop) - death*seasonality(t)*X

  #return list
  list(c(dX))
}

grid <- seq(0,100,100/1000)
indexes <- sapply(observations[,1], function(z) tail(which(z >= grid),1))

error.bd <- function(params0){
  
  # Solve ODEs
  out <- as.data.frame(ode(y=c(X=populationTop),times=grid,bd.model,parms=params0,hmax=1/120))
  
  # Error term
  er <- sum((out$X[indexes]-observations[,2])^2)

  # return
  return(er)
}

fit0 <- optim(c(0.01),error.bd); fit0$par
params0 <- fit0$par
out <- as.data.frame(ode(y=c(X=populationTop),times=grid,bd.model,parms=params0,hmax=1/120))
lines(out$time,out$X)


# Save all necessary for inference and a plot
#save(process,observations,out,file = "toplot.Rdata")

write.table(out,file=paste("ode_",populationTop,".txt",sep = ""),row.names = F,sep = " ",col.names = F)
write.table(observations,file=paste("data_",populationTop,".txt",sep = ""),row.names = F,sep = " ",col.names = F)

write.table(out,file=paste("ode2_",populationTop,".txt",sep = ""),row.names = F,sep = " ",col.names = F)
write.table(observations,file=paste("data2_",populationTop,".txt",sep = ""),row.names = F,sep = " ",col.names = F)


