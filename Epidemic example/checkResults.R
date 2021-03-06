

# Read files
trace1 <- read.csv("mcmc_chain_1.csv")
trace2 <- read.csv("mcmc_chain_2.csv")
trace3 <- read.csv("mcmc_chain_3.csv")
trace4 <- read.csv("mcmc_chain_4.csv")
trace5 <- read.csv("mcmc_chain_5.csv")
trace6 <- read.csv("mcmc_chain_6.csv")
trace7 <- read.csv("mcmc_chain_7.csv")
trace8 <- read.csv("mcmc_chain_8.csv")

burnin <- round(nrow(trace1)*0.1)
thinning <- 1
total <- nrow(trace1)

library(coda)

N <- 50 # 1000 # 200 #100 # 50 
beta <- 2/N
gamma <- 1

# Infection Rate
plot(density(trace1[-c(1:burnin),1]), main="beta"); lines(density(trace2[-c(1:burnin),1]), col="red"); 
lines(density(trace3[-c(1:burnin),1]), col="blue"); lines(density(trace4[-c(1:burnin),1]), col="orange"); 
lines(density(trace5[-c(1:burnin),1]), col="green"); lines(density(trace6[-c(1:burnin),1]), col="yellow"); 
lines(density(trace7[-c(1:burnin),1]), col="brown"); lines(density(trace8[-c(1:burnin),1]), col="pink"); abline(v=beta,lty =2)

# Removal Rate
plot(density(trace1[-c(1:burnin),2]), main="gamma"); lines(density(trace2[-c(1:burnin),2]), col="red"); 
lines(density(trace3[-c(1:burnin),2]), col="blue"); lines(density(trace4[-c(1:burnin),2]), col="orange"); 
lines(density(trace5[-c(1:burnin),2]), col="green"); lines(density(trace6[-c(1:burnin),2]), col="yellow");
lines(density(trace7[-c(1:burnin),2]), col="brown"); lines(density(trace8[-c(1:burnin),2]), col="pink"); abline(v=gamma,lty =2)

# Initial infection
plot(density(trace1[-c(1:burnin),3]), main="i_0"); lines(density(trace2[-c(1:burnin),3]), col="red"); 
lines(density(trace3[-c(1:burnin),3]), col="blue"); lines(density(trace4[-c(1:burnin),3]), col="orange");

# Effective samples - across 16 simulations (repeat twice)
acf(trace1); acf(trace2); acf(trace3); acf(trace4)
acf(trace5); acf(trace6); acf(trace7); acf(trace8)


effsSamps1<- rbind(effectiveSize(as.mcmc(trace1)), effectiveSize(as.mcmc(trace2)), effectiveSize(as.mcmc(trace3)),  effectiveSize(as.mcmc(trace4)),
               effectiveSize(as.mcmc(trace5)), effectiveSize(as.mcmc(trace6)), effectiveSize(as.mcmc(trace7)),  effectiveSize(as.mcmc(trace8)))


