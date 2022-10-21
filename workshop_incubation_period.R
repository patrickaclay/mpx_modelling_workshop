# Predictive Modeling Workshop
# N = 6 cases

# Code adapted from Lessler et al. Times to key events in Zika virus infection and implications for blood donation: a systematic review - PMC (nih.gov)
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5096355/


library(dplyr)
library(ggplot2)
library(gridExtra)
library(rjags)
library(knitr)



# Creating incubation period dataset



id <- seq(1, 6) # 6 cases
EL <- rep(0, 6)



# Exposure window
ER <- c() # FILL IN USING TABLE



ER <- ER -0.0007



# Max days to symptom onset
SR <- c() # FILL IN USING TABLE



SL <- SR - 1
SR <- SR -0.0007



monkeypox <- data.frame(id, EL, ER, SL, SR)



# Visual summary of interval-censored data



dat <- monkeypox %>%
  mutate(ELnew = EL-ER,
         ERnew = ER-ER,
         SLnew = SL-ER,
         SRnew = SR-ER)



IP_plot <- ggplot(dat, aes(y=factor(id))) +
  geom_segment(aes(x=ELnew, xend=ERnew, yend=factor(id)),
               color="blue", size=2) +
  geom_segment(aes(x=SLnew, xend=SRnew, yend=factor(id)),
               size=2, color="red", alpha=.5) +
  ggtitle("Incubation Period data") +
  xlab(NULL) +
  ylab("id") +
  coord_cartesian(xlim = c(-40, 20)) +
  theme(axis.text.y = element_text(size=6)) +
  annotate("text", x=-35, y="900.1", label="A")



IP_plot



# The intervals observed for each indivdual.
# Blue intervals represent windows of possible exposure.
# Red intervals represent windows of possible time of symptom onset



# Bayesian MCMC Framework for Estimating Key Distributions



# the time of infection for each individual, Ei as drawn from a uniform prior
# defined by the earliest and latest possible times of exposure (ELi and ERi)



# length of the incubation period, Y is interval censored random variable following a lognormal distribution



# SLi is the earliest possible time of symptom onset for case i
# SRi is the latest time of symptom onset for case i



# Make a data object for JAGS
monkeypox_jags <- monkeypox



jags_data <- list(       
  ER=monkeypox_jags$ER,
  SL=monkeypox_jags$SL,
  SR=monkeypox_jags$SR)



# Let jags know the censoring situation
# 1. means interval censored
# 2. means event occured after known time



jags_data$IPisCensored=rep(1,6) # no missing



# Define variables to hold the length of the
# time to event for symptoms
jags_data$Y_S <- rep(NA, 6)



#set the initial values for time to event (i.e., Y_S)
IPyInit =  jags_data$SL
IPyInit[which(is.na(jags_data$SL)==T)]=0
IPyInit[which(IPyInit==0)]=0.0000000011



# Set the parameters we want to track
parameters <- c("lm","lsd", "E")



# Specify JAGS model inside R
model1_string <- "
model {
  
  for(i in 1:6){
  
    IPisCensored[i] ~ dinterval( Y_S[i] , IPcensorLimitVec[i,1:2])
    IPcensorLimitVec[i,1]<-max(0.000000001,(SL[i]-E[i]))
    IPcensorLimitVec[i,2]<-max(0.000000001,(SR[i]-E[i]))
    
    Y_S[i] ~ dlnorm(lm,  tau1)
    
    E[i] ~ dunif(0,ER[i])
  }
  
  lm~dnorm(0,0.001)
  tau1 <-1/lsd^2
  lsd~dunif(0,3)
  
}
"



DistributionFitLWW <-textConnection(model1_string)



set.seed(12345)
# Initialization function for jags
jags_inits <-  function() {
  rc <-list(E=rep(0.0000000011,6), #start with a fixed E to avoid bad starting points
            lm=runif(1,log(2),log(10)),
            lsd = runif(1,.1,log(3)),
            Y_S=IPyInit)
  print(rc)
  return(rc)
}



# Initialize JAGS model
jagsfit_LWW <- jags.model(DistributionFitLWW,
                          data=jags_data,
                          inits=jags_inits,
                          n.chains=3, quiet=F,
                          n.adapt=10000)



iters<- 1000000
thin <- 50
full_fit_LWW <- coda.samples(jagsfit_LWW, parameters, n.iter=iters, thin=thin, n.chains=3)



# Make all of the chains a single matrix, removing the burn-in
ABC1=as.matrix(full_fit_LWW[[1]][,])
ABC2=as.matrix(full_fit_LWW[[2]][,])
ABC3=as.matrix(full_fit_LWW[[3]][,])
ABC1=ABC1[5001:(iters/thin),]
ABC2=ABC2[5001:(iters/thin),]
ABC3=ABC3[5001:(iters/thin),]



# Recreate MCMC object for diagnostics
full_fit_LWW <- list(as.mcmc(ABC1), as.mcmc(ABC2), as.mcmc(ABC3))
chains_LWW <- rbind(ABC1,ABC2,ABC3)
colnames(chains_LWW) <- varnames(full_fit_LWW[[1]])
chains_LWW <- as.data.frame(chains_LWW)



# R_hat statistic for primary results
print(gelman.diag(full_fit_LWW))



# Estimates of key distributions



# Distribution of parameters and key quantiles of the incubation period for MPXV infection.
inc_fit_jags_LWW <-  quantile(exp(chains_LWW$lm+chains_LWW$lsd^2/2), prob=c(0.5,0.025,0.975))
inc_fit_jags_LWW <- rbind(inc_fit_jags_LWW,
                          exp(c(median(chains_LWW$lm), quantile(chains_LWW$lm,prob=c(0.025,0.975)))))
inc_fit_jags_LWW <- rbind(inc_fit_jags_LWW,
                          exp(c(mean(chains_LWW$lsd), quantile(chains_LWW$lsd,prob=c(0.025,0.975))))) # https://www.cabdirect.org/cabdirect/abstract/19512202981
inc_fit_jags_LWW <- rbind(inc_fit_jags_LWW,
                          quantile(sqrt((exp(chains_LWW$lsd^2)-1)*exp(2*chains_LWW$lm + chains_LWW$lsd^2)),
                                   prob=c(0.5,0.025,0.975))) # https://en.wikipedia.org/wiki/Log-normal_distribution



for (q in c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  tmp <- qlnorm(q, chains_LWW$lm, chains_LWW$lsd)
  inc_fit_jags_LWW <- rbind(inc_fit_jags_LWW,
                            c(mean(tmp), quantile(tmp, prob=c(0.025, 0.975))))
}
colnames(inc_fit_jags_LWW) <- c("est","CIlow","CIhigh")
rownames(inc_fit_jags_LWW) <- c("mean","median",
                                "dispersion",
                                "sd",
                                "p5","p25","p50","p75","p95")



kable(inc_fit_jags_LWW, format="markdown", digits=2)



# Summarize log normal distribution for table in manuscript
# log mean
mean(chains_LWW$lm)
quantile(chains_LWW$lm, prob = c(0.025, 0.975))



# log sd
mean(chains_LWW$lsd)
quantile(chains_LWW$lsd, prob = c(0.025, 0.975))



inc_curve <- NULL



for (q in seq(0,60,.1)) {
  tmp <- plnorm(q, chains_LWW$lm, chains_LWW$lsd)
  tmp <- quantile(tmp, prob=c(0.025, .5, 0.975))
  inc_curve <- rbind(inc_curve, c(q=q,
                                  plow=tmp[1],
                                  pmid=tmp[2],
                                  phigh=tmp[3]))
  
}
inc_curve <- as.data.frame(inc_curve)
colnames(inc_curve) <- c("q","incplow","incpmid","incphigh")



inc_plt <- ggplot(inc_curve, aes(x=q)) +
  geom_ribbon(aes(ymin=incplow, ymax=incphigh), fill="blue", alpha=.2) +
  geom_line(aes(y=incpmid), col="blue") +theme_bw()+
  scale_x_continuous(limits=c(0, 60), expand = c(0, 0)) +
  ylab("Proportion") + xlab("Incubation period (days)")



inc_plt
