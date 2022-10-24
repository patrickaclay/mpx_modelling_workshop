# Here, we model the transmission of monkeypox through a population of 
# Gay, Bisexual, And other Men who have sex with Men using a 
# Stochastic Exponential Random Graph Approach. 

# Written by Emily Pollock and Patrick Clay
# Last updated 10/14/2022

#load in libraries, set seed
library(here)
library(EpiModel)
library(EpiModelHIV)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(lhs)
library(kableExtra)
library(GGally)
library(tidyr)
set.seed(12345)

##############################

### Run Stergm ###

# load network structure
load(here("mpx_network_object.rda"))

# source modules for stergm
source(here("workshop_network_modules.R"))


# Set parameters 
probability.infection.per.sex.act <- 0.9 
length.of.infectious.period <- 60  #in days
initial.number.infected <- 10
sex.act.prob.main.partner <- 1 #chance per day of having sex with main partner (set at 100% for demo purposes)
sex.act.prob.casual.partner <- 1 #chance per day of having sex with casual parter (set at 100% for demo purposes)
sex.act.prob.onetime.partner <- 1 #chance per day of having sex with casual parter (100% by definition)

# load in parameters
params <- param_msm()
# sequence of modules to run model, set duration and number simulations, and how many cores to use
controls <- control_msm(nsteps = 700, nsims=1, ncores=1)
# initial conditions
inits <- init_msm()
# run simulation
sim <- netsim(est, params, inits, controls)
#you can look at the various model outputs with sim$epi$ and then looking at autofill options
#extract medians and interquartile ranges (note, this extracts values for if you increase number of simulations)
#we currently only run one simulation so this will just give you output of that single sim. 
run_results_median <- extract_med_iqr(sim)

##########################
#make figures
##########################

#get data ready to show proportion of new cases coming from each relationship type
inst_vec <- run_results_median$si.flow.onetime.med
pers_vec <- run_results_median$si.flow.casual.med
main_vec <- run_results_median$si.flow.main.med
partnership_data <- data.frame("Day"=run_results_median$Day, "One-Time" = inst_vec,
                        "Casual" = pers_vec, "Main" = main_vec)
partnership_data <- partnership_data %>% pivot_longer(names(partnership_data)[-1], names_to="Partnership Type")
cols = c("darkolivegreen3", "goldenrod1", "steelblue3")
partnership_data$`Partnership Type` <- factor(partnership_data$`Partnership Type`, levels=c("Main", "Casual", "One.Time"))


#Visualize where partnerships are coming from
partnership_data %>%
  ggplot(aes(x=Day, y=value, fill=`Partnership Type`)) +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 1), panel.grid.minor = element_line(size = 0.5)) +
  theme(axis.title.x = element_text(size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=14), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_col(width = 1) +
  scale_fill_manual(values=cols) +
  ylab("New Infections (Daily)") +
  xlab("Time since importation of infection (days)") +
  theme(strip.background = element_rect(color="black", fill="white", linetype="solid"))


#Visualize epi curve
ggplot(data = run_results_median, aes(x = Day, y = prev.med)) +
  theme_bw() +
  theme(panel.grid.major = element_line(size = 1), panel.grid.minor = element_line(size = 0.5)) +
  theme(axis.title.x = element_text(size=14),axis.text.x = element_text(size=10)) +  #Size of x-axis title
  theme(axis.title.y = element_text(size=14),axis.text.y = element_text(size=10)) +  #size of y-axis title
  theme(title =element_text(size=14), plot.title = element_text(hjust = 0.5))  +    #Size of legend Title
  geom_ribbon(aes(ymin = prev.iqr1, ymax = prev.iqr3), alpha = 0.2) +
  geom_line(color = "darkblue", size = 2) 


####################################################
##Once you have run the code, here are some extra tasks
####################################################

# 1. How do probability of infection and infectious period
# alter the relative importance of main vs. casual relationships? 
# The following code calculates the proportion of 
# transmissions due to each relationship.

proportion.due.to.onetime <- sum(run_results_median$si.flow.onetime.med)/sum(run_results_median$si.flow.med)
proportion.due.to.main <- sum(run_results_median$si.flow.main.med)/sum(run_results_median$si.flow.med)
proportion.due.to.casual <- sum(run_results_median$si.flow.casual.med)/sum(run_results_median$si.flow.med)

# Alter probability of infection and the duration of infection and see how those parameters impact the proportion
# of cases resulting from each relationship type.

###################################################

# 2. You can remove one-time, casual, or main relationships from the network
# by setting the corresponding sex.act.prob parameter to 0.
# Examine the prevalence curve (figure 2) if you run the epidemic with only
# one-time, only main, or only casual relationships. Can any relationship type alone
# sustain transmission? (start with a high probability of infection and a long infectious period to test this)

# Consider the following details:
# Individuals in the model can have 0, 1, or 2 casual relationships at a give time, 
# and these relationships last for about 160 days on average.
# Individuals in the model can have 0 or 1 casual relationships at a given time
# and these relationships last for about 400 days on average. 
# Based on this information, do you expect for pathogens to spread more easily on 
# the main or the casual network?
# What would need to be true for a pathogen to spread successfully on the 
# main partner network without any support from one-time or casual relationships?

