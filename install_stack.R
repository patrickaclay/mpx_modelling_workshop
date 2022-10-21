# Install Stack for MPX Workshop
# 24 October 2022

# If you don't have remotes, install it. 
install.packages("remotes")

# Install EpiModel from CRAN & EpiModel HIV Stack from Github
install.packages("EpiModel", dependencies = TRUE)
remotes::install_github("statnet/tergmLite")
remotes::install_github("statnet/EpiModelHPC") # you need this even if you aren't using HPC 
remotes::install_github("statnet/EpiModelHIV")

# Other packages use
install.packages("here")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("doParallel")
install.packages("foreach")
install.packages("tidyr")
install.packages("gridExtra")
install.packages("rjags")
install.packages("knitr")