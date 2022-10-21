# Monkeypox Workshop 
# 24 October 2022
# Network Estimation
# Emily Pollock 

# Libraries ----------------------------------------------------------
library(EpiModel)
library(here)

# Target Statistics & Nodal Attributes -------------------------------

## Population Size
num <- 10000

# Overall Main/Pers Distribution (that many targets are based on)
deg <- matrix(c(0.471, 0.167, 0.074,
                 0.22, 0.047, 0.021), byrow = TRUE, nrow = 2) 

## Main Network Targets 
edges.main       <- 1440 # (# of people in a main relationships / 2)
pers.deg1        <- 470  # (# of people who have 1 main and 1 casual relationship)
pers.deg2        <- 210  # (# of people who have 1 main and 2 casual relationships)
absdiff.age.main <- 670  # (sum of sqrt of age difference between nodes in each relationship)

main.len <- 407 # (mean length of main partnerships in days)

## Casual Network Targets 
edges.cas       <- 2020
main.deg1       <- 890
conc            <- 950
absdiff.age.cas <- 1180

cas.len <- 166

## Instantaneous Network Targets
edges.inst <- 117

main1.pers0 <- 43
main0.pers1 <- 50
main1.pers1 <- 8.9
main0.pers2 <- 22
main1.pers2 <- 4

risk.1 <- 0
#risk.2 
risk.3 <- 10.3
risk.4 <- 19.3
risk.5 <- 60
risk.6 <- 140

absdiff.age.inst <- 64

inst.len <- 1


# Nodal Attributes 
## I, R, V role frequencies
role.prob <- c(0.242, 0.321, 0.437)

## Ages
ages <- seq(18, 39, 1)

# Generate Empty Networks w/ Initial Attributes --------------------------
nw <- network::network.initialize(num, directed = FALSE)

role.class <- sample(apportion_lr(num, c("I", "R", "V"), role.prob))
age <- seq(min(ages), max(ages) + 1, 1 / 365)
sqrt.age <- sqrt(age)
riskg <- sample(apportion_lr(num, 1:6, c(rep(0.19, 5), 0.05)))
deg.pers <- sample(apportion_lr(num, 0:2, colSums(deg)))

attr.names <- c("role.class", "age", "sqrt.age", "riskg", "deg.pers")
attr.values <- list(role.class, age, sqrt.age, riskg, deg.pers)

nw <- network::set.vertex.attribute(nw, attr.names, attr.values)

# let's view the network - it should have several vertex attributes but no edge attributes
nw 


# Estimation Workflow --------------------------------------------------
## Main Network 
### Formulas
formation.m <- ~edges +
  nodefactor("deg.pers") +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

diss.m <- dissolution_coefs(~offset(edges), duration = main.len, d.rate = 0)

targets.m <- c(edges.main, pers.deg1, pers.deg2, absdiff.age.main)

### Fit model
fit.m <- netest(nw,
                formation = formation.m,
                coef.form = c(-Inf, -Inf),
                target.stats = targets.m,
                coef.diss = diss.m,
                constraints = ~bd(maxout = 1)+sparse,
                set.control.ergm = control.ergm(MCMLE.maxit = 250))


# look at the fit object
summary(fit.m)

## Persistent Network
# Assign main degree based on previously fitted network
deg.main <- get_degree(fit.m$fit$newnetwork)
nw %v% "deg.main" <- deg.main

# Formulas
formation.p <- ~edges +
  nodefactor("deg.main") +
  concurrent +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

diss.p <- dissolution_coefs(~offset(edges), duration = cas.len, d.rate = 0)

targets.p <- c(edges.cas, main.deg1, conc, absdiff.age.cas)

# Fit model
fit.p <- netest(nw,
                formation = formation.p,
                coef.form = c(-Inf, -Inf),
                target.stats = targets.p,
                coef.diss = diss.p,
                constraints = ~bd(maxout = 2)+sparse,
                set.control.ergm = control.ergm(MCMLE.maxit = 250))





## Inst Network
### Update Pers Degree based on above fit
deg.pers <- get_degree(fit.p$fit$newnetwork)
nw %v% "deg.pers" <- deg.pers

### Formulas
formation.i <- ~edges +
  nodefactor(c("deg.main", "deg.pers")) +
  nodefactor("riskg", levels = -2) +
  absdiff("sqrt.age") +
  offset(nodematch("role.class", diff = TRUE, levels = 1:2))

targets.i <- c(edges.inst, 
               main1.pers0, main0.pers1, main1.pers1, main0.pers2, main1.pers2, 
               risk.1, risk.3, risk.4, risk.5, risk.6,
               absdiff.age.inst)

### Fit model
fit.i <- netest(nw,
                formation = formation.i,
                target.stats = targets.i,
                coef.form = c(-Inf, -Inf),
                coef.diss = dissolution_coefs(~offset(edges), 1),
                constraints = ~sparse,
                set.control.ergm = control.ergm(MPLE.max.dyad.types = 1e7,
                                                MCMLE.maxit = 250))

# Combine Networks into List
est <- list(fit.m, fit.p, fit.i) 
save(est, file=here("mpx_network_object.rda"))

# Diagnostics ----------------------------------------------------------

m.dx <- netdx(fit.m, dynamic=TRUE, nsims=10, nsteps=500)
p.dx <- netdx(fit.p, dynamic=TRUE, nsims=10, nsteps=500)
i.dx <- netdx(fit.i, dynamic=FALSE, nsims=1000)

