
# initialize_msm builds the initial network and assigns attributes to nodes
# param_msm enters all model parameters
# control_msm runs all modules in a particular order, determines number and length of runs
# simnet_msm controls the breaking and creating of sexual contact at each timestep
# infect_msm2 determines new infections at each timestep
# discord_edgelist  determine susceptible/infected pairs
# progress2 determines state transition for infected individuals
# prevalence_msm records prevalence and other epidemiological trackers



##################### Netsim Setup Functions #####################################

# Parameters ####
param_msm <- function(init.inf = initial.number.infected,                                # initial number of infected people 
                      init.hiv.age = c(0.044, 0.109, 0.154, 0.183), # prev rate among 18:24, 25:29, 30:34, 35:39
                      population.size = 10000,

                      act.rate.main = sex.act.prob.main.partner, 
                      act.rate.casual = sex.act.prob.casual.partner, 
                      act.rate.instant = sex.act.prob.onetime.partner, 

                      inf.prob = probability.infection.per.sex.act,       # probability of infection upon sexual contact
                      i.to.r.rate = 1/length.of.infectious.period,  # natural clearance rate
                      
                      ...) {
  
  ## Process parameters
  
  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))
  
  class(p) <- "param.net"
  return(p)
  
}


# Control Settings ####
control_msm <- function(simno = 1,
                        nsteps = 200,
                        start = 1,
                        nsims = 1,
                        ncores = 4,
                        cumulative.edgelist = FALSE,
                        truncate.el.cuml = 0,
                        initialize.FUN = initialize_msm,
                        progress.FUN = progress_msm,
                        infection.FUN = infect_msm2,
                        resim_nets.FUN = simnet_msm,
                        prev.FUN = prevalence_msm,
                        verbose.FUN = verbose.net,
                        module.order = NULL,
                        save.nwstats = FALSE,
                        save.other = c("el", "attr"),
                        tergmLite = TRUE,
                        tergmLite.track.duration = FALSE, # CPN2
                        set.control.ergm = control.simulate.formula(MCMC.burnin = 2e5),
                        set.control.stergm = control.simulate.network(),
                        verbose = TRUE,
                        skip.check = TRUE,
                        ...) {
  
  formal.args <- formals(sys.function())
  dot.args <- list(...)
  p <- get_args(formal.args, dot.args)
  
  p$skip.check <- TRUE
  p$save.transmat <- FALSE
  
  bi.mods <- grep(".FUN", names(formal.args), value = TRUE)
  bi.mods <- bi.mods[which(sapply(bi.mods, function(x) !is.null(eval(parse(text = x))),
                                  USE.NAMES = FALSE) == TRUE)]
  p$bi.mods <- bi.mods
  p$user.mods <- grep(".FUN", names(dot.args), value = TRUE)
  p[["f.names"]] <- c(p[["bi.mods"]], p[["user.mods"]])
  p$save.other <- c("attr", "temp", "el", "p", "mpx_degdist")
  
  p$save.network <- FALSE
  if (is.null(p$verbose.int)) {
    p$verbose.int <- 1
  }
  
  p <- set.control.class("control.net", p) # CPN2
  return(p)
}


##################### Critical Network Modules #####################################

#### Initialize Module #####
# gets called once at the beginning of each simulation to construct networks
# and master dat object 

initialize_msm <- function(x, param, init, control, s) { #So this is what sets up the network that is then changed by simnet
  
  dat <- create_dat_object(param, init, control) #Ah, this is why everything depends on init
  
  #### Network Setup ####
  # Initial network setup
  # Simulate each network on the basis of their fits
  # Add to dat object 
  
  dat[["nw"]] <- list()
  nnets <- 3
  for (i in 1:nnets) {
    dat[["nw"]][[i]] <- simulate( 
      x[[i]][["fit"]],
      basis = x[[i]][["fit"]][["newnetwork"]],
      dynamic = FALSE
    )
  }
  nw <- dat[["nw"]]
  
  # Pull Network parameters
  dat[["nwparam"]] <- list()
  for (i in 1:nnets) {
    dat[["nwparam"]][i] <- list(x[[i]][-which(names(x[[i]]) == "fit")])
  }
  
  # Convert to tergmLite method
  dat <- init_tergmLite(dat)
  
  #### Nodal Attributes Setup ####
  dat[["attr"]] <- param[["netstats"]][["attr"]]
  
  num <- network.size(nw[[1]])
  dat <- append_core_attr(dat, 1, num)
  
  # Pull in attributes on network. 
  # We pull from the one-time network because this has both deg.main and deg.pers attributes 
  # (see estimation script)
  nwattr.all <- names(nw[[3]][["val"]][[3]])
  nwattr.use <- nwattr.all[!nwattr.all %in% c("na", "vertex.names")]
  for (i in seq_along(nwattr.use)) {
    dat$attr[[nwattr.use[i]]] <- get.vertex.attribute(nw[[3]], nwattr.use[i])
  }
  
  # Add other attributes 
  
  # First we need some initial conditions parameters
  init.inf <- get_param(dat, "init.inf")
  
  
  # Generate status vector based on num init infected 
  # starting values s / i
  # initial infections among highest-activity groups unless init size is > than size of those groups 
  
  status <- rep("s", num)
  riskg <- get_attr(dat, "riskg")
  
  
  # initial infecteds
  risk.high <- which(status == "s" & (riskg == "5" | riskg == "6"))
  
  if (length(risk.high) > init.inf) {
    infected <- sample(risk.high, init.inf, replace=FALSE)
  }
  else {infected <- sample(which(status=="s"), init.inf, replace=FALSE)}
  
  status[infected] <- "i"
  
  
  
  # Pull sqrt age to get age 
  sqrt.age <- get_attr(dat, "sqrt.age")
  age <- sqrt.age^2
  

  
  # set attributes
  dat <- set_attr(dat, "status", status)       # initial status
  dat <- set_attr(dat, "age", sqrt.age^2)      # age attribute to accompany the sqrt.age attribute

  #### Other Setup ####
  dat[["stats"]] <- list()
  dat[["stats"]][["nwstats"]] <- list()
  dat[["temp"]] <- list()
  dat[["epi"]] <- list()
  dat[["mpx_degdist"]] <- list()
  
  dat <- set_epi(dat, "num", at = 1,  num)
  dat <- set_epi(dat, "cuml.infs", at = 1, init.inf)
  dat <- set_epi(dat, "prev", at = 1, init.inf)
  
  
  # Setup Partner List for all 3 networks 
  # (only gets updated if tracking turned on in control function)
  for (n_network in seq_len(3)) {
    dat <- update_cumulative_edgelist(dat, n_network)
  }
  
  # Network statistics
  if (dat[["control"]][["save.nwstats"]]) {
    for (i in seq_along(x)) {
      nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])
      nwstats <- summary(
        dat[["control"]][["nwstats.formulas"]][[i]],
        basis = nwL,
        term.options = dat[["control"]][["mcmc.control"]][[i]][["term.options"]],
        dynamic = i < 3
      )
      
      dat[["stats"]][["nwstats"]][[i]] <- matrix(
        nwstats, 
        nrow = 1, ncol = length(nwstats),
        dimnames = list(NULL, names(nwstats))
      )
      
      dat[["stats"]][["nwstats"]][[i]] <- 
        as.data.frame(dat[["stats"]][["nwstats"]][[i]])
    }
  }
  
  class(dat) <- "dat"
  return(dat)
  
}



#### Simnet -- Simulate Networks & Update Coefs based on pop size changes ###############
simnet_msm <- function(dat, at) {
  

  
  
  ## Grab parameters from dat object 
  cumulative.edgelist <- get_control(dat, "cumulative.edgelist") # are we tracking the cumulative edgelist (T/F)
  truncate.el.cuml <- get_control(dat, "truncate.el.cuml")       # how long in the past do we keep edgelist
  set.control.stergm <- get_control(dat, "set.control.stergm")   # specific control settings for network simulation
  set.control.ergm <- get_control(dat, "set.control.ergm")       # specific control settings for network simulation
  
  

  ## Main network
  for (i in 1:length(dat$el)) {    #I believe this is where loops through overlapping networks
    nwparam <- EpiModel::get_nwparam(dat, network = i)   #get parameters of this network
    isTERGM <- ifelse(nwparam$coef.diss$duration > 1, TRUE, FALSE) #set isTERGM to true is relationships non-instantaneous
    
    nwL <- networkLite(dat[["el"]][[i]], dat[["attr"]])
    
    
    if (get_control(dat, "tergmLite.track.duration") == TRUE) { #figure out how long relationships have lasted
      nwL %n% "time" <- dat[["nw"]][[i]] %n% "time"
      nwL %n% "lasttoggle" <- dat[["nw"]][[i]] %n% "lasttoggle"
    }
    
    if (isTERGM == TRUE) { #The following happens only if tergm (relationships last)
      dat[["nw"]][[i]] <- simulate( #forms, dissolves contacts. generic function. simulate distribution corresponding to fitted model object
        nwL, #what are the attributes of each node
        formation = nwparam[["formation"]],
        dissolution = nwparam[["coef.diss"]][["dissolution"]],
        coef.form = nwparam[["coef.form"]],
        coef.diss = nwparam[["coef.diss"]][["coef.adj"]],
        constraints = nwparam[["constraints"]],
        time.start = at - 1,
        time.slices = 1,
        time.offset = 1,
        control = set.control.stergm,
        output = "final"
      )
    } else {  #The following happens if tergm is false, e.g. instantaneous relationships
      dat[["nw"]][[i]] <- simulate( #forms contacts
        basis = nwL,
        object = nwparam[["formation"]],
        coef = nwparam[["coef.form"]],
        constraints = nwparam[["constraints"]],
        control = set.control.ergm,
        dynamic = FALSE,
        nsim = 1,
        output = "network"
      )
    }
    
    dat[["el"]][[i]] <- as.edgelist(dat[["nw"]][[i]])
    
    if (get_control(dat, "save.nwstats") == TRUE) {
      term.options <- if (isTERGM == TRUE) {
        set.control.stergm$term.options
      } else {
        set.control.ergm$term.options
      }
      dat$stats$nwstats[[i]] <- rbind(dat$stats$nwstats[[i]],
                                      summary(dat$control$nwstats.formulas[[i]],
                                              basis = nwL,
                                              term.options = term.options,
                                              dynamic = isTERGM))
    }
    
  }
  
  if (get_control(dat, "cumulative.edgelist") == TRUE) {
    for (n_network in seq_len(3)) {
      dat <- update_cumulative_edgelist(dat, n_network, truncate.el.cuml)
    }
  }
  
  # update main degree (nodal attribute based on network status)
  dat$attr$deg.main <- rep(0,length(dat$attr$deg.main))
  el <- get_edgelist(dat, 1)
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    for(j in 1:nrow(el)){
      dat$attr$deg.main[el[j,1]] <- 1
      dat$attr$deg.main[el[j,2]] <- 1
    }
  }
  
  # adjust pers degree
  dat$attr$deg.pers <- rep(0,length(dat$attr$deg.pers))
  el <- get_edgelist(dat, 2)
  if (nrow(el) > 0) {
    el <- el[sample(1:nrow(el)), , drop = FALSE]
    for(j in 1:nrow(el)){
      dat$attr$deg.pers[el[j,1]] <- dat$attr$deg.pers[el[j,1]] + 1
      if(dat$attr$deg.pers[el[j,1]] > 2){dat$attr$deg.pers[el[j,1]] <- 2}
      dat$attr$deg.pers[el[j,2]] <- dat$attr$deg.pers[el[j,2]] + 1
      if(dat$attr$deg.pers[el[j,2]] > 2){dat$attr$deg.pers[el[j,2]] <- 2}
    }
  }
  
  
  
  return(dat)
}






################### DISEASE RELATED MODULES ############################
# Infection ####
##Infection
infect_msm2 <- function(dat, at) {
  
  # model-specific discordant edgelist function
  discord_edgelist_mpx <- function (dat, at, network){ 
    status <- get_attr(dat, "status")
    active <- get_attr(dat, "active")
    el <- get_edgelist(dat, network)
    del <- NULL
    if (nrow(el) > 0) {
      el <- el[sample(1:nrow(el)), , drop = FALSE]
      stat <- matrix(status[el], ncol = 2)
      isInf <- matrix(stat %in% c("i", "im", "a"), ncol = 2)
      isSus <- matrix(stat %in% c("s", "v"), ncol = 2)
      SIpairs <- el[isSus[, 1] * isInf[, 2] == 1, , drop = FALSE]
      ISpairs <- el[isSus[, 2] * isInf[, 1] == 1, , drop = FALSE]
      pairs <- rbind(SIpairs, ISpairs[, 2:1])
      if (nrow(pairs) > 0) {
        sus <- pairs[, 1]
        inf <- pairs[, 2]
        del <- data.frame(at, sus, inf)
        keep <- rowSums(matrix(c(active[del$sus], active[del$inf]), 
                               ncol = 2)) == 2
        del <- del[keep, ]
        if (nrow(del) < 1) {
          del <- NULL
        }
      }
    }
    return(del)
  }
  
  
  #### Setup ####
  
  # Get attributes and parameters 
  active    <- get_attr(dat, "active")
  status    <- get_attr(dat, "status")
  riskgroup <- get_attr(dat, "riskg")
  
  inf.prob         <- get_param(dat, "inf.prob")
  act.rate.main    <- get_param(dat, "act.rate.main")
  act.rate.casual  <- get_param(dat, "act.rate.casual")
  act.rate.instant <- get_param(dat, "act.rate.instant")

  # Set up trackers 
  nInf <- 0
  idsInf <- NULL
  
  # New infections in each network
  nInfsMain <- 0
  nInfsPers <- 0
  nInfsInst   <- 0
  
  #### Transmissions ####
  # Loop through discordant edgelists in each network
  for(nw.transmit in 1:3){ 
    
    if(nw.transmit == 1){act.rate <- act.rate.main}
    if(nw.transmit == 2){act.rate <- act.rate.casual}
    if(nw.transmit == 3){act.rate <- act.rate.instant}
    
    # get discordant edgelist
    del <- discord_edgelist_mpx(dat, at, network = nw.transmit)
    
    if (!(is.null(del))) {
      
      # add status column for sus 
      del$statusSus <- status[del$sus]
      
      # set inf prob
      del$transProb <- inf.prob

      # act rate 
      del$actRate <- act.rate
      
      # final transmission probability 
      del$finalProb <- 1 - (1 - del$transProb)^del$actRate
      
      # filter down to which pairs transmitted infection
      transmit <- rbinom(nrow(del), 1, del$finalProb)
      del <- del[which(transmit == 1), ]
      
      idsNewInf <- unique(del$sus)
      idsOldInf <- unique(del$inf)
      
      if (length(idsNewInf) > 0) {
        
        status[idsNewInf]  <- "i"

        
        # in case any susceptible ids get infected by more than one person at this time step
        # we pick the first time they show up in the edgelist
        if (length(idsNewInf) < length(idsOldInf)) {
          
          tab <- table(del$sus)
          ids <- as.numeric(names(which(tab>1)))
          
          for (i in 1:length(ids)){
            #if(length(ids)>1) {          browser()}
            d <- del[del$sus %in% ids[i],]
            if(i==1){ 
              first <- d[1,]
            } else {
              first <- rbind(first, d[1,])
            }
          }
          
          single_del <- del[-c(which(del$sus %in% ids)),]
          newdel <- rbind(single_del, first)
          
          idsNewInf <- newdel$sus
          idsOldInf <- newdel$inf
          
        }
        
        
        
        if(nw.transmit == 1){nInfsMain <- length(idsNewInf)}
        if(nw.transmit == 2){nInfsPers <- length(idsNewInf)}
        if(nw.transmit == 3){nInfsInst <- length(idsNewInf)}
        
        idsInf <- c(idsInf, idsNewInf)
        
        nInf <- nInf + length(idsNewInf)

      }
    } 
  }
  

  #update attributes and trackers
  
  # nodal attributes 
  dat <- set_attr(dat, "status", status)

  # update epi trackers 
  dat <- set_epi(dat, "si.flow", at, nInf)
  
  prev.infs <- get_epi(dat, "cuml.infs", at=at-1)
  dat <- set_epi(dat, "cuml.infs", at, prev.infs + nInf)
  
  dat <- set_epi(dat, "si.flow.main", at, nInfsMain)
  dat <- set_epi(dat, "si.flow.casual", at, nInfsPers)
  dat <- set_epi(dat, "si.flow.onetime", at, nInfsInst)
  
  return(dat)
}


# Disease Progression ####
progress_msm <- function(dat, at) {
  
  # i  - symptomatic and infectious
  # r  - recovered and immune 
  
  ## Get Params & Attributes 
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")

    

  i.to.r.rate <-     get_param(dat, "i.to.r.rate")

  

    ## Natural Recovery among i
    nRec       <- 0
    idsEligRec <- which(active == 1 & status == "i")
    nEligRec   <- length(idsEligRec)
    
    if (nEligRec > 0) {
      vecRec <- which(rbinom(nEligRec, 1, i.to.r.rate) == 1) # vector of who recovers 
      if (length(vecRec) > 0) {
        idsRec <- idsEligRec[vecRec]
        nRec   <- length(idsRec)
        status[idsRec] <- "r"
      }
    }
    

    
  
  

  
  
  # Update attributes
  dat <- set_attr(dat, "status", status)

  # Update epidemiology trackers 

  dat <- set_epi(dat, "ir.flow", at, nRec)
  
  return(dat)
}


# Track Prevalence & Other Metrics ####
prevalence_msm <- function(dat, at) {
  
  active <- get_attr(dat, "active")
  status <- get_attr(dat, "status")
  
  nsteps <- dat$control$nsteps
  rNA <- rep(NA, nsteps)
  
  dat <- set_epi(dat, "num",   at, sum(active == 1))
  dat <- set_epi(dat, "num.i", at, sum(active == 1 & status == "i", na.rm = TRUE))
  dat <- set_epi(dat, "num.r", at, sum(active == 1 & status == "r", na.rm = TRUE))
  dat <- set_epi(dat, "num.s", at, sum(active == 1 & status == "s", na.rm = TRUE))
  dat <- set_epi(dat, "prev",  at, sum(active == 1 & status=="e") + 
                   sum(active == 1 & status=="a") + 
                   sum(active == 1 & status=="i"))
  
  
 
  
  return(dat)
}



# create function to extract median and IQR at each time step for all outputs
extract_med_iqr <- function(sim){
  epi <- sim$epi
  timesteps <- sim$control$nsteps
  vars <- names(epi)
  x <- NULL
  d <- matrix(NA, nrow=timesteps, ncol=3)
  
  for (i in vars) {
    #browser()
    for (time in 1:timesteps){
      t <- summary(as.numeric(epi[[i]][time,]))
      vals <- cbind(t[3], t[2], t[5])
      d[time,] <- vals
    }
    #browser()
    colnames(d) <- c(paste0(i, ".med"), paste0(i, ".iqr1"), paste0(i, ".iqr3"))
    x <- cbind(x, d)
    d <- matrix(NA, nrow=timesteps, ncol=3)
  }
  
  x <- x[-1,]
  x <- as.data.frame(x)
  x$Day <- 1:nrow(x)
  
  return(x)
}


