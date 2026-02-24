Simu_stData.surface <- function(para =NULL, loc = NULL, W_ts = NULL, alphat = NULL, X = 0){
  time <- seq(0, 1,, para$Nt)
  TRUE.Y_ts <- sim.Wts <- NULL

  if(is.null(loc)){
    Y_ts <-  X_ts <- NULL
    load("./data/smoothed_surface.rds")

    for(r in 1:5){
      temp <- simData.DataBase[[r]]$reg.Site
      mean.Intensity <- simData.DataBase[[r]]$reg.Site$log.log.mean.Intensity#/100
      var.Intensity <- simData.DataBase[[r]]$reg.Site$log.log.var.Intensity#/100
      f <- function(t){
        z <- 0.3 * t^11 * (10 * (1 - t))^6 + 10 * (10 * t)^3 * (1 - t)^10
        z <- z/2e1
        return(z)
      }

      f_list <- list()
      mean.In <- (max(mean.Intensity) - mean.Intensity)/(max(mean.Intensity) - min(mean.Intensity))
      var.In <- (max(var.Intensity) - var.Intensity)/(max(var.Intensity) - min(var.Intensity))
      Ds <- fields::rdist(data.frame(temp$LAT_Y, temp$LON_X, temp$Distance))/1e3
      Dt <- fields::rdist(data.frame(temp$time.scale))
      f <- exp(-(mean.Intensity))*exp(-(var.Intensity))*exp(-Ds/20)*exp(-Dt/0.5)

      eig <- eigen(f)
      vectors <- eig$vectors
      K <- 50
      lambda <- runif(K, 1E-1, 1e0)
      for(k in 1:K){
        f_list[[k]] <- vectors[, k]
      }
      xi <- rnorm(K, mean = 0, sd = lambda)
      w <- Reduce(`+`, Map(`*`, f_list, xi))#/K


      simData.DataBase[[r]]$reg.Site$W_ts <- as.vector(w) - mean(w)
      if(r == 1){
        sim.Wts <- simData.DataBase[[r]]$reg.Site
      }else{
        sim.Wts <- rbind(sim.Wts, simData.DataBase[[r]]$reg.Site)
      }
    }
    n  <- length(unique(sim.Wts$Village_ID))
    Nt <- length(unique(sim.Wts$time))
    da <- sim.Wts %>% reshape2::dcast(time~Village_ID, value.var = "W_ts")

    W_ts <- as.matrix(da[, -1]) - mean(as.matrix(da[, -1]))
    # W_ts <-  W_ts - mean(W_ts)
    loc <- unique(sim.Wts[, c(1:7)]) %>% data.table::setorderv(c("Village_ID"))
  }else{
    id <- NULL
    for(r in 1:length(as.character(unique(loc$flag)))){
      Village_ID <- loc[loc$flag %in% unique(loc$flag)[r], ]$Village_ID

      id <- c(id, sample(Village_ID, para$n*length(Village_ID)/nrow(loc), replace = F))
      if(r == length(as.character(unique(loc$flag)))){
        id <- c(id, sample(Village_ID, para$n - length(id), replace = F))
      }
    }



    loc$Simu <- ifelse(loc$Village_ID %in%id, "Train", "Test")

    cat(unique(loc$Simu), "\n------\n")
    para$n <- nrow(loc)
    #######################################################################
    #######################################################################
    TRUE.Y_ts <- Y_ts <- matrix(0, nrow = para$Nt, ncol = para$n)
    colnames(Y_ts) <- loc$Village_ID
    px   <-  para$px
    X_ts <- array(0, dim = c(px, para$Nt, para$n),
                  dimnames = list(c(paste0("X", 1:px)),
                                  c(1:para$Nt),
                                  1:para$n))
    setDF(loc)
    mu <- matrix(0, nrow = para$Nt, ncol = para$n)

    setDF(X)

    for(t in 1:para$Nt){
      X_ts[1, t, ] <- 1
      X_ts[2, t, ] <- X[X$Year == (2010 + t), 2]
      X_ts[3, t, ] <- X[X$Year == (2010 + t), 3]
      X_ts[4, t, ] <- X[X$Year == (2010 + t), 4]
      X_ts[5, t, ] <- X[X$Year == (2010 + t), 5]
      X_ts[6, t, ] <- X[X$Year == (2010 + t), 6]
      X_ts[7, t, ] <- X[X$Year == (2010 + t), 7]
      X_ts[8, t, ] <- X[X$Year == (2010 + t), 8]
      X_ts[9, t, ] <- X[X$Year == (2010 + t), 9]
      X_ts[10, t, ] <- X[X$Year == (2010 + t), 10]
      X_ts[11, t, ] <- X[X$Year == (2010 + t), 11]
      mu[t, ]      <-  as.vector(t(X_ts[1:length(para$alpha), t, ]) %*% para$alpha)
      Y_ts[t, ]    <- TRUE.Y_ts[t, ]   <- mu[t, ] +  W_ts[t, ]  +
        rnorm(length(W_ts[t, ]), sd = sqrt(para$nugget))
    }
  }

  re <- list(
    Y_ts      = Y_ts,
    TRUE.Y_ts = TRUE.Y_ts,
    X_ts    = X_ts,
    W_ts    = W_ts,
    sim.Wts = sim.Wts,
    loc     = loc,
    locid   = "Village_ID"
  )
  return(re)
}


Simu_stData <- function(para =NULL, loc = NULL, W_ts = NULL, alphat = NULL, X = 0){


  time <- seq(0, 1,, para$Nt)
  TRUE.Y_ts <- sim.Wts <- NULL

  if(is.null(loc)){
    Y_ts <-  X_ts <- NULL
    load("./data/sim.Cov.Data.rds")
    for(r in 1:5){

      D <- chol(sim.Cov.Data[[r]]$reg.cov_matrix)
      w <-  t(D) %*% rnorm(nrow(sim.Cov.Data[[r]]$reg.cov_matrix))
      sim.Cov.Data[[r]]$reg.Site$W_ts <- w[, 1]# - mean(w[, 1])

      temp <- sim.Cov.Data[[r]]$reg.Site
      w <- as.matrix((temp %>% reshape2::dcast(time~Village_ID, value.var = "W_ts")))[, -1]
      sim.Cov.Data[[r]]$reg.Site$W_ts <- as.vector(w)
      if(r == 1){
        sim.Wts <- sim.Cov.Data[[r]]$reg.Site
      }else{
        sim.Wts <- rbind(sim.Wts, sim.Cov.Data[[r]]$reg.Site)
      }
    }
    n  <- length(unique(sim.Wts$Village_ID))
    Nt <- length(unique(sim.Wts$time))
    da <- sim.Wts %>% reshape2::dcast(time~Village_ID, value.var = "W_ts")

    W_ts <- as.matrix(da[, -1]) - mean(as.matrix(da[, -1]))
    # W_ts <-  W_ts - mean(W_ts)
    loc <- unique(sim.Wts[, c(1:7)]) %>% data.table::setorderv(c("Village_ID"))

  }else{
    id <- NULL#para$n
    for(r in 1:length(as.character(unique(loc$flag)))){
      Village_ID <- loc[loc$flag %in% unique(loc$flag)[r], ]$Village_ID
      # cat("...", length(Village_ID), "\n")
      id <- c(id, sample(Village_ID, para$n*length(Village_ID)/nrow(loc), replace = F))
      if(r == length(as.character(unique(loc$flag)))){
        id <- c(id, sample(Village_ID, para$n - length(id), replace = F))
      }
    }



    loc$Simu <- ifelse(loc$Village_ID %in%id, "Train", "Test")

    cat(unique(loc$Simu), "\n------\n")
    para$n <- nrow(loc)
    #######################################################################
    #######################################################################
    TRUE.Y_ts <- Y_ts <- matrix(0, nrow = para$Nt, ncol = para$n)
    colnames(Y_ts) <- loc$Village_ID
    px   <-  para$px
    X_ts <- array(0, dim = c(px, para$Nt, para$n),
                  dimnames = list(c(paste0("X", 1:px)),
                                  c(1:para$Nt),
                                  1:para$n))


    library(mvnfast)
    setDF(loc)
    mu <- matrix(0, nrow = para$Nt, ncol = para$n)

    setDF(X)

    for(t in 1:para$Nt){
      X_ts[1, t, ] <- 1
      X_ts[2, t, ] <- X[X$Year == (2010 + t), 2]
      X_ts[3, t, ] <- X[X$Year == (2010 + t), 3]
      X_ts[4, t, ] <- X[X$Year == (2010 + t), 4]
      X_ts[5, t, ] <- X[X$Year == (2010 + t), 5]
      X_ts[6, t, ] <- X[X$Year == (2010 + t), 6]
      X_ts[7, t, ] <- X[X$Year == (2010 + t), 7]
      X_ts[8, t, ] <- X[X$Year == (2010 + t), 8]
      X_ts[9, t, ] <- X[X$Year == (2010 + t), 9]
      X_ts[10, t, ] <- X[X$Year == (2010 + t), 10]
      X_ts[11, t, ] <- X[X$Year == (2010 + t), 11]
      mu[t, ]      <-  as.vector(t(X_ts[1:length(para$alpha), t, ]) %*% para$alpha)
      Y_ts[t, ]    <- TRUE.Y_ts[t, ]   <- mu[t, ] +  W_ts[t, ]  +
                      rnorm(length(W_ts[t, ]), sd = sqrt(para$nugget))
    }
  }

  re <- list(
    Y_ts      = Y_ts,
    TRUE.Y_ts = TRUE.Y_ts,
    X_ts    = X_ts,
    W_ts    = W_ts,
    sim.Wts = sim.Wts,
    loc     = loc,
    locid   = "Village_ID"
  )
  return(re)
}


data.collect <- function(simData, sim.para){
  simData.DataBase <- data.table(
    Village_ID = as.vector(t(matrix(simData$loc$Village_ID,
                                    nrow = nrow(simData$loc),
                                    ncol = sim.para$Nt))),
    flag = as.vector(t(matrix(simData$loc$flag,
                              nrow = nrow(simData$loc),
                              ncol = sim.para$Nt))),
    region = as.vector(t(matrix(simData$loc$flag,
                                nrow = nrow(simData$loc),
                                ncol = sim.para$Nt))),
    Study_Arm = as.vector(t(matrix(simData$loc$Study_Arm,
                                   nrow = nrow(simData$loc),
                                   ncol = sim.para$Nt))),
    LON = as.vector(t(matrix(simData$loc$LON,
                             nrow = nrow(simData$loc),
                             ncol = sim.para$Nt))),
    LAT = as.vector(t(matrix(simData$loc$LAT,
                             nrow = nrow(simData$loc),
                             ncol = sim.para$Nt))),
    time.scale = rep(seq(0, 1, length = sim.para$Nt), times = nrow(simData$loc)),
    time.index = rep(1:sim.para$Nt, times = nrow(simData$loc)),
    Year = rep(2010 + 1:sim.para$Nt, times = nrow(simData$loc)),
    Y_ts = as.vector(simData$Y_ts),
    TRUE.Y_ts = as.vector(simData$TRUE.Y_ts),
    W_ts = as.vector(simData$W_ts),
    Intercept = as.vector(simData$X_ts[1,,]),
    CWT_1 = as.vector(simData$X_ts[2,,]),
    CWT_2 = as.vector(simData$X_ts[3,,]),
    CWT_3 = as.vector(simData$X_ts[4,,]),
    CWT_4 = as.vector(simData$X_ts[5,,]),
    SBT_1 = as.vector(simData$X_ts[6,,]),
    SBT_2 = as.vector(simData$X_ts[7,,]),
    SBT_3 = as.vector(simData$X_ts[8,,]),
    SBT_4 = as.vector(simData$X_ts[9,,]),
    IEt.CWT = as.vector(simData$X_ts[10,,]),
    IEt.SBT = as.vector(simData$X_ts[11,,]),
    Simu = as.vector(t(matrix(simData$loc$Simu,
                              nrow = nrow(simData$loc),
                              ncol = sim.para$Nt)))
  )
  return(simData.DataBase)
}
