#' Run VB with EnKF updates
#'
#' @param data Input training data object.
#' @param test Optional test data object.
#' @param prior Prior specification list.
#' @param Para.List Initial parameter list.
#' @param center.y Logical; whether to center the response.
#' @param scale.y Logical; whether to scale the response.
#' @param Object Grouping variable name.
#' @param verbose Logical; whether to print progress messages.
#' @param verbose.VB Logical; whether to print VB-specific progress messages.
#' @param transf.Response Response transformation option.
#' @param Ne Number of ensemble members.
#' @param save.Predict Logical; whether to save predictions.
#' @param cs Spatial tapering range.
#' @param ct Temporal lag parameter.
#' @param n.cores Number of CPU cores to use.
#' @param itMax Maximum number of iterations.
#' @param itMin Minimum number of iterations.
#' @param seed Random seed.
#' @param tol.real Convergence tolerance.
#' @param plot Logical; whether to create plots.
#' @param positive Logical; whether to enforce positivity.
#' @param MCMC Logical; whether to use the MCMC routine.
.VB.EnKF <- function (data,
                      test            = NULL,
                      prior           = NULL,
                      Para.List       = NULL,
                      center.y        = FALSE,
                      scale.y         = FALSE,
                      Object          = "ALL",
                      verbose         = TRUE,
                      verbose.VB      = FALSE,
                      transf.Response = c("normal"),
                      Ne              = 100,
                      save.Predict    = FALSE,
                      cs              = 1e16,
                      ct              = 0,
                      n.cores         = 1,
                      itMax           = 100,
                      itMin           = 10,
                      seed            = 1234,
                      tol.real        = 0.001,
                      plot            = TRUE,
                      positive        = TRUE,
                      MCMC            = FALSE)
{
  compute_good_of_fitting <- function(obs, eta.samples, sigma.info) {
    obs.vec <- as.vector(obs)
    keep <- which(!is.na(obs.vec))
    obs.vec <- obs.vec[keep]
    eta.mat <- matrix(eta.samples,
                      nrow = length(as.vector(obs)),
                      ncol = dim(eta.samples)[3])[keep, , drop = FALSE]

    nsamp <- ncol(eta.mat)
    sigma.mean <- sigma.info$mu.sigma.sq
    sigma.draws <- rep(sigma.mean, nsamp)
    if (!is.null(sigma.info$a) && !is.null(sigma.info$b)) {
      sigma.draws <- 1 / rgamma(nsamp, shape = sigma.info$a, rate = sigma.info$b)
    }

    loglik.mat <- sapply(seq_len(nsamp), function(s) {
      dnorm(obs.vec, mean = eta.mat[, s], sd = sqrt(sigma.draws[s]), log = TRUE)
    })
    if (!is.matrix(loglik.mat)) {
      loglik.mat <- matrix(loglik.mat, ncol = 1)
    }

    eta.mean <- rowMeans(eta.mat)
    loglik.mean <- dnorm(obs.vec, mean = eta.mean, sd = sqrt(sigma.mean), log = TRUE)
    deviance.mean <- -2 * sum(loglik.mean)
    deviance.samples <- -2 * colSums(loglik.mat)
    mean.deviance <- mean(deviance.samples)
    p.d <- mean.deviance - deviance.mean
    dic <- mean.deviance + p.d

    lppd <- sum(log(rowMeans(exp(loglik.mat))))
    p.waic <- sum(apply(loglik.mat, 1, var))
    waic <- -2 * (lppd - p.waic)

    list(
      DIC = dic,
      p.eff.DIC = p.d,
      mean.deviance = mean.deviance,
      deviance.of.mean = deviance.mean,
      WAIC = waic,
      p.eff.WAIC = p.waic,
      lppd = lppd,
      loglik.matrix = loglik.mat
    )
  }
  compute_current_fitted_random_effect_draw <- function() {
    current_draw <- vector("list", PG)
    names(current_draw) <- names(G.Mat$data.G.Mat)

    for (Pg in seq_len(PG)) {
      region_name <- names(G.Mat$data.G.Mat)[Pg]
      current_draw[[region_name]] <- matrix(
        0,
        nrow = Nt,
        ncol = n,
        dimnames = list(time, ID)
      )

      if (Pg == 1) {
        if (is.null(ini.IEnKF.Input$train.Z.tsv) & is.null(ini.IEnKF.Input$train.sZ.tsv)) {
          temp <- sapply(seq_len(Nt), function(t) {
            G.Mat$data.G.Mat[[region_name]][t, , ] %*%
              slope.stRF.ens[[region_name]]$Vt.mean[t + 1, ]
          }, simplify = "array")
        } else if (!is.null(ini.IEnKF.Input$train.Z.tsv) & is.null(ini.IEnKF.Input$train.sZ.tsv)) {
          temp <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(
              G.Mat$data.G.Mat[[region_name]][t, , ],
              ini.IEnKF.Input$train.Z.tsv[[t]]
            ) %*% slope.stRF.ens[[region_name]]$Vt.mean[t + 1, ]
          }, simplify = "array")
        } else if (is.null(ini.IEnKF.Input$train.Z.tsv) & !is.null(ini.IEnKF.Input$train.sZ.tsv)) {
          temp <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(
              G.Mat$data.G.Mat[[region_name]][t, , ],
              ini.IEnKF.Input$train.sZ.tsv[[t]]
            ) %*% slope.stRF.ens[[region_name]]$Vt.mean[t + 1, ]
          }, simplify = "array")
        } else {
          temp <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(
              G.Mat$data.G.Mat[[region_name]][t, , ],
              ini.IEnKF.Input$train.Z.tsv[[t]],
              ini.IEnKF.Input$train.sZ.tsv[[t]]
            ) %*% slope.stRF.ens[[region_name]]$Vt.mean[t + 1, ]
          }, simplify = "array")
        }
      } else {
        temp <- sapply(seq_len(Nt), function(t) {
          spam::as.spam(G.Mat$data.G.Mat[[region_name]][t, , ]) %*%
            slope.stRF.ens[[region_name]]$Vt.mean[t + 1, ]
        }, simplify = "array")
      }

      for (t in seq_len(Nt)) {
        current_draw[[region_name]][t, ] <- temp[, , t]
      }
    }

    current_draw
  }
  {
    {
      run_time <- paste0(if_else(month(Sys.Date()) > 9
                                 , as.character(month(Sys.Date()))
                                 , paste0("0",  as.character(month(Sys.Date()))))
                         , "_", if_else(day(Sys.Date()) > 9
                                        , as.character(day(Sys.Date()))
                                        , paste0("0",  as.character(day(Sys.Date()))))
                         , "_", if_else(hour(Sys.time()) > 9
                                        , as.character(hour(Sys.time()))
                                        , paste0("0",  as.character(hour(Sys.time())))))
      digist.num <- 5
      if (!is.null(seed)) {
        set.seed(seed)
      }
      Py        <- length(data)
      true.Y_ts <- list()

      for(py in 1:Py){
        true.Y_ts[[py]] <- data[[py]]$Y_ts
        if (transf.Response %in% c("SQRT", "sr", "sqrt", "Sqrt")) {
          data[[py]]$Y_ts <- (data[[py]]$Y_ts)^(1/2)
        }
        if (transf.Response %in% c("log", "Log", "LOG")) {
          data[[py]]$Y_ts <- log(data[[py]]$Y_ts)
        }
        if (transf.Response %in% c("cubic")) {
          data[[py]]$Y_ts <- (data[[py]]$Y_ts)^(1/3)
        }
        if (transf.Response %in% c("loglog")) {
          # data[[py]]$Y_ts <- -log(5 - log(data[[py]]$Y_ts))

          data[[py]]$Y_ts <- log(-log(1 - data[[py]]$Y_ts))
        }

      }
      mean.y <-  rep(0, Py)
      sd.y   <- rep(1, Py)
      for(py in 1:Py){
        if(center.y&scale.y){
          mean.y[py]      <- mean(data[[py]]$Y_ts, na.rm = TRUE)
          sd.y[py]        <- sd(data[[py]]$Y_ts, na.rm = TRUE)
          data[[py]]$Y_ts <- (data[[py]]$Y_ts- mean.y[py])/sd.y[py]
        }else if(!center.y&scale.y){
          sd.y[py]        <- sd(data[[py]]$Y_ts, na.rm = TRUE)
          data[[py]]$Y_ts <- data[[py]]$Y_ts/sd.y[py]
        }else if(center.y&!scale.y){
          mean.y[py]      <- mean(data[[py]]$Y_ts, na.rm = TRUE)
          data[[py]]$Y_ts <- data[[py]]$Y_ts - mean.y[py]
        }
      }

      Miss.ini <- c(0.65, 0.8, 0.9)#c(0.65, 0.8, 0.9)
      # Initialize missing data ----
      Fill.Res.Miss.Values <- test.Res.miss.index <- train.Res.miss.index <- list()
      for(py in 1:Py){
        train.Res.miss.index[[py]] <- which(is.na(data[[py]]$Y_ts), arr.ind = TRUE)

        if (nrow(train.Res.miss.index[[py]]) > 0) {
          cat("There are some missing data ...\n")
          cat("Missing for training datasets are imputed by using median at the initial step.\n")
          if(!is.null(data[[py]]$ini.Miss_ts[1])){
            data[[py]]$Y_ts[train.Res.miss.index[[py]][, 1], train.Res.miss.index[[py]][, 2]] <-
              data[[py]]$ini.Miss_ts[train.Res.miss.index[[py]][, 1], train.Res.miss.index[[py]][, 2]] +
              rnorm(nrow(train.Res.miss.index[[py]]), 0, 1E-3)
          }else{
            for(i in 1:nrow( train.Res.miss.index[[py]])){
              t     <-  train.Res.miss.index[[py]][i, 1]
              s.ind <-  train.Res.miss.index[[py]][i, 2]
              if(t %in% c(2)){
                ini.miss.value <- quantile(data[[py]]$Y_ts[t, ], probs = Miss.ini[1], na.rm = TRUE)#max(data[[py]]$Y_ts, na.rm = TRUE)
              }
              if(t %in% c(3)){
                ini.miss.value <- quantile(data[[py]]$Y_ts[t, ], probs = Miss.ini[2], na.rm = TRUE)#max(data[[py]]$Y_ts, na.rm = TRUE)
              }
              if(t %in% c(4)){
                ini.miss.value <- quantile(data[[py]]$Y_ts[t, ], probs = Miss.ini[3], na.rm = TRUE)
              }

              data[[py]]$Y_ts[t, s.ind] <-  ifelse(is.na(data[[py]]$Y_ts[t, s.ind]),
                                                   ini.miss.value,
                                                   data[[py]]$Y_ts[t, s.ind])
            }

          }

        }


        if(!is.null(test)){
          test.Res.miss.index[[py]] <- which(is.na(test[[py]]$Y_ts), arr.ind = TRUE)
          if (nrow(test.Res.miss.index[[py]]) > 0) {
            cat("Missing for test datasets are imputed by using median at the initial step.\n")
            test[[py]]$Y_ts[test.Res.miss.index[[py]][, 1], test.Res.miss.index[[py]][, 2]] <-
              median(test[[py]]$Y_ts, na.rm = TRUE)
          }
        }
      }

      data.trans.tsv <- tranform.list.to.matrix(data)
      if(!is.null(test)){
        test.trans.tsv <- tranform.list.to.matrix(test)
      }


      n <- sPg <- 0
      Nt        <- data[[1]]$Nt
      Px        <- data[[1]]$Px
      # Pg        <- data[[1]]$Pg# + 1
      ID        <- ID.test <- NULL
      time      <- dimnames(data[[1]]$Y_ts)[[1]]
      time.test <- c(dimnames(test[[1]]$Y_ts)[[1]])


      for(py in 1:Py){
        n   <- n + data[[py]]$n
        sPg <- sPg + data[[py]]$sPg
        ID  <- c(ID, dimnames(data[[py]]$Y_ts)[[2]])
        if(!is.null(test)){
          ID.test <- c(ID.test, dimnames(test[[py]]$Y_ts)[[2]])
        }
      }

      yts.minus.xts <- data.trans.tsv$tsv.y
      if(!is.null(data[[1]]$X_ts)){
        if(is.na(Para.List[[names(data)[[1]]]]$beta$mu.beta[1])){
          public.beta.name   <- dimnames(data[[1]]$X_ts)[[1]]
          X_CuSum            <- XYXi_CuSum <- 0
          E_inverse_sigma.sq <- NULL
          for(py in 1:Py){
            E_inverse_sigma.sq <- c(E_inverse_sigma.sq, rep(1/Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq, data[[py]]$n))

          }

          for (t in 1:Nt) {
            X_ts <- matrix(as.matrix(data.trans.tsv$tsv.x[[t]]), nrow = dim(data[[1]]$X_ts)[[1]], ncol = n)
            X_CuSum    <- X_CuSum    + X_ts %*% (t(X_ts) *E_inverse_sigma.sq)
            XYXi_CuSum <- XYXi_CuSum + X_ts %*% (data.trans.tsv$tsv.y[t,]*E_inverse_sigma.sq)
          }
          post_betaX_sigma.sq <- tryCatch(
            solve(X_CuSum),
            error = function(e) ginv(X_CuSum)
          )

          eta                 <- XYXi_CuSum
          post_betaX_mu       <- post_betaX_sigma.sq %*% eta %>% as.vector()

          for(py in 1:Py){
            Para.List[[names(data)[[py]]]]$beta$mu.beta           <- matrix(NA, nrow = length(public.beta.name), ncol = 1)
            rownames(Para.List[[names(data)[[py]]]]$beta$mu.beta) <- public.beta.name
            colnames(Para.List[[names(data)[[py]]]]$beta$mu.beta) <- "Sharing coffficients"
            Para.List[[names(data)[[py]]]]$beta$mu.beta[, 1]      <- post_betaX_mu
            Para.List[[names(data)[[py]]]]$beta$sigma.sq          <- post_betaX_sigma.sq
          }

          yts.minus.xts <- yts.minus.xts - Xts.beta.fun(tsv.x = data.trans.tsv$tsv.x,
                                                        alpha = Para.List[[names(data)[[1]]]]$beta$mu.beta,
                                                        self  = FALSE)
        }
      }
      print(Para.List[[names(data)[[1]]]]$beta$mu.beta)



      fitting.data.frame <- list()
      for(py in 1:Py){
        fitting.data <- data.table(Y_ts       = as.vector(data[[py]]$Y_ts),
                                   DATE_TIME  = rep(dimnames(data[[py]]$Y_ts)[[1]], times = data[[py]]$n),
                                   time.index = rep(1:length(dimnames(data[[py]]$Y_ts)[[1]]), times = data[[py]]$n),
                                   LON        = as.vector(t(matrix(data[[py]]$spCoordinates$LON, nrow = data[[py]]$n, ncol = data[[py]]$Nt))),
                                   LAT        = as.vector(t(matrix(data[[py]]$spCoordinates$LAT, nrow = data[[py]]$n, ncol = data[[py]]$Nt))))
        if(!is.null(data[[py]]$X_ts)){
          for(px in 1:dim(data[[py]]$X_ts)[1]){
            fitting.data <- cbind(fitting.data, as.vector(data[[py]]$X_ts[px,,]))
          }
        }

        colnames(fitting.data) <- c(paste0(names(data)[py], ".Yts"), "DATE_TIME", "time.index", "LON", "LAT", dimnames(data[[py]]$X_ts)[[1]])

        fitting.data.frame[[names(data)[py]]] <- fitting.data
      }


      if(!is.null(test)) {
        pred.data.frame <- list()
        for(py in 1:Py){
          pred.data <- data.table(Y_ts       = as.vector(test[[py]]$Y_ts),
                                  # test.Y_ts  = NA,
                                  DATE_TIME  = rep(dimnames(test[[py]]$Y_ts)[[1]], times = test[[py]]$n),
                                  time.index = rep(1:length(dimnames(test[[py]]$Y_ts)[[1]]), times = test[[py]]$n),
                                  LON        = as.vector(t(matrix(test[[py]]$spCoordinates$LON, nrow = test[[py]]$n, ncol = test[[py]]$Nt))),
                                  LAT        = as.vector(t(matrix(test[[py]]$spCoordinates$LAT, nrow = test[[py]]$n, ncol = test[[py]]$Nt))))
          if(!is.null(test[[py]]$X_ts)){
            for(px in 1:dim(test[[py]]$X_ts)[1]){
              pred.data <- cbind(pred.data, as.vector(test[[py]]$X_ts[px,,]))
            }
          }
          colnames(pred.data) <- c( paste0(names(test)[py], ".Yts"), "DATE_TIME", "time.index", "LON", "LAT", dimnames(test[[py]]$X_ts)[[1]])

          pred.data.frame[[names(test)[py]]] <- pred.data

        }
      }
      spTaper <- spatial.taper(cs = cs, data)
    }
    # Initialize input ----
    ini.IEnKF.Input <- IniEnKF(data = data, Para.List = Para.List, test = test, IDE = IDE)
    cat("..........................\n")
    slope.stRF.ens <- slope.Hv.Zg <- list()

    PG <- ifelse(is.null(data[[1]]$Pg), 0, data[[1]]$Pg)
    for(py in 1:Py){
      PG <- PG + ifelse(is.null(data[[py]]$sPg), 0, data[[py]]$sPg)
    }


    G.Mat <- make.G.Mat(data = data, test = test)#, ini.IEnKF.Input = ini.IEnKF.Input)
    mcmc.random.effect.burnin <- 1e4
    mcmc.random.effect.draws <- setNames(vector("list", PG), names(G.Mat$data.G.Mat))
    mcmc.random.effect.iterations <- integer(0)
    GRFs.process <- NULL
    rmrf.Updating.order <- NULL

    # Sorting MRFs  ----
    if((PG > 1)){
      for(Pg in 1:PG){
        cat("\n ************************************************************\n")
        cat("* \n")
        cat(paste0("*  To test random fields in the [", names(G.Mat$data.G.Mat)[Pg], "]\n"))
        cat("* \n")
        index <- G.Mat$Y.index[[names(G.Mat$data.G.Mat)[Pg]]]
        py.inx <- G.Mat$Res.indic[[names(G.Mat$data.G.Mat)[Pg]]]
        mu.sigma.sq <- inv.mu.sigma.sq <- vector()
        if(py.inx == 0){
          for (py in 1:Py){
            mu.sigma.sq[data[[py]]$index.y]     <- Para.List[[names(data)[py]]]$obs.sigma.sq$mu.sigma.sq
            inv.mu.sigma.sq[data[[py]]$index.y] <- Para.List[[names(data)[py]]]$obs.sigma.sq$a/Para.List[[names(data)[py]]]$obs.sigma.sq$b
            # n          <- n + data[[py]]$n
            # n.diff[py] <- data[[py]]$n
          }
          spTaper0 <- spTaper[[1]]
        }else{
          spTaper0 <- spTaper[[py.inx]]
          mu.sigma.sq[1:length(index)]     <- Para.List[[names(data)[py.inx]]]$obs.sigma.sq$mu.sigma.sq
          inv.mu.sigma.sq[1:length(index)] <- Para.List[[names(data)[py.inx]]]$obs.sigma.sq$a/Para.List[[names(data)[py.inx]]]$obs.sigma.sq$b
        }


        slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- EnKF(data          = data,
                                                              yts.minus.xts = ini.IEnKF.Input$yts.minus.xts[, index],
                                                              Para.List     = Para.List,
                                                              H             = G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][, index, ],
                                                              # Mb = ini.IEnKF.Input$Mb,
                                                              slope.G.mat   = TRUE,
                                                              Mv            = ini.IEnKF.Input$Mv[[names(G.Mat$data.G.Mat)[Pg]]],
                                                              spTaper       = spTaper0,
                                                              ct            = ct,
                                                              Ne            = Ne,
                                                              var.name        = names(G.Mat$data.G.Mat)[Pg],
                                                              mu.sigma.sq     = mu.sigma.sq,
                                                              inv.mu.sigma.sq = inv.mu.sigma.sq,
                                                              Py              = py.inx)

        if(is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t, ,] %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                             ini.IEnKF.Input$train.Z.tsv[[t]]) %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }else if(is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                             ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                             ini.IEnKF.Input$train.Z.tsv[[t]],
                             ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }
        slope.Hv.Zg.sum <- slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]

        rmrf.Updating.order <- rbind(rmrf.Updating.order,
                                     data.table(Pg = Pg,
                                                Variable = names(G.Mat$data.G.Mat)[Pg],
                                                Loglik   =  loglik(data  = data,
                                                                   Para.List     = Para.List,
                                                                   yts.minus.xts = ini.IEnKF.Input$yts.minus.xts,
                                                                   slope.Hv.Zg   = slope.Hv.Zg.sum),
                                                w.max    = max(abs(slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[-1, ]))
                                     ))

      }
      rmrf.Updating.order.0 <- rmrf.Updating.order[, -c(1, 4)]
      rmrf.Updating.order$loglik_w.max <- rmrf.Updating.order$Loglik#/sum(rmrf.Updating.order$Loglik) #+ 0.2*rmrf.Updating.order$w.max/sum(rmrf.Updating.order$w.max)
      setorder(rmrf.Updating.order, -loglik_w.max) #rmrf.Updating.order$w.max/sum(rmrf.Updating.order$w.max)#rmrf.Updating.order$w.max/sum(rmrf.Updating.order$w.max)#
      cat("\n ************************************************************\n\n")
      cat("The order for updating random fields ... \n")
      print(rmrf.Updating.order)
      cat("\n\n")

      Order.list <- rmrf.Updating.order$Pg
    }else{
      Order.list <- 1
    }

    rmrf.Updating.order <- NULL
    # Recovering MRFs  ----
    for(Pg in Order.list){
      cat("\n ************************************************************\n")
      cat("* \n")
      cat(paste0("*  To recover random fields in the [", names(G.Mat$data.G.Mat)[Pg], "]\n"))
      cat("* \n")
      index  <- G.Mat$Y.index[[names(G.Mat$data.G.Mat)[Pg]]]
      py.inx <- G.Mat$Res.indic[[names(G.Mat$data.G.Mat)[Pg]]]
      mu.sigma.sq <- inv.mu.sigma.sq <- vector()
      if(py.inx == 0){
        for (py in 1:Py){
          mu.sigma.sq[data[[py]]$index.y]     <- Para.List[[names(data)[py]]]$obs.sigma.sq$mu.sigma.sq
          inv.mu.sigma.sq[data[[py]]$index.y] <- Para.List[[names(data)[py]]]$obs.sigma.sq$a/Para.List[[names(data)[py]]]$obs.sigma.sq$b
        }
        spTaper0 <- spTaper[[1]]
      }else{
        spTaper0                   <- spTaper[[py.inx]]
        mu.sigma.sq[1:length(index)]     <- Para.List[[names(data)[py.inx]]]$obs.sigma.sq$mu.sigma.sq
        inv.mu.sigma.sq[1:length(index)] <- Para.List[[names(data)[py.inx]]]$obs.sigma.sq$a/Para.List[[names(data)[py.inx]]]$obs.sigma.sq$b
      }
      if(Pg == Order.list[1]){
        slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- EnKF(data          = data,
                                                              yts.minus.xts = ini.IEnKF.Input$yts.minus.xts[, index],
                                                              H             = G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][, index, ],
                                                              Mv            = ini.IEnKF.Input$Mv[[names(G.Mat$data.G.Mat)[Pg]]],
                                                              Para.List     = Para.List,
                                                              slope.G.mat   = TRUE,
                                                              spTaper       = spTaper0,
                                                              ct            = ct,
                                                              Ne            = Ne,
                                                              var.name      = names(G.Mat$data.G.Mat)[Pg],
                                                              mu.sigma.sq     = mu.sigma.sq,
                                                              inv.mu.sigma.sq = inv.mu.sigma.sq,
                                                              Py              = py.inx)


        if(is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,] %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                             ini.IEnKF.Input$train.Z.tsv[[t]]) %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }else if(is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                             ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                             ini.IEnKF.Input$train.Z.tsv[[t]],
                             ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
              slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]
          }, simplify = "matrix") %>% t()
        }
        slope.Hv.Zg.sum <- slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
      }else{
        slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- EnKF(data         = data,
                                                              yts.minus.xts = ini.IEnKF.Input$yts.minus.xts[, index] - slope.Hv.Zg.sum[, index],
                                                              H             = G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][, index, ],
                                                              Mv            = ini.IEnKF.Input$Mv[[names(G.Mat$data.G.Mat)[Pg]]],
                                                              Para.List     = Para.List,
                                                              slope.G.mat   = TRUE,
                                                              spTaper       = spTaper0,
                                                              ct            = ct,
                                                              Ne            = Ne,
                                                              var.name      = names(G.Mat$data.G.Mat)[Pg],
                                                              mu.sigma.sq     = mu.sigma.sq,
                                                              inv.mu.sigma.sq = inv.mu.sigma.sq,
                                                              Py              = py.inx)

        slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
          spam::as.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,]) %*%
            slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ]}, simplify = "matrix") %>% t()
        slope.Hv.Zg.sum <- slope.Hv.Zg.sum + slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
      }

      rmrf.Updating.order <- rbind(rmrf.Updating.order,
                                   data.table(Pg = Pg,
                                              Variable = names(G.Mat$data.G.Mat)[Pg],
                                              Vt.max    = sum(abs(slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[-1, ]))))

    }

    Remove.Other.MRFs <- list()
    for(Pg in Order.list){
      Remove.Other.MRFs[[names(G.Mat$data.G.Mat)[Pg]]] <- ini.IEnKF.Input$yts.minus.xts -
        slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
    }
    if((PG > 1)){
      rmrf.Updating.order <- rmrf.Updating.order %>% left_join(rmrf.Updating.order.0, by = c("Variable"))
      rmrf.Updating.order$Vt.max <-
        rmrf.Updating.order$Vt.max/sum(rmrf.Updating.order$Vt.max)
      setorder(rmrf.Updating.order, -Vt.max)
      cat("\n ************************************************************\n\n")
      cat("The order for updating random fields ... \n")
      print(rmrf.Updating.order)
      cat("\n\n")

      Order.list <- rmrf.Updating.order$Pg

    }else{
      Order.list <- 1
    }
    # }

    slope.Au.Rg.sum <- 0

    # Likelihood ----
    ELBO.old <- loglik(data          = data,
                       Para.List     = Para.List,
                       yts.minus.xts = ini.IEnKF.Input$yts.minus.xts,
                       slope.Hv.Zg   = slope.Hv.Zg.sum)


    criterion <- TRUE; iter <- 0; ini.Para.List <- Para.List;
    public.beta.name <- NULL
    if(!is.null(data[[1]]$X_ts)){
      public.beta.name <- paste0("pub.", dimnames(data[[1]]$X_ts)[[1]])
    }



    beta.name         <- public.beta.name
    obs.sigma.sq.name <- paste0(names(data), ".obs.sigma.sq")
    RowName           <- vector()

    # PX <- ifelse(is.null(public.beta.name), length(self.beta.name) , length(public.beta.name))
    vb.EnKS.Iter <- matrix(NA, nrow = 1, ncol = (Py + length(beta.name)  + 5*data[[1]]$Grid.infor$summary$res))


    colnames(vb.EnKS.Iter) <- c("Iter",
                                beta.name,
                                obs.sigma.sq.name,
                                paste0("rho.v_",      1:data[[1]]$Grid.infor$summary$res),
                                paste0("Phi.v_",      1:data[[1]]$Grid.infor$summary$res),
                                paste0("tau.sq_",     1:data[[1]]$Grid.infor$summary$res),
                                # paste0("ini.tau.sq_", 1:data[[1]]$Grid.infor$summary$res),
                                "ELBO")
    public.beta <- NULL
    if(!is.null(data[[1]]$X_ts)){
      public.beta <- as.vector(Para.List[[names(data)[[1]]]]$beta$mu.beta)
    }



    mu.sigma.sq <- NULL
    for(py in 1:Py){
      mu.sigma.sq <- c(mu.sigma.sq, Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq)
    }
    beta <- c(public.beta)
    if(is.null(Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]])){
      show.list <- Para.List[[names(data)[[1]]]]$st.sRF[[names(G.Mat$data.G.Mat)[1]]]
    }else{
      show.list <- Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]]
    }
    vb.EnKS.Iter[1, ] <- c(0,
                           round(beta, digist.num),
                           round(mu.sigma.sq, digist.num),
                           round(show.list$rho.v$mu.rho.v, digist.num),
                           round(show.list$Phi.v$mu.Phi.v[1], digist.num),
                           round(show.list$proc.tau.sq$mu.tau.sq, digist.num),
                           # round(show.list$proc.ini.tau.sq$mu.tau.sq, digist.num),
                           round(ELBO.old, 3))

  }

  #--plot
  {
    public.beta <- NULL
    if(!is.null(data[[1]]$X_ts)){
      public.beta <- as.vector(Para.List[[names(data)[[1]]]]$beta$mu.beta)
    }

    mu.sigma.sq <- NULL
    for(py in 1:Py){
      mu.sigma.sq <- c(mu.sigma.sq, Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq)
    }
    beta   <- c(public.beta)
    n.grid <- length(beta)  + 4
    if(plot){
      nr <- floor(sqrt(n.grid))
      for(i in 0:2){
        nc <- nr + i
        if(nr*nc > (n.grid - 1)) break;
      }
      #mar: bottom-left-top-right #mgp: space between xlab and plot- axis ticks - axis
      par(mar = c(2, 2, 2, 1) + 0.5, mfrow = c(2, 3), mgp = c(1.5, 0.5, 0)) #c(nr, nc)
      # par(cex = 1.0)
      par(cex      = ifelse(n.grid > 20, 0.7, 1.0),
          cex.axis = ifelse(n.grid > 20, 0.7, 1.0),
          cex.lab  = ifelse(n.grid > 20, 1.0, 1.3),
          cex.main = ifelse(n.grid > 20, 1.0, 1.3),
          lwd      = 2)
    }

  }

  fitted.slope.Hv.Zg.sum <- Relative.ELBO.error <- NULL
  while (criterion & iter < itMax) {
    cat("..........................\n")
    start_time <- proc.time()

    pb         <- progress::progress_bar$new(format = "S.Mat: (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                             , total  = PG
                                             , clear  = FALSE)
    S <- Sw <-list()
    for(Pg in 1:PG){
      py.inx <- G.Mat$Res.indic[[names(G.Mat$data.G.Mat)[Pg]]]
      if(py.inx == 0){
        spTaper0 <- spTaper[[1]]
      }else{
        spTaper0 <- spTaper[[py.inx]]
      }
      pb$tick()
      Sys.sleep(1 / (2*PG))
      S[[names(G.Mat$data.G.Mat)[Pg]]] <- sfunc(data,
                                                Ks      = slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]],
                                                spTaper = spTaper0,
                                                cores   = n.cores,
                                                py = py.inx) #n.cores





    }

    # VB.Laplace ----
    VB.Laplace <- .VB(data            = data,
                      Hv.Zg           = slope.Hv.Zg.sum,
                      data.trans.tsv  = data.trans.tsv,
                      Ks              = slope.stRF.ens,
                      S               = S,
                      Para.List       = Para.List,
                      prior           = prior,
                      verbose.VB      = verbose.VB,
                      test            = test,
                      Ne              = Ne,
                      iter            = iter + 1)
    Para.List <- VB.Laplace$Para.List
    data <- VB.Laplace$data
    test <- VB.Laplace$test
    G.Mat <- make.G.Mat(data = data, test = test)

    t1 <- proc.time()

    ini.IEnKF.Input <- IniEnKF(data = data, Para.List = VB.Laplace$Para.List, test = test, IDE = IDE)

    rmrf.Updating.order <- NULL

    # Update MRFs ----
    for(Pg in Order.list){
      cat("************************************************************\n")
      cat("* \n")
      cat(paste0("*  To recover spatiotemporal GMRFs in the [", names(G.Mat$data.G.Mat)[Pg], "]\n"))
      cat("* \n")
      # cat("*-----------------------------------------------------------\n")
      slope.Hv.Zg.sum <- slope.Hv.Zg.sum - slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
      index           <- G.Mat$Y.index[[names(G.Mat$data.G.Mat)[Pg]]]
      py.inx          <- G.Mat$Res.indic[[names(G.Mat$data.G.Mat)[Pg]]]
      mu.sigma.sq     <- inv.mu.sigma.sq <- vector()
      if(py.inx == 0){
        for (py in 1:Py){
          mu.sigma.sq[data[[py]]$index.y]     <- Para.List[[names(data)[py]]]$obs.sigma.sq$mu.sigma.sq
          inv.mu.sigma.sq[data[[py]]$index.y] <- Para.List[[names(data)[py]]]$obs.sigma.sq$a/Para.List[[names(data)[py]]]$obs.sigma.sq$b
        }
        spTaper0 <- spTaper[[1]]
      }else{
        spTaper0 <- spTaper[[py.inx]]
        mu.sigma.sq[1:length(index)]     <- Para.List[[names(data)[py.inx]]]$obs.sigma.sq$mu.sigma.sq
        inv.mu.sigma.sq[1:length(index)] <- Para.List[[names(data)[py.inx]]]$obs.sigma.sq$a/Para.List[[names(data)[py.inx]]]$obs.sigma.sq$b
      }
      if(Pg == Order.list[1]){
        var.name  <- names(G.Mat$data.G.Mat)[Pg]
        rf.gamma  <- VB.Laplace$Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[Pg]]]$gamma
        post.prob <- 1
        for(py in 1:Py){
          if(var.name %in% names(VB.Laplace$Para.List[[names(data)[[py]]]]$st.sRF)){
            rf.gamma <- VB.Laplace$Para.List[[names(data)[[py]]]]$st.sRF[[names(G.Mat$data.G.Mat)[Pg]]]$gamma
            post.prob <- 1
          }
        }

        slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- EnKF(data          = data,
                                                              yts.minus.xts = ini.IEnKF.Input$yts.minus.xts[, index] -
                                                                slope.Hv.Zg.sum[, index],
                                                              H             = G.Mat$data.G.Mat[[var.name]][, index, ],
                                                              Mv            = ini.IEnKF.Input$Mv[[names(G.Mat$data.G.Mat)[Pg]]],
                                                              Para.List     = VB.Laplace$Para.List,
                                                              slope.G.mat   = TRUE,
                                                              spTaper       = spTaper0,
                                                              ct            = ct,
                                                              Ne            = Ne,
                                                              var.name      = var.name,
                                                              mu.sigma.sq     = mu.sigma.sq,
                                                              inv.mu.sigma.sq = inv.mu.sigma.sq,
                                                              Py              = py.inx)


        if(is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
            G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,] %*% (slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ])*post.prob
          }, simplify = "matrix") %>% t()
        }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <-   sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,], ini.IEnKF.Input$train.Z.tsv[[t]]) %*%
              (slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ])*post.prob

          }, simplify = "matrix") %>% t()
        }else if(is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <-  sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,], ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
              (slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ])*post.prob
          }, simplify = "matrix") %>% t()
        }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
          slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <-  sapply(seq_len(Nt), function(t) {
            spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,], ini.IEnKF.Input$train.Z.tsv[[t]],
                             ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
              (slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ])*post.prob
          }, simplify = "matrix") %>% t()
        }
        slope.Hv.Zg.sum <- slope.Hv.Zg.sum + slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
      }else{
        var.name <-  names(G.Mat$data.G.Mat)[Pg]
        rf.gamma <- VB.Laplace$Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[Pg]]]$gamma
        post.prob <- 1
        for(py in 1:Py){
          if(var.name %in% names(VB.Laplace$Para.List[[names(data)[[py]]]]$st.sRF)){
            rf.gamma <- VB.Laplace$Para.List[[names(data)[[py]]]]$st.sRF[[names(G.Mat$data.G.Mat)[Pg]]]$gamma
            post.prob <- 1
          }
        }
        slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- EnKF(data            = data,
                                                              yts.minus.xts = ini.IEnKF.Input$yts.minus.xts[, index] -
                                                                slope.Hv.Zg.sum[, index], #- slope.Au.Rg.sum,
                                                              H             = G.Mat$data.G.Mat[[var.name]][, index, ],
                                                              Mv            = ini.IEnKF.Input$Mv[[names(G.Mat$data.G.Mat)[Pg]]],
                                                              Para.List     = VB.Laplace$Para.List,
                                                              slope.G.mat   = TRUE,
                                                              spTaper       = spTaper0,
                                                              ct            = ct,
                                                              Ne            = Ne,
                                                              var.name      = var.name,
                                                              mu.sigma.sq     = mu.sigma.sq,
                                                              inv.mu.sigma.sq = inv.mu.sigma.sq,
                                                              Py              = py.inx)
        slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]] <- sapply(seq_len(Nt), function(t) {
          spam::as.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,]) %*%
            (slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, ])*post.prob
        }, simplify = "matrix") %>% t()
        slope.Hv.Zg.sum <- slope.Hv.Zg.sum + slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
      }


      rmrf.Updating.order <- rbind(rmrf.Updating.order,
                                   data.table(Pg = Pg,
                                              Variable  = names(G.Mat$data.G.Mat)[Pg],
                                              Vt.max    = max(abs(slope.stRF.ens[[var.name]]$Vt.mean[-1, ]))))
    }

    Remove.Other.MRFs <- list()
    for(Pg in Order.list){
      Remove.Other.MRFs[[names(G.Mat$data.G.Mat)[Pg]]] <- ini.IEnKF.Input$yts.minus.xts -
        slope.Hv.Zg[[names(G.Mat$data.G.Mat)[Pg]]]
    }

    if((PG > 1)){
      rmrf.Updating.order        <- rmrf.Updating.order %>% left_join(rmrf.Updating.order.0, by = c("Variable"))
      rmrf.Updating.order$Vt.max <- #0.8*rmrf.Updating.order$Loglik/sum(rmrf.Updating.order$Loglik) +
        rmrf.Updating.order$Vt.max/sum(rmrf.Updating.order$Vt.max)
      setorder(rmrf.Updating.order, -Vt.max)
      cat("\n ************************************************************\n\n")
      cat("The order for updating random fields ... \n")
      print(rmrf.Updating.order)
      cat("\n\n")
      Order.list <- rmrf.Updating.order$Pg
    }else{
      Order.list <-  1
    }

    v.pre <- Matrix::Diagonal(0)
    for(py in 1:Py){
      mu.alpha <- Para.List[[names(data)[py]]]$st.sRF[[paste0(names(data)[[py]], ".Intercept")]]$LR.coef$mu.alpha
      mat.x    <- as.matrix(data[[py]]$Grid.infor$var.covariable)
      mu.var   <- exp(mat.x %*% mu.alpha)
      rf.corr  <- exp(-data[[py]]$Grid.infor$level[[1]]$BAUs.Dist/Para.List[[names(data)[py]]]$st.sRF[[paste0(names(data)[[py]], ".Intercept")]]$Phi.v$mu.Phi.v[1])
      rf.corr  <-  diag(as.vector(mu.var^(0.5))) %*% rf.corr %*% diag(as.vector(mu.var^(0.5)))
      Q        <- Para.List[[py]]$st.sRF[[paste0(names(data)[[py]], ".Intercept")]]$proc.tau.sq$mu.tau.sq[1]*solve(rf.corr) #
      v.pre <- Matrix::bdiag(v.pre,  Q)
    }
    S11 <- S10 <- S00 <-Matrix::Diagonal(0)
    q.sigma.sq <- 0
    for(py in 1:Py){
      q.sigma.sq <- q.sigma.sq -lgamma(VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$a) - log(VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$b) + (VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$a + 1) * digamma(VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$a) -  VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$a

      rf.name     <- dimnames(data[[names(data)[[py]]]]$sG_ts)[[1]]
      sub.rf.name <- paste0(names(data)[[py]], ".", rf.name[1])
      mu.rho.v  <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[1]
      var.rho.v <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[1]

      S11 <- Matrix::bdiag(S11,  as.matrix(S[[sub.rf.name]]$S11[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[1]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[1]]]))
      S10 <- Matrix::bdiag(S10,  mu.rho.v*as.matrix(S[[sub.rf.name]]$S10[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[1]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[1]]]))
      S00 <- Matrix::bdiag(S00,  (mu.rho.v^2 + var.rho.v)*as.matrix(S[[sub.rf.name]]$S00[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[1]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[1]]]))
    }

    q.beta <- -0.5* as.vector(determinant(VB.Laplace$Para.List[[names(data)[[1]]]]$beta$sigma.sq)$modulus) -
      0.5*sum(diag(VB.Laplace$Para.List[[names(data)[[1]]]]$beta$sigma.sq %*% (VB.Laplace$Para.List[[names(data)[[1]]]]$beta$mu.beta[, 1] %*%
                                                                                 t(VB.Laplace$Para.List[[names(data)[[1]]]]$beta$mu.beta[, 1])))) +
      as.vector(VB.Laplace$Para.List[[names(data)[[1]]]]$beta$mu.beta[, 1] %*% VB.Laplace$Para.List[[names(data)[[1]]]]$beta$sigma.sq %*%
                  VB.Laplace$Para.List[[names(data)[[1]]]]$beta$mu.beta[, 1])

    q.v        <- -0.5*(sum(diag(v.pre %*% S11)) - 2*sum(diag(v.pre %*% S10)) + sum(diag(v.pre %*% S00)))  +  0.5*(Nt + 1) *Matrix::determinant(v.pre)$modulus

    # Update likelihhod ----
    ELBO.new <- loglik(data          = data,
                       Para.List     = VB.Laplace$Para.List,
                       yts.minus.xts = ini.IEnKF.Input$yts.minus.xts,
                       slope.Hv.Zg   = slope.Hv.Zg.sum) - q.beta - q.v - q.sigma.sq


    mid_time <- proc.time()
    if (verbose) {
      cat("The running time for the spEnKS algorithm: \n")
      print(mid_time - start_time)
      cat("................................................................. \n")
    }


    # Update ELBO.error ----
    ELBO.error <- abs((ELBO.new -  ELBO.old)/ ELBO.old)
    end_time <- proc.time()

    if(iter >= itMin){
      criterion <- (ELBO.error > tol.real)
    }

    public.beta <- NULL
    if(!is.null(data[[1]]$X_ts)){
      public.beta <- as.vector(VB.Laplace$Para.List[[names(data)[[1]]]]$beta$mu.beta)
    }

    mu.sigma.sq <- NULL
    for(py in 1:Py){
      mu.sigma.sq <- c(mu.sigma.sq, VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq)
    }
    beta <- c(public.beta)


    if((is.logical(VB.Laplace$Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]]$rho.v$var.rho.v))|(
      is.null(VB.Laplace$Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]]$rho.v$var.rho.v)
    )){
      show.list <-  VB.Laplace$Para.List[[names(data)[[1]]]]$st.sRF[[names(G.Mat$data.G.Mat)[1]]]

    }else{
      show.list <- VB.Laplace$Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]]
    }

    # for(py in 1:Py){
    vb.EnKS.Iter <- rbind(as.data.frame(vb.EnKS.Iter), c(iter + 1,
                                                         round(beta, digist.num),
                                                         mu.sigma.sq,
                                                         round(show.list$rho.v$mu.rho.v, digist.num),
                                                         round(show.list$Phi.v$mu.Phi.v[1], digist.num),
                                                         round(show.list$proc.tau.sq$mu.tau.sq, digist.num),
                                                         # round(show.list$proc.ini.tau.sq$mu.tau.sq, digist.num),
                                                         round(ELBO.new, 3)))

    # }

    Name <- colnames(vb.EnKS.Iter)
    inx  <- as.vector(which(sapply(vb.EnKS.Iter, anyNA)))

    out <- vb.EnKS.Iter
    if(length(inx) >= 1){
      out <- as.data.table(out[, -inx])
    }

    print(tail(out, 10))

    setDF(out)

    setDF(vb.EnKS.Iter)
    if(plot){
      for(i in c(2:3, 4:5, 10, ncol(out))){ #c(2:5, ncol(out))c(2:(length(beta) + 4), ncol(out) - 2, ncol(out))
        x <- 1:(iter + 2)
        y <- out[, i]
        plot(x, y[x], ylab = colnames(out)[i], type = "l", xlab = "Iteration")
        if(i == 2){
          abline(h = 5)
        }else{
          abline(h = -1)
        }
      }
    }


    # Fit model ----
    RMSE.1 <- CRPS.1 <- Coef.1 <- FAC2.1 <- MAE.1 <-  Coverage.1 <-
      RMSE.2 <- CRPS.2 <- Coef.2 <- FAC2.2 <- MAE.2 <-  Coverage.2 <- NA
    if((iter == 0 ) | (save.Predict) | (iter == (itMax - 1))|(isFALSE(criterion))){
      cat("************************************************************\n")
      cat("* \n")
      cat("* \n")
      cat("*                      Fit model                            \n")
      cat("* \n")
      cat("* \n")
      cat("************************************************************\n")
      fitted.Pred.ens <- list()
      fitted.Pred.trans.ens <- list()
      Goodness.of.fit.by.response <- vector("list", Py)
      names(Goodness.of.fit.by.response) <- names(data)
      {

        fitted.slope.Hv.Zg.sum <- 0
        fitted.slope.Hv.Zg.ens <- list()
        for(Pg in 1:(PG)){
          fitted.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- array(NA,
                                                                         dim = c(Nt, n, Ne),
                                                                         dimnames = list(time, ID, as.character(paste0("Ens.", 1:Ne))))
          if(Pg == 1){
            if(is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
              temp <- sapply(seq_len(Nt), function(t) {
                G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,] %*%
                  slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
              }, simplify = "array")
            }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&is.null(ini.IEnKF.Input$train.sZ.tsv)){
              temp <- sapply(seq_len(Nt), function(t) {
                spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                                 ini.IEnKF.Input$train.Z.tsv[[t]]) %*%
                  slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
              }, simplify = "array")
            }else if(is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
              temp <- sapply(seq_len(Nt), function(t) {
                spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                                 ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
                  slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
              }, simplify = "array")
            }else if(!is.null(ini.IEnKF.Input$train.Z.tsv)&!is.null(ini.IEnKF.Input$train.sZ.tsv)){
              temp <- sapply(seq_len(Nt), function(t) {
                spam::cbind.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                                 ini.IEnKF.Input$train.Z.tsv[[t]],
                                 ini.IEnKF.Input$train.sZ.tsv[[t]]) %*%
                  slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
              }, simplify = "array")
            }

            for(t in 1:Nt){
              fitted.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]][t,,] <- temp[,,t]
            }
            fitted.slope.Hv.Zg.sum <- fitted.slope.Hv.Zg.sum + fitted.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]]
          }else{
            temp <- sapply(seq_len(Nt), function(t) {
              spam::as.spam(G.Mat$data.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,]) %*%
                slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
            }, simplify = "array")
            for(t in 1:Nt){
              fitted.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]][t,,] <- temp[,,t]
            }
            fitted.slope.Hv.Zg.sum <- fitted.slope.Hv.Zg.sum + fitted.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]]
          }
        }


        fitted.Xts.ens <- array(0, dim = c(Nt, n, Ne), dimnames = list(time, ID, paste0("En_", 1:Ne)))
        if(!is.null(data[[1]]$X_ts)){
          sigma <- VB.Laplace$Para.List[[names(data)[[py]]]]$beta$sigma.sq
          chol.sigma <- tryCatch({chol(sigma)},
                                 error = function(e)
                                   print("ERROR"))
          if(length(chol.sigma) == 1){
            diag(sigma) <- diag(sigma) + 1e5
          }
          beta.sample <- mvnfast::rmvn(Ne, mu = VB.Laplace$Para.List[[names(data)[[py]]]]$beta$mu.beta,
                                       sigma =  sigma)

          Xts.beta.ens <- mcmapply(t = seq_len(Nt), SIMPLIFY = F,
                                   FUN = function(t){
                                     X_ts <- data.trans.tsv$tsv.x[[t]];
                                     Matrix::crossprod(x = t(beta.sample),
                                                       y = X_ts);}
                                   , mc.cores = 1)
          for(t in 1:Nt){
            fitted.Xts.ens[t,, ] <- as.matrix(t(Xts.beta.ens[[t]]))

          }
        }

        for(py in 1:Py){
          e <- rnorm(Nt*length(data[[py]]$index.y)*Ne, sd = sqrt(VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq))
          fitted.Pred.ens[[py]] <- fitted.Xts.ens[, data[[py]]$index.y,] +
            # ifelse(data$Pr > 0, slope.sRF.Au.Rg[[dimnames(data$R_ts)[[1]][Pr]]]$fitted.AU.Ens, 0) +
            fitted.slope.Hv.Zg.sum[, data[[py]]$index.y, ] #+ e

          fitted.Pred.trans.ens[[py]] <- fitted.Pred.ens[[py]]



          cat(paste0("Fitted values of the response ***", names(data)[py], "*** [", Round(min(fitted.Pred.ens[[py]]), 5), " ",
                     Round(max(fitted.Pred.ens[[py]]), 5), "]\n"))
          cat(paste0("---\n"))

          if(center.y&scale.y){
            fitted.Pred.ens[[py]] <-  fitted.Pred.ens[[py]]*sd.y[py] + mean.y[py]
          }else if(!center.y&scale.y){
            fitted.Pred.ens[[py]] <-  fitted.Pred.ens[[py]]*sd.y[py]
          }else if(center.y&!scale.y){
            fitted.Pred.ens[[py]] <-  fitted.Pred.ens[[py]] + mean.y[py]
          }

          if (transf.Response %in% c("SQRT", "sr", "sqrt", "Sqrt")) {
            ind.len <- length(which(fitted.Pred.ens[[py]] < 0))
            fitted.Pred.ens[[py]] <- ifelse(fitted.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), fitted.Pred.ens[[py]])
            fitted.Pred.ens[[py]] <- fitted.Pred.ens[[py]]^2
          }else if (transf.Response %in% c("log", "Log", "LOG")) {
            fitted.Pred.ens[[py]] <- exp(fitted.Pred.ens[[py]]) #- 1
            ind.len <- length(which(fitted.Pred.ens[[py]] < 0))
            fitted.Pred.ens[[py]] <- ifelse(fitted.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), fitted.Pred.ens[[py]])
          }else if (transf.Response %in% c("cubic")) {
            fitted.Pred.ens[[py]] <- (fitted.Pred.ens[[py]])^3
            ind.len <- length(which(fitted.Pred.ens[[py]] < 0))
            fitted.Pred.ens[[py]] <- ifelse(fitted.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), fitted.Pred.ens[[py]])
          }else if (transf.Response %in% c("NORMAL", "Normal", "normal")) {
            ind.len <- length(which(fitted.Pred.ens[[py]] < 0))
            if(positive){
              fitted.Pred.ens[[py]] <- ifelse(fitted.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), fitted.Pred.ens[[py]])
            }
          }else if(transf.Response %in% c("loglog")){

            fitted.Pred.ens[[py]] <- 1 - exp(-exp(fitted.Pred.ens[[py]]))



            ind.len <- length(which(fitted.Pred.ens[[py]] < 0))
            fitted.Pred.ens[[py]] <- ifelse(fitted.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), fitted.Pred.ens[[py]])
          }

          up <- log(max(true.Y_ts[[py]]))^(1/1) + 2

          Goodness.of.fit.by.response[[py]] <- compute_good_of_fitting(
            obs = data[[py]]$Y_ts,
            eta.samples = fitted.Pred.trans.ens[[py]],
            sigma.info = VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq
          )

          cat(paste0("Transform fitted values ***", names(data)[py],
                     "*** to original scale [", Round(min(fitted.Pred.ens[[py]]), 5), " ",
                     Round(max(fitted.Pred.ens[[py]]), 5), "]\n"))
          cat(paste0("---\n"))
          fitted.Pred <- apply(fitted.Pred.ens[[py]], c(1, 2),
                               quantile, prob = c(0.025, 0.5, 0.975))

          fitted.Pred.sd <- apply(fitted.Pred.ens[[py]], 1:2, sd)

          fitted.spT <- spT_validation(z        = true.Y_ts[[py]] ,
                                       zhat     = fitted.Pred[2,,],
                                       sigma    = fitted.Pred.sd,#,
                                       zhat.Ens = fitted.Pred.ens[[py]],
                                       names    = F,
                                       CC       = F)
          # }
          true.da <- as.vector(true.Y_ts[[py]])
          inx <- which(true.da > 0)

          if(plot&(py == 1)){
            plot(true.da, as.vector(fitted.Pred[2,,]), pch = 19, cex = 0.5,
                 xlab = "True values", ylab = "Fitted values",
                 xlim = c(min(true.da, as.vector(fitted.Pred[2,,]), na.rm = T), max(true.da, as.vector(fitted.Pred[2,,]), na.rm = T)),
                 ylim = c(min(true.da, as.vector(fitted.Pred[2,,]), na.rm = T), max(true.da, as.vector(fitted.Pred[2,,]), na.rm = T)))

          }

          # for(py in 1:Py){
          if(nrow(train.Res.miss.index[[py]]) > 0) {
            cat("Missing for training datasets will be updated via JSTVC ...\n")
            for(i in 1:nrow(train.Res.miss.index[[py]])){
              t     <-  train.Res.miss.index[[py]][i, 1]
              s.ind <-  train.Res.miss.index[[py]][i, 2]
              if(t %in% c(2)){
                ini.miss.value <- quantile(as.vector(data[[py]]$Y_ts[t, ]), probs = Miss.ini[1], na.rm = TRUE)#max(data[[py]]$Y_ts, na.rm = TRUE)
              }
              if(t %in% c(3)){
                ini.miss.value <- quantile(as.vector(data[[py]]$Y_ts[t, ]), probs = Miss.ini[2], na.rm = TRUE)#max(data[[py]]$Y_ts, na.rm = TRUE)
              }
              if(t %in% c(4)){
                ini.miss.value <- quantile(as.vector(data[[py]]$Y_ts[t, ]), probs = Miss.ini[3], na.rm = TRUE)
              }

              fill.value <- as.vector(fitted.Pred[2, t, s.ind])
              ifelse(fill.value > max(data[[py]]$Y_ts, na.rm = TRUE),
                     ini.miss.value,
                     fill.value)
            }


            Fill.Res.Miss.Values[[names(data)[py]]] <-
              fitted.Pred[2,train.Res.miss.index[[py]][, 1], train.Res.miss.index[[py]][, 2]]
          }else {
            Fill.Res.Miss.Values <- NULL
          }
          # }

          # print(fitted.spT)
          fitted.spT <- as.vector(fitted.spT)
          RMSE.1     <- fitted.spT[1]
          Coef.1     <- fitted.spT[2]
          FAC2.1     <- fitted.spT[3]
          CRPS.1     <- fitted.spT[4]
          MAE.1      <- fitted.spT[5]
          if(verbose){
            cat(paste0("Mean absolute error (MAE) [", Round(MAE.1, 4),"] \n"
                       , "Root mean squared error (RMSE) [", RMSE.1, "] \n"
                       , "Continuous rank probability score (CRPS) [", CRPS.1, "] \n"
                       , "Fraction of predictions within a factor of two (FAC2) [", FAC2.1,"] \n"
                       , "Pearson correlation between the prediction and obs [", Coef.1,"] \n\n"
            ))
          }
        }
      }
    }

    Goodness.of.fit.loglik <- do.call(
      rbind,
      lapply(Goodness.of.fit.by.response, function(x) x$loglik.matrix)
    )
    if (!is.null(Goodness.of.fit.loglik)) {
      overall.mean.deviance <- mean(-2 * colSums(Goodness.of.fit.loglik))
      overall.deviance.of.mean <- sum(vapply(Goodness.of.fit.by.response, function(x) x$deviance.of.mean, numeric(1)))
      overall.p.eff.DIC <- overall.mean.deviance - overall.deviance.of.mean
      overall.lppd <- sum(log(rowMeans(exp(Goodness.of.fit.loglik))))
      overall.p.eff.WAIC <- sum(apply(Goodness.of.fit.loglik, 1, var))
      Goodness.of.fit.overall <- list(
        DIC = overall.mean.deviance + overall.p.eff.DIC,
        p.eff.DIC = overall.p.eff.DIC,
        mean.deviance = overall.mean.deviance,
        deviance.of.mean = overall.deviance.of.mean,
        WAIC = -2 * (overall.lppd - overall.p.eff.WAIC),
        p.eff.WAIC = overall.p.eff.WAIC,
        lppd = overall.lppd
      )
    } else {
      Goodness.of.fit.overall <- list(
        DIC = NA_real_,
        p.eff.DIC = NA_real_,
        mean.deviance = NA_real_,
        deviance.of.mean = NA_real_,
        WAIC = NA_real_,
        p.eff.WAIC = NA_real_,
        lppd = NA_real_
      )
    }
    Goodness.of.fit <- list(
      Overall = Goodness.of.fit.overall,
      By.response = lapply(Goodness.of.fit.by.response, function(x) {
        x[names(x) != "loglik.matrix"]
      })
    )

    # Test model ----
    if((iter == 0) | (save.Predict) |(iter == (itMax - 1))|(isFALSE(criterion))){
      if (!is.null(test)) {
        cat("************************************************************\n")
        cat("* \n")
        cat("* \n")
        cat("*                       Test model                          \n")
        cat("* \n")
        cat("* \n")
        cat("************************************************************\n")

        test.Pred.ens <- list()
        {

          pred.slope.Hv.Zg.sum <- 0
          Pred.slope.Hv.Zg.ens <- list()
          for(Pg in 1:(PG)){
            Pred.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]] <- array(NA, dim  =  c(Nt, length(ID.test), Ne),
                                                                         dimnames = list(time, ID.test, as.character(paste0("Ens.", 1:Ne))))

            if(Pg == 1){
              if(is.null(ini.IEnKF.Input$test.Z.tsv)&is.null(ini.IEnKF.Input$test.sZ.tsv)){
                temp <- sapply(seq_len(Nt), function(t) {
                  G.Mat$test.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,] %*%
                    slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
                }, simplify = "array")
              }else if(!is.null(ini.IEnKF.Input$test.Z.tsv)&is.null(ini.IEnKF.Input$test.sZ.tsv)){
                temp <- sapply(seq_len(Nt), function(t) {
                  spam::cbind.spam(G.Mat$test.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                                   ini.IEnKF.Input$test.Z.tsv[[t]]) %*%
                    slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
                }, simplify = "array")
              }else if(is.null(ini.IEnKF.Input$test.Z.tsv)&!is.null(ini.IEnKF.Input$test.sZ.tsv)){
                temp <- sapply(seq_len(Nt), function(t) {
                  spam::cbind.spam(G.Mat$test.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                                   ini.IEnKF.Input$test.sZ.tsv[[t]]) %*%
                    slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
                }, simplify = "array")
              }else if(!is.null(ini.IEnKF.Input$test.Z.tsv)&!is.null(ini.IEnKF.Input$test.sZ.tsv)){
                temp <- sapply(seq_len(Nt), function(t) {
                  spam::cbind.spam(G.Mat$test.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,],
                                   ini.IEnKF.Input$test.Z.tsv[[t]],
                                   ini.IEnKF.Input$test.sZ.tsv[[t]]) %*%
                    slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
                }, simplify = "array")
              }


              for(t in 1:Nt){
                Pred.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]][t,,] <- temp[,,t]
              }
              pred.slope.Hv.Zg.sum <- pred.slope.Hv.Zg.sum + Pred.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]]
            }else{
              temp <- sapply(seq_len(Nt), function(t) {
                spam::as.spam(G.Mat$test.G.Mat[[names(G.Mat$data.G.Mat)[Pg]]][t,,]) %*%
                  slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.ens[t + 1, ,]
              }, simplify = "array")
              for(t in 1:Nt){
                Pred.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]][t,,] <- temp[,,t]
              }
              pred.slope.Hv.Zg.sum <- pred.slope.Hv.Zg.sum + Pred.slope.Hv.Zg.ens[[names(G.Mat$data.G.Mat)[Pg]]]
            }
          }


          test.Xts.ens <- array(0, dim   = c(Nt, length(ID.test), Ne),
                                dimnames = list(time, ID.test, paste0("En_", 1:Ne)))

          if(!is.null(test[[1]]$X_ts)){
            sigma <- VB.Laplace$Para.List[[names(data)[[py]]]]$beta$sigma.sq
            chol.sigma <- tryCatch({chol(sigma)},
                                   error = function(e)
                                     print("ERROR"))
            if(length(chol.sigma) == 1){
              diag(sigma) <- diag(sigma) + 1e5
            }
            beta.sample <- mvnfast::rmvn(Ne, mu = VB.Laplace$Para.List[[names(data)[[py]]]]$beta$mu.beta,
                                         sigma  = sigma)


            Xts.beta.ens <- mcmapply(t          = seq_len(Nt),
                                     SIMPLIFY   = F,
                                     FUN        = function(t){
                                       X_ts <- test.trans.tsv$tsv.x[[t]];
                                       Matrix::crossprod(x = t(beta.sample),
                                                         y = X_ts);}
                                     , mc.cores = 1)

            for(t in 1:Nt){
              test.Xts.ens[t,, ] <- as.matrix(t(Xts.beta.ens[[t]]))
            }
          }


          for(py in 1:Py){
            e <- rnorm(Nt*length(test[[py]]$index.y)*Ne, sd = sqrt(VB.Laplace$Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq))
            test.Pred.ens[[py]] <- test.Xts.ens[, test[[py]]$index.y,] +
              pred.slope.Hv.Zg.sum[, test[[py]]$index.y, ] + e

            cat(paste0("Predicted values of the response ***", names(data)[py], "*** [", Round(min(test.Pred.ens[[py]]), 5), " ",
                       Round(max(test.Pred.ens[[py]]), 5), "]\n"))
            cat(paste0("---\n"))


            if(center.y&scale.y){
              test.Pred.ens[[py]] <-  test.Pred.ens[[py]]*sd.y[py] + mean.y[py]
            }else if(!center.y&scale.y){
              test.Pred.ens[[py]]   <-  test.Pred.ens[[py]]*sd.y[py]
            }else if(center.y&!scale.y){
              test.Pred.ens[[py]]   <-  test.Pred.ens[[py]] + mean.y[py]
            }

            if (transf.Response %in% c("SQRT", "sr", "sqrt", "Sqrt")) {
              ind.len <- length(which(test.Pred.ens[[py]] < 0))
              test.Pred.ens[[py]] <- ifelse(test.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), test.Pred.ens[[py]])
              test.Pred.ens[[py]] <- test.Pred.ens[[py]]^2
            }else if (transf.Response %in% c("log", "Log", "LOG")) {
              test.Pred.ens[[py]] <- exp(test.Pred.ens[[py]]) #- 1
              ind.len <- length(which(test.Pred.ens[[py]] < 0))
              test.Pred.ens[[py]] <- ifelse(test.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), test.Pred.ens[[py]])
            }else if (transf.Response %in% c("cubic")) {
              test.Pred.ens[[py]] <- (test.Pred.ens[[py]])^3
              ind.len <- length(which(test.Pred.ens[[py]] < 0))
              test.Pred.ens[[py]] <- ifelse(test.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), test.Pred.ens[[py]])
            }else if (transf.Response %in% c("NORMAL", "Normal", "normal")) {
              ind.len <- length(which(test.Pred.ens[[py]] < 0))
              if(positive){
                test.Pred.ens[[py]] <- ifelse(test.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), test.Pred.ens[[py]])
              }

            }else if(transf.Response %in% c("loglog")){
              # test.Pred.ens[[py]] <- exp(5 - exp(-test.Pred.ens[[py]]))

              test.Pred.ens[[py]] <- 1 - exp(-exp(test.Pred.ens[[py]]))


              ind.len <- length(which(test.Pred.ens[[py]] < 0))
              test.Pred.ens[[py]] <- ifelse(test.Pred.ens[[py]] < 0, runif(ind.len, 1e-8, 1e-5), test.Pred.ens[[py]])
            }


            cat(paste0("Transform predictions ***",  names(data)[py],
                       "*** to original scale [", Round(min(test.Pred.ens[[py]]), 5), " ",
                       Round(max(test.Pred.ens[[py]]), 5), "]\n"))
            cat(paste0("---\n"))
            test.Pred    <- apply(test.Pred.ens[[py]], c(1, 2), quantile, prob = c(0.025, 0.5, 0.975))
            test.Pred.sd <- apply(test.Pred.ens[[py]], 1:2, sd)
            # fitted.Pred.sd <- t(matrix(apply(fitted.Pred.ens, 2, sd),
            #                            n = data$n, ncol = Nt))
            if (nrow(test.Res.miss.index[[py]]) > 0) {
              cat("Missing for test datasets will be updated via JSTVC ...\n")
              test[[py]]$Y_ts[test.Res.miss.index[[py]][, 1], test.Res.miss.index[[py]][, 2]] <-
                ifelse(test.Pred[2, test.Res.miss.index[[py]][, 1], test.Res.miss.index[[py]][, 2]] >
                         max(data[[py]]$Y_ts, na.rm = TRUE),
                       median(data[[py]]$Y_ts, na.rm = TRUE),
                       test.Pred[2, test.Res.miss.index[[py]][, 1], test.Res.miss.index[[py]][, 2]]
                )
            }
            # plot(test.Pred[2,,])
            test.spT <- spT_validation(z        = test[[py]]$Y_ts,
                                       zhat     = test.Pred[2,,],
                                       sigma    = test.Pred.sd,#,
                                       zhat.Ens = test.Pred.ens[[py]],
                                       names    = F,
                                       CC       = F)

            Y.true <- as.vector(test[[py]]$Y_ts)
            Y.pred <- as.vector(test.Pred[2,,])
            true.da <- as.vector(test[[py]]$Y_ts)
            inx <- which(true.da > 0)
            if(plot){
              plot(true.da, as.vector(test.Pred[2,,]), pch = 19, cex = 0.5,
                   xlab = "True values", ylab = "Predicted values",
                   xlim = c(min(true.da, as.vector(test.Pred[2,,]), na.rm = T), max(true.da, as.vector(test.Pred[2,,]), na.rm = T)),
                   ylim = c(min(true.da, as.vector(test.Pred[2,,]), na.rm = T), max(true.da, as.vector(test.Pred[2,,]), na.rm = T)))
            }

            test.spT <- as.vector(test.spT)
            RMSE.2   <- test.spT[1]
            Coef.2   <- test.spT[2]
            FAC2.2   <- test.spT[3]
            CRPS.2   <- test.spT[4]
            MAE.2    <- test.spT[5]

            if(verbose){cat(paste0("Testing object [", Object, "]\n"))
              cat(paste0("Mean absolute error (MAE) [", Round(MAE.2, 4),"] \n"
                         , "Root mean squared error (RMSE) [", RMSE.2, "] \n"
                         , "Continuous rank probability score (CRPS) [", CRPS.2, "] \n"
                         , "Fraction of predictions within a factor of two (FAC2) [", FAC2.2,"] \n"
                         , "Pearson correlation between the predictions and obs [", Coef.2,"] \n\n"
              ))
            }
          }
        }
      }
    }

    # plot ----
    if((iter == 0) |(iter %%1 == 0) |(iter == (itMax - 1))|(isFALSE(criterion))){
      Fitted.RF <- Pred.RF <- list()
      for(Pg in 1:(PG)){
        py.inx <- G.Mat$Res.indic[[names(G.Mat$data.G.Mat)[Pg]]]
        Fitted.slope.Hv <- sapply(seq_len(Nt), function(t) {

          ini.IEnKF.Input$train.H.basis[[py.inx]] %*% slope.stRF.ens[[names(G.Mat$data.G.Mat)[Pg]]]$Vt.mean[t + 1, 1:data[[py.inx]]$nKnots]

        }, simplify = "matrix") %>% t()
        # for(py in 1:Py){
        Fitted.RF[[names(data)[py.inx]]] <- Fitted.slope.Hv
        # }
        if (!is.null(test)) {
          Pred.slope.Hv <- sapply(seq_len(Nt), function(t) {
            # G.Mat$test.G.Mat[[names(G.Mat$test.G.Mat)[Pg]]][t,,]
            ini.IEnKF.Input$test.H.basis[[py.inx]] %*% slope.stRF.ens[[names(G.Mat$test.G.Mat)[Pg]]]$Vt.mean[t + 1, 1:data[[py.inx]]$nKnots]
          }, simplify = "matrix") %>% t()
          # for(py in 1:Py){
          Pred.RF[[names(data)[py.inx]]] <- Pred.slope.Hv
          # }
        }
      }


      for(py in 1:Py){
        colnames(Fitted.RF[[names(data)[py]]]) <- colnames(data[[py]]$Y_ts)
      }
      if (!is.null(test)) {
        for(py in 1:Py){
          colnames(Pred.RF[[names(data)[py]]]) <- colnames(test[[py]]$Y_ts)
        }
      }
      for(py in 1:Py){
        # G.names <- c(dimnames(data[[py]]$sG_ts)[[1]])
        G.names <- paste0(names(data)[[py]], ".", dimnames(data[[py]]$sG_ts)[[1]])
        if(!is.null(G.names)){
          G.Ind <- which(names(G.Mat$data.G.Mat) %in% G.names)
          for(Pg in G.Ind){
            cat(paste0("Fitted GMRF wt(s) in the [", names(G.Mat$data.G.Mat)[Pg], "] of ", names(data)[py], ": "),
                Round(range(Fitted.RF[[names(data)[py]]]), 3), "\n")
            if (!is.null(test)) {
              cat(paste0("Predicted GMRF wt(s) in the [", names(G.Mat$test.G.Mat)[Pg], "] of ", names(data)[py], ": "),
                  Round(range(Pred.RF[[names(data)[py]]]), 3), "\n")
            }
            cat(paste0("...\n"))
          }
        }
      }

      all.slope.Hv <- Fitted.slope.Hv <- Pred.slope.Hv <- rf.plot.data <- list()
      for(py in 1:Py){
        Fitted.slope.Hv[[names(data)[py]]]     <- cbind(fitting.data.frame[[py]], intercept.GRF = as.vector(Fitted.RF[[py]]))

        setDT(Fitted.slope.Hv[[names(data)[py]]]);
        Fitted.slope.Hv[[names(data)[py]]]$Type <- "Fitted GRFs"


        if (!is.null(test)) {
          Pred.slope.Hv[[names(data)[py]]] <- cbind(pred.data.frame[[py]], intercept.GRF = as.vector(Pred.RF[[py]]))
          setDT(Pred.slope.Hv[[names(data)[py]]])
          Pred.slope.Hv[[names(data)[py]]]$Type <- "Interpolated GRFs"
          all.slope.Hv[[names(data)[py]]] <- rbind(Fitted.slope.Hv[[names(data)[py]]], Pred.slope.Hv[[names(data)[py]]])
        }else{
          all.slope.Hv[[names(data)[py]]] <- Fitted.slope.Hv[[names(data)[py]]]
        }

      }
    }
    Para.List <- VB.Laplace$Para.List
    {
      # Print results ----
      para.new <- cbind(c(iter + 1, round(c(beta,
                                            mu.sigma.sq,
                                            show.list$rho.v$mu.rho.v,
                                            show.list$Phi.v$mu.Phi.v[1],
                                            show.list$proc.tau.sq$mu.tau.sq,
                                            # show.list$proc.ini.tau.sq$mu.tau.sq,
                                            (end_time - start_time)[[3]]), 4), Object)) %>% t() %>% as.data.frame() %>% setDT()

      true.Para.List <- NULL
      true.temp <- NULL
      current.para <- rbind(true.temp, para.new)
      colnames(para.new) <- colnames(current.para) <- c("Iter", beta.name,
                                                        obs.sigma.sq.name,
                                                        paste0("rho.v_", 1:data[[1]]$Grid.infor$summary$res),
                                                        paste0("Phi.v_", 1:data[[1]]$Grid.infor$summary$res),
                                                        paste0("tau.sq_", 1:data[[1]]$Grid.infor$summary$res),
                                                        # paste0("ini.tau.sq_", 1:data[[1]]$Grid.infor$summary$res),
                                                        "elapsed", "Object")

      cat("************************************************************\n")
      cat("* \n")
      cat("* \n")
      cat("*              Currently estimated parameters               \n")
      cat("* \n")
      cat("* \n")
      cat("************************************************************\n")
      if(data[[1]]$Grid.infor$summary$res > 1){
        index.1 <- which(colnames(para.new) %in% obs.sigma.sq.name)
        print(current.para[1:nrow(para.new), 1:index.1])

        # rho.v
        index.2 <- which(colnames(para.new) %in% "Phi.v_1")
        cat("\n")
        print(current.para[1:nrow(para.new), (index.1 + 1): (index.2 - 1)])


        # zeta0_1
        cat("\n")
        index.5 <- which(colnames(para.new) %in% "tau.sq_1")
        print(current.para[1:nrow(para.new), (index.4): (index.5 - 1)])

        #proc.ini.tau.sq_1
        cat("\n")
        index.7 <- which(colnames(para.new) %in% "elapsed")
        print(current.para[1:nrow(para.new), (index.5): (index.7 - 1)])

        #other
        cat("\n")
        index.8 <- which(colnames(para.new) %in% "Object")
        print(current.para[1:nrow(para.new), (index.7): (index.8)])
      }else{
        out <- as.data.frame(current.para)
        Name <- colnames(out)
        inx <- as.vector(which(sapply(out, anyNA)))
        if(length(inx) >= 1){
          print(as.data.table(out[, -inx]))
        }else{
          print(as.data.table(out))
        }
      }

      cat(paste0("••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• \n"))
      cat(paste0("***************************************************************** \n"))
      cat(paste0("\n Iteration: ", iter + 1))
      cat(paste0("\n Relative error of log(likelihood) in the data model [",
                 Round(ELBO.error, 6), "]\n\n"))
      cat(paste0("***************************************************************** \n"))
      cat(paste0("••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••••• \n"))

    }

    Relative.ELBO.error <- c(Relative.ELBO.error, ELBO.error)


    # Saving results ----
    if((iter == 0) | (save.Predict) | (iter == (itMax - 1))|(isFALSE(criterion))){
      Tab.name <- c("iter", "Fitting_RMSE", "Testing_RMSE",
                    obs.sigma.sq.name,
                    "ELBO_Error",
                    "ELBO",
                    beta.name,
                    paste0("rho_v", 1:data[[1]]$Grid.infor$summary$res),
                    paste0("phi_v", 1:data[[1]]$Grid.infor$summary$res),
                    paste0("Proc_tau_sq_", 1:data[[1]]$Grid.infor$summary$res),
                    "Fitting_CRPS", "Testing_CRPS",
                    "Fitting_MAE", "Testing_MAE",
                    "Corr_Before", "Corr_After",
                    "Fitting_FAC2", "Testing_FAC2",
                    "elapsed", "run_time", "Object")
      {

        if(iter == 0){
          true.temp <- NULL

          if(is.null(VB.Laplace$Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]])){
            ini.show.list <-  ini.Para.List[[names(data)[[1]]]]$st.sRF[[names(G.Mat$data.G.Mat)[1]]]
          }else{
            ini.show.list <- ini.Para.List[[names(data)[[1]]]]$st.RF[[names(G.Mat$data.G.Mat)[1]]]
          }
          temp <- cbind(rbind(true.temp,
                              c(0, rep(NA, 2),
                                mu.sigma.sq,
                                rep(NA, 1),    ELBO.old,
                                beta,
                                ini.show.list$rho.v$mu.rho.v,
                                ini.show.list$Phi.v$mu.Phi.v[1],
                                ini.show.list$proc.tau.sq$mu.tau.sq,
                                # ini.show.list$proc.ini.tau.sq$mu.tau.sq,
                                cs, ct, rep(NA, 7)),
                              c(iter + 1, RMSE.1, RMSE.2,
                                mu.sigma.sq,
                                ELBO.error,
                                ELBO.new,
                                beta,
                                show.list$rho.v$mu.rho.v,
                                show.list$Phi.v$mu.Phi.v[1],
                                show.list$proc.tau.sq$mu.tau.sq,
                                # show.list$proc.ini.tau.sq$mu.tau.sq,
                                CRPS.1, CRPS.2,
                                MAE.1, MAE.2,
                                Coef.1, Coef.2, FAC2.1, FAC2.2,
                                round((end_time - start_time)[[3]], 2))) %>% as.data.frame(),
                        run_time, Object)
          colnames(temp) <- Tab.name
          # RowName[1] <- "True:"
          RowName[1] <-  "Initi:"
          RowName[2] <-  "Iter1:"
          temp[, 2:(ncol(temp) - 3)] <- round(temp[, 2:(ncol(temp) - 3)], 5)

          details.est <- NULL
        }else{
          temp0 <- cbind(rbind(NULL, c(iter + 1, RMSE.1, RMSE.2,
                                       mu.sigma.sq,
                                       ELBO.error,
                                       ELBO.new,
                                       beta,
                                       show.list$rho.v$mu.rho.v,
                                       show.list$Phi.v$mu.Phi.v[1],
                                       show.list$proc.tau.sq$mu.tau.sq,
                                       # show.list$proc.ini.tau.sq$mu.tau.sq,
                                       CRPS.1, CRPS.2,
                                       MAE.1, MAE.2,
                                       Coef.1, Coef.2, FAC2.1, FAC2.2,
                                       round((end_time - start_time)[[3]], 2)))%>% as.data.frame(),
                         run_time, Object)
          colnames(temp0) <- Tab.name

          temp0[, 2:(ncol(temp0) - 3)] <- round(temp0[, 2:(ncol(temp0) - 3)], 5)

          details.est <- rbind(details.est, temp0)
          RowName[iter + 3] <- paste0("Iter", iter + 1, ":")
        }
      }
    }
    ELBO.old <- ELBO.new

    indx <- c(which(colnames(all.slope.Hv[[py]]) %in% c("time.index", "LON", "LAT")),
              ncol(all.slope.Hv[[py]]) - 1, ncol(all.slope.Hv[[py]]))

    if(MCMC){
      if(iter > mcmc.random.effect.burnin){
        # }else{
        current.random.effect.draw <- compute_current_fitted_random_effect_draw()
        for (region_name in names(current.random.effect.draw)) {
          mcmc.random.effect.draws[[region_name]][[length(mcmc.random.effect.draws[[region_name]]) + 1]] <-
            current.random.effect.draw[[region_name]]
        }
        mcmc.random.effect.iterations <- c(mcmc.random.effect.iterations, iter)

        if((exists("sim_Data"))){
          for(py in 1:Py){
            temp <- all.slope.Hv[[py]][all.slope.Hv[[py]]$Type == "Fitted GRFs", ] #c(3, 4, 5, 17, 18)
            setDF(temp)
            temp <- temp[, indx]
            setDT(temp)
            temp$iter <- iter
            GRFs.process <- rbind(GRFs.process, temp)
          }
        }
      }
    }
    iter <- iter + 1
    # save(vb.EnKS.Iter,
    #      file = paste0(Tab, "/Result/", real, "_smoothed_JSTVC_mcmc_", MCMC, ".RData"))
  }

  if((exists("sim_Data")) & (MCMC)){
    GRFs.process <- GRFs.process %>% left_join(sim_Data[, c(5:6, 8, 12)], by = c("LON", "LAT", "time.index"))

  }

  mcmc.random.effect.history <- NULL
  if (MCMC && length(mcmc.random.effect.iterations) > 0) {
    mcmc.random.effect.history <- lapply(mcmc.random.effect.draws, function(draw_list) {
      draw_array <- simplify2array(draw_list)
      if (length(dim(draw_array)) == 2) {
        draw_array <- array(
          draw_array,
          dim = c(nrow(draw_list[[1]]), ncol(draw_list[[1]]), 1),
          dimnames = list(
            rownames(draw_list[[1]]),
            colnames(draw_list[[1]]),
            paste0("Iter_", mcmc.random.effect.iterations[1])
          )
        )
      } else {
        dimnames(draw_array)[[3]] <- paste0("Iter_", mcmc.random.effect.iterations)
      }
      draw_array
    })
  }



  if(!MCMC){
    # mcmc.random.effect.iterations <- NULL
    for(py in 1:Py){
      temp <- all.slope.Hv[[py]][all.slope.Hv[[py]]$Type == "Fitted GRFs", ]  #c(3, 4, 5, 17, 18)
      setDF(temp)
      temp <- temp[, indx]
      setDT(temp)
      temp$iter <- iter
      GRFs.process <- rbind(GRFs.process, temp)
    }
  }
  # save(vb.EnKS.Iter, GRFs.process,
  #      file = paste0(Tab, "/Result/", real, "_JSTVC_smoothed_mcmc_", MCMC, ".RData"))

  for(py in 1:Py){
    data[[py]]$Y_ts <- true.Y_ts[[py]]
  }
  if (!is.null(test)) {
    Re <- list(data                   = data,
               test.data              = test,
               ini.Para.List          = ini.Para.List,
               update.Para.List       = VB.Laplace$Para.List,
               detailed.Para.List     = details.est,
               augEnKS                = slope.stRF.ens,
               all.slope.Hv           = all.slope.Hv,
               fitted.slope.Hv.Zg.ens = fitted.slope.Hv.Zg.ens,
               fitted.Xts.ens         = fitted.Xts.ens,
               Pred.slope.Hv.Zg.ens   = Pred.slope.Hv.Zg.ens,
               test.Xts.ens           = test.Xts.ens,
               fitted.Pred.ens        = fitted.Pred.ens,
               test.Pred.ens          = test.Pred.ens,
               data.G.Mat             = G.Mat,
               rf.plot.dat            = rf.plot.data,
               mean.y                 = mean.y,
               sd.y                   = sd.y,
               Goodness.of.fit        = Goodness.of.fit,
               Testing_validation     = test.spT,
               Tuning.parameter       = list(cs = cs, ct = ct, Ne = Ne),
               Process.monitoring     = list(ELBO  = ELBO.new,
                                             DIC            = Goodness.of.fit$Overall$DIC,
                                             WAIC           = Goodness.of.fit$Overall$WAIC,
                                             tol.real       = tol.real,
                                             iteration      = iter,
                                             itMin          = itMin,
                                             itMax          = itMax,
                                             cores          = n.cores,
                                             VB.EnKS.Iter = as.data.table(vb.EnKS.Iter),
                                             MCMC.random.effect.iterations = mcmc.random.effect.iterations,
                                             MCMC.random.effect.history = mcmc.random.effect.history,
                                             GRFs.process   = GRFs.process,
                                             Relative.ELBO.error = Relative.ELBO.error),
               Object                   = Object)
  }else{
    Re <- list(
      data                   = data,
      ini.Para.List          = ini.Para.List,
      update.Para.List       = VB.Laplace$Para.List,
      detailed.Para.List     = details.est,
      augEnKS                = slope.stRF.ens,
      all.slope.Hv           = all.slope.Hv,
      fitted.slope.Hv.Zg.ens = fitted.slope.Hv.Zg.ens,
      fitted.Pred.ens        = fitted.Pred.ens,
      data.G.Mat             = G.Mat,
      rf.plot.dat            = rf.plot.data,
      mean.y                 = mean.y,
      sd.y                   = sd.y,
      Goodness.of.fit        = Goodness.of.fit,
      Tuning.parameter       = list(cs = cs, ct = ct, Ne = Ne),
      Process.monitoring     = list(ELBO  = ELBO.new,
                                    DIC            = Goodness.of.fit$Overall$DIC,
                                    WAIC           = Goodness.of.fit$Overall$WAIC,
                                    tol.real       = tol.real,
                                    iteration      = iter,
                                    itMin          = itMin,
                                    itMax          = itMax,
                                    cores          = n.cores,
                                    VB.EnKS.Iter   = as.data.table(vb.EnKS.Iter),
                                    MCMC.random.effect.iterations = mcmc.random.effect.iterations,
                                    MCMC.random.effect.history = mcmc.random.effect.history,
                                    GRFs.process   = GRFs.process,
                                    Relative.ELBO.error = Relative.ELBO.error),
      Object                = Object)
  }
  return(Re)
}
