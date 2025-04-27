.VB.Laplace <- function (data              = NULL,
                         Hv.Zg           = 0,
                         data.trans.tsv  = NULL,
                         Ks              = NULL,
                         S               = NULL,
                         prior           = NULL,
                         Para.List       = NULL,
                         verbose.VB      = FALSE,
                         nu              = c(1e-1, 1e-1, 1e-1),
                         var.select      = TRUE,
                         threshold       = c(2, 2),
                         test            = test,
                         iter            = 0){
  options(warn = -1)
  cat(paste0("\n***************************************************************** \n"))
  cat(paste0("*\n"))
  cat(paste0("*     Start to execute the VB and Laplace approximate ... \n"))
  cat(paste0("*\n"))
  cat(paste0("***************************************************************** \n"))
  n  <- 0
  Nt <- data[[names(data)[[1]]]]$Nt
  Py <- length(names(data))
  residuals <- residuals.x <- data.trans.tsv$tsv.y - Hv.Zg
  k <- 0

  for(py in 1:Py){ n <- n + data[[names(data)[[py]]]]$n }

  #  beta ----
  Sigma.alpha <- NULL
  if(!is.null(data[[names(data)[[1]]]]$X_ts)){
    #update varPhi
    if(is.null(Para.List[[names(data)[[1]]]]$beta$gamma[1])){
      Para.List[[names(data)[[1]]]]$beta$gamma <- vector()
      for(px in 1:data[[names(data)[[1]]]]$Px){
        Para.List[[names(data)[[1]]]]$beta$gamma[px] <- 1
      }
    }
    if(is.null(Para.List[[names(data)[[1]]]]$beta$sigma.sq[1])){
      Para.List[[names(data)[[1]]]]$beta$sigma.sq <- vector()
      for(px in 1:data[[names(data)[[1]]]]$Px){
        Para.List[[names(data)[[1]]]]$beta$sigma.sq[px] <- 1/1E4
      }
    }
    #beta distribution
    a  <- 1; b <- 1

    X.names <- dimnames(data[[names(data)[[1]]]]$X_ts)[[1]]

    Sigma.alpha <- prior[[names(data)[[1]]]]$beta$sigma.sq

    # if(is.na(Para.List[[names(data)[[1]]]]$beta$mu.beta[1])){
    beta.name             <- dimnames(data[[names(data)[[1]]]]$X_ts)[[1]]
    X_CuSum               <- XYXi_CuSum <- 0
    E_inverse_sigma.sq    <- NULL
    # Yts <- data.trans.tsv$tsv.y
    for(py in 1:Py){
      E_inverse_sigma.sq <- c(E_inverse_sigma.sq, rep(1/Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq, data[[names(data)[[py]]]]$n))
      if(!is.null(data[[names(data)[[py]]]]$sX_ts)){
        residuals.x[, data[[names(data)[[py]]]]$index.y]  <- residuals[, data[[names(data)[[py]]]]$index.y] -
          Xts.beta.fun(tsv.x = data[[names(data)[[py]]]]$sX_ts,
                       alpha = Para.List[[names(data)[[py]]]]$sbeta$mu.beta,
                       self  = TRUE) #- Hv.Zg[, data[[names(data)[[py]]]]$index.y]
      }
    }

    for (t in 1:Nt) {
      X_ts <- matrix(as.matrix(data.trans.tsv$tsv.x[[t]]), nrow = dim(data[[names(data)[[1]]]]$X_ts)[[1]], ncol = n)
      rownames(X_ts) <- dimnames(data[[names(data)[[1]]]]$X_ts)[[1]]
      X_CuSum        <- X_CuSum + X_ts %*% (t(X_ts) * E_inverse_sigma.sq)
      XYXi_CuSum     <- XYXi_CuSum + X_ts %*% (residuals.x[t,] * E_inverse_sigma.sq)
    }
  }


  if(!is.null(data[[names(data)[[1]]]]$X_ts)){
    post_betaX_sigma.sq <- solve(X_CuSum + solve(Sigma.alpha))
    if(var.select){eta <- XYXi_CuSum}else{
      eta <- XYXi_CuSum + solve(Sigma.alpha) %*% prior[[names(data)[[1]]]]$beta$mu.beta
    }
    post_betaX_mu <- post_betaX_sigma.sq %*% eta %>% as.vector()
    for(py in 1:Py){
      Para.List[[names(data)[[py]]]]$beta$mu.beta           <- matrix(NA, nrow = length(beta.name), ncol = 1)
      rownames(Para.List[[names(data)[[py]]]]$beta$mu.beta) <- beta.name
      colnames(Para.List[[names(data)[[py]]]]$beta$mu.beta) <- "Sharing coffficients"
      Para.List[[names(data)[[py]]]]$beta$mu.beta[, 1]      <- post_betaX_mu
      Para.List[[names(data)[[py]]]]$beta$sigma.sq          <- post_betaX_sigma.sq

    }

    cov.sd <- data.frame(round(as.data.frame(post_betaX_sigma.sq), 5), "|",
                         SD = round(sqrt(as.vector(diag(post_betaX_sigma.sq))), 5))
    colnames(cov.sd)[c(ncol(cov.sd) - 1)] <- c("|")
    Para.List[[names(data)[[1]]]]$beta$cov.prob <- cov.sd

    if (verbose.VB) {
      cat("\n.................................\n")
      cat(paste0("Posterior expectation of the beta (public): \n"))
      print(round(Para.List[[names(data)[[py]]]]$beta$mu.beta, 5))
      cat(paste0("Posterior variance-covariance of the beta (public): \n"))
      print(cov.sd)
    }
    # }
    residuals <- residuals - Xts.beta.fun(tsv.x = data.trans.tsv$tsv.x,
                                          alpha = Para.List[[names(data)[[1]]]]$beta$mu.beta,
                                          self  = FALSE)
  }

  # sigma.sq
  t1 <- proc.time()
  n.diff <- NULL
  for(py in 1:length(names(data))){
    Para.List[[names(data)[[py]]]]$obs.sigma.sq$a <- a_sigma.sq <- prior[[names(data)[[py]]]]$obs.sigma.sq$a + data[[names(data)[[py]]]]$n * data[[names(data)[[py]]]]$Nt/2
    f_sigma.sq                     <-  t(as.vector(residuals[, data[[names(data)[[py]]]]$index.y])) %*% as.vector(residuals[, data[[names(data)[[py]]]]$index.y])
    Para.List[[names(data)[[py]]]]$obs.sigma.sq$b <- b_sigma.sq <- prior[[names(data)[[py]]]]$obs.sigma.sq$b + f_sigma.sq/2
    Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq <- ifelse(is.na(as.vector(b_sigma.sq/(a_sigma.sq - 1))), 1,
                                                                      as.vector(b_sigma.sq/(a_sigma.sq - 1)))

    if(Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq > 1e4){
      # cat("**************************", "\n")
      Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq <- 1
    }

    if (verbose.VB) {
      cat("---------------------------------", "\n")
      cat(paste0("obs.sigma.sq of ", names(data)[[py]], " [", round(Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq, 5), "]\n\n"))

    }
    n.diff[py] <- data[[names(data)[[py]]]]$n
  }

  # if(corr.within.res)
  {
    if(length(diff(n.diff)) == 0){
      Para.List$common.obs.sigma.sq$mu.sigma.sq <-
        matrix(Para.List[[names(data)[[1]]]]$obs.sigma.sq$mu.sigma.sq, 1, 1)
    }else if(sum(diff(n.diff)) == 0){
      temp <- NULL
      for(py in 1:length(names(data))){
        # for(t in 1: data[[names(data)[[py]]]]$Nt){
        temp <- cbind(temp, (
          as.vector(residuals[, data[[names(data)[[py]]]]$index.y])
        ))
        # }
      }
      Para.List$common.obs.sigma.sq$mu.sigma.sq <- cov(temp)
    }else if(sum(diff(n.diff)) != 0){
      Para.List$common.obs.sigma.sq$mu.sigma.sq <- diag(length(names(data)))
      for (py in 1:Py){
        Para.List$common.obs.sigma.sq$mu.sigma.sq[py, py] <- Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq
      }
    }
  }

  # st.sRF ----
  for(py in 1:length(names(data))){
    if(!is.null(data[[names(data)[[py]]]]$sG_ts)){
      rf.name   <- dimnames(data[[names(data)[[py]]]]$sG_ts)[[1]]
      for(sPg in 1:dim(data[[names(data)[[py]]]]$sG_ts)[[1]]){
        sub.rf.name <- paste0(names(data)[[py]], ".", rf.name[sPg])
        cat("\n")
        cat("************************************************************\n")
        cat("* \n")
        cat(paste0("* To find parameters of dynamic GMRFs (self) in the [", paste0(names(data)[[py]], ".", rf.name[sPg]), "] of ", names(data)[py], "\n"))
        cat("* \n")
        # cat("*---\n")
        for(g in 1:data[[names(data)[[py]]]]$Grid.infor$summary$res){
          if (verbose.VB) {
            # cat("************************************************************\n")
            # cat("* \n")
            cat("*---\n")
            cat(paste0("* GMRFs (self) on the ", g, "th subdomain ...\n"))  
            # cat("* \n")
            cat("************************************************************\n\n")
          }

          {
            g.index <- data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]
            #---other parameters

            rf.corr <-
              exp(-sqrt(data[[py]]$Grid.infor$level[[g]]$BAUs.Dist^2/Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1]^2))

            mu.alpha <-  Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$LR.coef$mu.alpha
            assign("var.covariable", data[[names(data)[[py]]]]$Grid.infor$var.covariable, envir = .GlobalEnv)
            #

            mat.x <- as.matrix(var.covariable)
            mu.var <- exp(mat.x %*% mu.alpha)
            # cat("Variance: ", range(mu.var), "\n\n\n\n\n")
            rf.corr <- diag(as.vector(mu.var^(0.5))) %*% rf.corr %*% diag(as.vector(mu.var^(0.5)))
            Q <- spam::as.spam(solve(rf.corr))
            #----other parameters

            Q.St_1_0 <- sum(spam::diag.spam((Q) %*%S[[sub.rf.name]]$S10[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]],
                                                                        data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]]))

            Q.St_1_1 <- sum(spam::diag.spam((Q) %*%S[[sub.rf.name]]$S00[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]],
                                                                        data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]]))

            #
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g] <- solve(Q.St_1_1*Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq + 1/prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$sigma.sq[g]) %>%
              as.vector()
            theta1_eta <- Q.St_1_0*Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq
              prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu[g]/prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$sigma.sq[g]

            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g] %*% theta1_eta %>% as.vector()

            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] <- ifelse(abs(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g]) < 1e-10,
                                                                                             rnorm(1, mean = prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu[g],
                                                                                                   sd = sqrt(prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$sigma.sq[g])),
                                                                                             Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g])


            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] <- ifelse(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] > 1, 1,
                                                                                             Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g])
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] <- ifelse(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] < -1, -1,
                                                                                             Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g])


            if(verbose.VB){
              cat(paste0("Expectation of the rho.v in the ", g, "th subdomain [",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g], 5), "]\n"))
              cat(paste0("Variance of the rho.v in the ", g, "th subdomain [",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g], 5), "]\n\n"))
            }


            #--- tau.sq
            previous.est <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g]

            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g] <- prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g] + data[[names(data)[[py]]]]$Grid.infor$level[[g]]$nKnots * (Nt + 1)/2
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g] <- prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g] + (sum(spam::diag.spam(Q %*% S[[sub.rf.name]]$S11[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]],
                                                                                                                                                                                                         data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]])) -
                                                                                                                                                            2*Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] * Q.St_1_0 +
                                                                                                                                                            (Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g] +
                                                                                                                                                               Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g]^2) * Q.St_1_1)/2 +
              sum(spam::diag.spam((Q) %*% S[[sub.rf.name]]$V.t00[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]],data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]]))/2
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g]/(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g])


            range.tau.sq <- c(1e-2, 1E2)
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.ini.tau.sq$mu.tau.sq[g] <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] <-
              ifelse((Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g]   < range.tau.sq[1])|
                       (Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] > range.tau.sq[2]),
                     previous.est, Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g])
            if(verbose.VB){
              cat(paste0("Expectation of the tau.sq in the ",
                         g, "th subdomain [", Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g], 5), "]\n"))
              # cat("............................s............................\n")
            }
            assign("mu.tau.sq", Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g], envir = .GlobalEnv)
            assign("Nt", Nt, envir = .GlobalEnv)
            assign("Rho", Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g], envir = .GlobalEnv)
            assign("Rho.sq", (Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g] + Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g]^2), envir = .GlobalEnv)
            assign("S00", S[[sub.rf.name]]$S00[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]], envir = .GlobalEnv)
            assign("S10", S[[sub.rf.name]]$S10[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]], envir = .GlobalEnv)
            assign("S11", S[[sub.rf.name]]$S11[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]], envir = .GlobalEnv)
            assign("V.t00", S[[sub.rf.name]]$V.t00[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]], data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]], envir = .GlobalEnv)

            assign("G", data[[py]]$Grid.infor$level[[g]]$BAUs.Dist, envir = .GlobalEnv)

            assign("Adj.Mat.0",  data[[names(data)[[py]]]]$Grid.infor$level[[g]]$Adj.Mat, envir = .GlobalEnv)

            assign("var.dist", data[[py]]$Grid.infor$level[[g]]$var.dist, envir = .GlobalEnv)
            assign("mu.Phi.v", Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v, envir = .GlobalEnv)
            assign("prior.Phi.v.a", prior[[names(data)[py]]]$st.sRF[[sub.rf.name]]$Phi.v$a, envir = .GlobalEnv)
            assign("prior.Phi.v.b", prior[[names(data)[py]]]$st.sRF[[sub.rf.name]]$Phi.v$b, envir = .GlobalEnv)
            assign("LR.coef.prior", prior[[names(data)[py]]]$st.sRF[[sub.rf.name]]$LR.coef, envir = .GlobalEnv)
            assign("residuals", residuals, envir = .GlobalEnv)
            assign("Hdist", data[[names(data)[[py]]]]$Hdist[, data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]], envir = .GlobalEnv)
            assign("sigma.sq", Para.List$common.obs.sigma.sq$mu.sigma.sq[py, py], envir = .GlobalEnv)
            assign("Vt.mean", Ks[[sub.rf.name]]$Vt.mean[-1, ], envir = .GlobalEnv)
            assign("index.y", data[[names(data)[[py]]]]$index.y, envir = .GlobalEnv)
            op.var_delta <- optim(c(mu.alpha),
                                   fn = non_f_var_alpha,
                                   method = "L-BFGS-B",
                                   lower   = c(-1e2),
                                   upper   = c(1e0),
                                   control = list(trace = verbose.VB, iter.max = 30))
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$LR.coef$mu.alpha <- op.var_delta$par
            assign("mu.alpha", op.var_delta$par, envir = .GlobalEnv)
            if(verbose.VB){
              cat(paste0("Expectation of the delta in the ",
                         g, "th subdomain [", Round(op.var_delta$par, 5), "]\n"))
            }
            # Hs ----
            op.zeta <- optim(c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1])),
                              fn      = non_f_range, 
                              method  = "L-BFGS-B",
                              lower   = c(log(1e0)),
                              upper   = c(log(1e2))
                           )
            hessc <- numDeriv::hessian(func = non_f_range, x = op.zeta$par)
            mu.phi <- exp(op.zeta$par[1])
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1] <- mu.phi
            Hs.temp <- exp(-(Hdist/mu.phi))
            data[[names(data)[[py]]]]$Grid.infor$Hs <- as(Hs.temp/apply(Hs.temp, 1, sum), "sparseMatrix")
             if(!is.null(test)){
               Hs.temp <-  exp(-(test[[names(test)[[py]]]]$Hdist[, test[[names(test)[[py]]]]$Grid.infor$summary$g.index[[g]]]/mu.phi))
               test[[names(test)[[py]]]]$Grid.infor$Hs <- as(Hs.temp/apply(Hs.temp, 1, sum), "sparseMatrix")
             }

            if(verbose.VB){
              cat(paste0("Expectation of the Phi.v in the ",
                         g, "th subdomain [",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1], 3), "]\n"))
            }
          }
        }
      }
    }
  }
  return(list(Para.List = Para.List, prior = prior, data = data, test = test))
}


