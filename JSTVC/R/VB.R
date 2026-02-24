.VB <- function (data            = NULL,
                 Hv.Zg           = 0,
                 data.trans.tsv  = NULL,
                 Ks              = NULL,
                 S               = NULL,
                 prior           = NULL,
                 Para.List       = NULL,
                 verbose.VB      = FALSE,
                 test            = NULL,
                 Ne              = 100,
                 iter            = 0){
  options(warn = -1)
  cat(paste0("\n***************************************************************** \n"))
  cat(paste0("*\n"))
  cat(paste0("*     Start to execute VB ... \n"))
  cat(paste0("*\n"))
  cat(paste0("***************************************************************** \n"))
  n  <- 0
  Nt <- data[[names(data)[[1]]]]$Nt
  assign("sub.Nt", Nt, envir = .GlobalEnv)


  Py <- length(names(data))
  residuals <- residuals.x <- data.trans.tsv$tsv.y - Hv.Zg
  k <- 0
  ww <- 1#(iter + 1)^(-0.5)

  for(py in 1:Py){ n <- n + data[[names(data)[[py]]]]$n }

  #  beta ----
  if(!is.null(data[[names(data)[[1]]]]$X_ts)){
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
    beta.name             <- dimnames(data[[names(data)[[1]]]]$X_ts)[[1]]
    X_CuSum               <- XYXi_CuSum <- 0
    E_inverse_sigma.sq    <- NULL
    # Yts <- data.trans.tsv$tsv.y
    for(py in 1:Py){
      E_inverse_sigma.sq <- c(E_inverse_sigma.sq,
                              rep(Para.List[[names(data)[[py]]]]$obs.sigma.sq$a/Para.List[[names(data)[[py]]]]$obs.sigma.sq$b,
                                  data[[names(data)[[py]]]]$n))
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


    # ensemble for beta
     w <- array(NA, dim  =  c(Nt, n, Ne),
                 dimnames = list(1:Nt, 1:n, as.character(paste0("Ens.", 1:Ne))))
     for(py in 1:Py){
       for(t in 1:Nt){
         w[t, data[[names(data)[[py]]]]$index.y,] <-  data[[names(data)[[py]]]]$Hs %*%
           Ks[[paste0(names(data)[[py]], ".Intercept")]]$Vt.ens[t + 1,,]
       }
     }
     XYXi_CuSum.ens <- matrix(0, nrow = data[[1]]$Px, ncol = Ne)
     for(e in 1:Ne){
      for (t in 1:Nt) {
       X_ts <- matrix(as.matrix(data.trans.tsv$tsv.x[[t]]), nrow = dim(data[[names(data)[[1]]]]$X_ts)[[1]], ncol = n)
       rownames(X_ts) <- dimnames(data[[names(data)[[1]]]]$X_ts)[[1]]
       XYXi_CuSum.ens[, e]     <- XYXi_CuSum.ens[, e] + X_ts %*% ((data.trans.tsv$tsv.y[t,] - w[t,, e]) * E_inverse_sigma.sq)
       }
     }
  }

# A <- apply( w, 1:2, mean)



  if(!is.null(data[[names(data)[[1]]]]$X_ts)){
    post_betaX_sigma.sq <- solve(X_CuSum + solve(prior[[names(data)[[1]]]]$beta$sigma.sq))

    eta <- XYXi_CuSum + solve(prior[[names(data)[[1]]]]$beta$sigma.sq) %*% prior[[names(data)[[1]]]]$beta$mu.beta

    post_betaX_mu     <- post_betaX_sigma.sq %*% eta %>% as.vector()
    post_betaX_mu.ens <- post_betaX_sigma.sq %*% XYXi_CuSum.ens

    Sigma_beta_LR <- cov(t(post_betaX_mu.ens), t(post_betaX_mu.ens))


    for(py in 1:Py){
      Para.List[[names(data)[[py]]]]$beta$mu.beta           <- matrix(Para.List[[names(data)[[py]]]]$beta$mu.beta[, 1], nrow = length(beta.name), ncol = 1)
      rownames(Para.List[[names(data)[[py]]]]$beta$mu.beta) <- beta.name
      colnames(Para.List[[names(data)[[py]]]]$beta$mu.beta) <- "Sharing coffficients"

      old.mu <- Para.List[[names(data)[[py]]]]$beta$mu.beta[, 1]
      old.cov <- ifelse(is.null(Para.List[[names(data)[[py]]]]$beta$sigma.sq[3]),
                        post_betaX_sigma.sq, Para.List[[names(data)[[py]]]]$beta$sigma.sq)


      Para.List[[names(data)[[py]]]]$beta$mu.beta[, 1]      <-  (1 - ww)*old.mu + ww*post_betaX_mu
      Para.List[[names(data)[[py]]]]$beta$sigma.sq          <- (1 - ww)*old.cov + ww*(post_betaX_sigma.sq + Sigma_beta_LR)
    }


    cov.sd <- data.frame(round(as.data.frame(Para.List[[names(data)[[1]]]]$beta$sigma.sq), 5), "|",
                         SD = round(sqrt(as.vector(diag(Para.List[[names(data)[[1]]]]$beta$sigma.sq))), 5), "|",
                         old.SD = round(sqrt(as.vector(diag(as.matrix(post_betaX_sigma.sq)))), 5))
    colnames(cov.sd)[c(ncol(cov.sd) - 3, ncol(cov.sd) - 1)] <- c("*", "**")
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
 mean.residuals <- data.trans.tsv$tsv.y - Xts.beta.fun(tsv.x = data.trans.tsv$tsv.x,
                                                       alpha = Para.List[[names(data)[[1]]]]$beta$mu.beta,
                                                       self  = FALSE)
  t1 <- proc.time()
  #  sigma.sq ----
  n.diff <- NULL
  for(py in 1:length(names(data))){

    old.a <- Para.List[[names(data)[[py]]]]$obs.sigma.sq$a
    old.b <- Para.List[[names(data)[[py]]]]$obs.sigma.sq$b


    Para.List[[names(data)[[py]]]]$obs.sigma.sq$a <- a_sigma.sq <-
     (prior[[names(data)[[py]]]]$obs.sigma.sq$a + data[[names(data)[[py]]]]$n * data[[names(data)[[py]]]]$Nt/2)



    temp1 <-  sapply(seq_len(Nt), function(t) {
      (t(data[[names(data)[[py]]]]$Y_ts[t, ]) %*% data[[names(data)[[py]]]]$Y_ts[t, ] - 2 * t(Para.List[[names(data)[[py]]]]$beta$mu.beta)%*%
         (data[[names(data)[[py]]]]$X_ts[, t, ]) %*% (data[[names(data)[[py]]]]$Y_ts[t, ]) )[1] +
        sum(diag((data[[names(data)[[py]]]]$X_ts[, t, ] %*% t(data[[names(data)[[py]]]]$X_ts[, t, ]) %*%
                    (Para.List[[names(data)[[py]]]]$beta$sigma.sq + Para.List[[names(data)[[py]]]]$beta$mu.beta %*% t(Para.List[[names(data)[[py]]]]$beta$mu.beta)))))
    }, simplify = "matrix") %>% sum()


    # temp1 <-  sapply(seq_len(Nt), function(t) {
    #   (t(residuals.x[t, data[[names(data)[[py]]]]$index.y]) %*% residuals.x[t, data[[names(data)[[py]]]]$index.y] - 2 * t(Para.List[[names(data)[[1]]]]$beta$mu.beta)%*%
    #      (data[[names(data)[[py]]]]$X_ts[, t, ]) %*% (residuals.x[t, data[[names(data)[[py]]]]$index.y]) )[1] +
    #     sum(diag((data[[names(data)[[py]]]]$X_ts[, t, ] %*% t(data[[names(data)[[py]]]]$X_ts[, t, ]) %*%
    #                 (Para.List[[names(data)[[1]]]]$beta$sigma.sq + Para.List[[names(data)[[1]]]]$beta$mu.beta %*% t(Para.List[[names(data)[[1]]]]$beta$mu.beta)))))
    # }, simplify = "matrix") %>% sum()


    temp2 <- sapply(seq_len(Nt), function(t) {
      -(2 * spam::t(Hv.Zg[t, data[[names(data)[[py]]]]$index.y]) %*% mean.residuals[t, data[[names(data)[[py]]]]$index.y])[1]
      # -(2 * spam::t(apply(xi[t, data[[names(data)[[py]]]]$index.y,], 1, mean)) %*%
      #     mean.residuals[t, data[[names(data)[[py]]]]$index.y])[1]
    }, simplify = "matrix") %>% sum()




    temp3 <-  sum(spam::diag.spam((t(data[[names(data)[[py]]]]$Hs) %*% data[[names(data)[[py]]]]$Hs %*% S[[paste0(names(data)[[py]], ".Intercept")]]$S11)))

    f_sigma.sq <- temp1 + temp2 + temp3

    if((f_sigma.sq < 0)){
       f_sigma.sq  <-  t(as.vector(residuals[, data[[names(data)[[py]]]]$index.y])) %*% as.vector(residuals[, data[[names(data)[[py]]]]$index.y])
    }

     Para.List[[names(data)[[py]]]]$obs.sigma.sq$b <- b_sigma.sq <-

     (prior[[names(data)[[py]]]]$obs.sigma.sq$b + f_sigma.sq/2)



    # if(VB){
    # cat("obs.sigma.sq---", b_sigma.sq, a_sigma.sq, "\n")

    old.sigma.sq  <-  Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq
    ww0 <- (iter + 1)^(-0.5)

    Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq <- (1 - ww0) * old.sigma.sq + ww0*ifelse(is.na(as.vector(b_sigma.sq/(a_sigma.sq - 1))), old.sigma.sq,
                                                                                                   as.vector(b_sigma.sq/(a_sigma.sq - 1)))
#
#     if(Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq > 1e2){
#       # cat("**************************", "\n")
#       # cat(paste0("obs.sigma.sq [", round(Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq, 5), "]\n"))
#       Para.List[[names(data)[[py]]]]$obs.sigma.sq$mu.sigma.sq <- old.sigma.sq#var(as.vector(data[[names(data)[[py]]]]$Y_ts))
#     }

    if (verbose.VB) {
      cat("---------------------------------", "\n")
      # cat("*********************************", "\n")
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
            cat(paste0("* GMRFs (self) on the ", g, "th subdomain ...\n"))  # with respect to the [", sub.rf.name, "] of ", names(data)[py], "
            # cat("* \n")
            cat("************************************************************\n\n")
          }

          {
            g.index <- data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]
            #---other parameters

            rf.corr <-
              exp(-data[[py]]$Grid.infor$level[[g]]$BAUs.Dist/Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1]# +
                          # data[[py]]$Grid.infor$level[[g]]$var.dist^2/Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[2]^2
                        )


            # rf.corr <- fields::Wendland(data[[py]]$Grid.infor$level[[g]]$BAUs.Dist,
            #                  theta = Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1],
            #                  dimension = 1,
            #                  k = 1)



            mu.alpha <-  Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$LR.coef$mu.alpha
            assign("var.covariable", data[[names(data)[[py]]]]$Grid.infor$var.covariable, envir = .GlobalEnv)
            #

            mat.x <- as.matrix(var.covariable)
            # mu.var <- vector()
            # for(l in 1:nrow(mat.x)){
            #   mu.var[l] <- as.vector(mat.x[l, ] %*% mu.alpha)^2
            # }
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


            QS.Sum <- (sum(spam::diag.spam(Q %*% S[[sub.rf.name]]$S11[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]],
                                                                  data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]])) -
                     2*Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g] * Q.St_1_0 +
                     (Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$var.rho.v[g] +
                        Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$rho.v$mu.rho.v[g]^2) * Q.St_1_1)/2 +
              sum(spam::diag.spam((Q) %*% S[[sub.rf.name]]$V.t00[data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]],data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]]))/2

            # rate
            # Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g] <- prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g] + data[[names(data)[[py]]]]$Grid.infor$level[[g]]$nKnots * (Nt + 1)/2
            #
            #
            #  Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g] <- prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g] + QS.Sum
            #
            # Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g]/(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g])


            # tau.sq <- range.tau.sq[1] + diff(range.tau.sq) / (1 + exp(-eta))
            #  tau.sq ----
            Low <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$L[g]#1e-5
            Upp <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$H[g]#1e-2

            assign("sub.n", data[[names(data)[[py]]]]$n, envir = .GlobalEnv)
            # assign("sub.Nt", data[[names(data)[[py]]]]$Nt, envir = .GlobalEnv)
            assign("Low", Low, envir = .GlobalEnv)
            assign("Upp", Upp, envir = .GlobalEnv)
            assign("prior.tau.a", prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g], envir = .GlobalEnv)
            assign("prior.tau.b", prior[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g], envir = .GlobalEnv)
            assign("QS.Sum", QS.Sum, envir = .GlobalEnv)

            if(is.null(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g])){
              Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g] <- 2
            }
            if(is.null(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g])){
              Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g] <- 1
            }
            op.tau.sq <- optim(c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g]),
                                 log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g])), #
                             fn = non_f_tau.sq.kl,
                             method  = "L-BFGS-B", #L-BFGS-B
                             lower   = c(log(1e-2), log(1e-2)),
                             upper   = c(log(1e2), log(1e2)),
                             control = list(maxit = 100,
                                            factr = 1e0,  # Less precise convergence
                                            lmm = 5,
                                            trace = verbose.VB)     # Smaller history
            )


            # op.tau.sq <- optim(c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g]),
            #                      log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g])),
            #                    non_f_tau.sq.kl,
            #                    lower   = c(log(1e-3)),
            #                    upper   = c(log(1e3)))

            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g] <- exp(op.tau.sq$par[1]) + 2.1
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g] <- exp(op.tau.sq$par[2])
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] <- Low + (Upp - Low)/(1 + (exp(op.tau.sq$par[2]))/(exp(op.tau.sq$par[1]) + 1))


            # range.tau.sq <- c(L,  H)
            # old.tau.sq <-  Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.ini.tau.sq$mu.tau.sq[g]
            #
            # Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.ini.tau.sq$mu.tau.sq[g] <-
            #     Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] <-
            #   (1 - ww)*old.tau.sq + ww*
            #   ifelse((Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g]   < range.tau.sq[1])|
            #            (Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] > range.tau.sq[2]),
            #          previous.est, Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g])



            # previous.est <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g]

            # Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g] <- 1
            if(verbose.VB){
              cat(paste0("Expectation of the tau.sq in the ",
                         g, "th subdomain [",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$a[g], 3), "; ",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$b[g], 3), "; ",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g], 5), "]\n"))
              # cat("............................s............................\n")
            }
            assign("mu.tau.sq", Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$proc.tau.sq$mu.tau.sq[g], envir = .GlobalEnv)
            # assign("Nt", Nt, envir = .GlobalEnv)
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



            assign("residuals", mean.residuals, envir = .GlobalEnv) #mean.residuals
            assign("Hdist", data[[names(data)[[py]]]]$Hdist[, data[[names(data)[[py]]]]$Grid.infor$summary$g.index[[g]]], envir = .GlobalEnv)

            assign("sigma.sq", Para.List$common.obs.sigma.sq$mu.sigma.sq[py, py], envir = .GlobalEnv)
            assign("Vt.mean", Ks[[sub.rf.name]]$Vt.mean[-1, ], envir = .GlobalEnv)
            assign("index.y", data[[names(data)[[py]]]]$index.y, envir = .GlobalEnv)



            # op.var_alpha <- nlminb(c(mu.alpha),
            #                   non_f_var_alpha, #optim#pso::psoptim#nlminb
            #                   # method = "BFGS",
            #                   lower   = c(1e-5, -10, -10, -10, -10),
            #                   upper   = c(5, 5, 5, 5, 5, 5),
            #                   control = list(trace = verbose.VB, iter.max = 30))
            #  alpha ----
            op.var_delta <- optim(c(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$LR.coef$mu.alpha),
                                   fn = non_f_var_alpha, #optim#pso::psoptim#nlminb
                                   method = "L-BFGS-B", #Nelder-Mead
                                   lower   = c(-1e0),
                                   upper   = c(5e0),
                                   control = list(trace = verbose.VB
                                                  # , iter.max = 1e2
                                                  ))

            # op.var_delta$par <- c(0, 0)
            # # print(op.var_delta)
            # # op.var_delta$par <- c(1.5, -1, -0.1)
            # hessc <- numDeriv::hessian(func = non_f_var_alpha, x = op.var_delta$par)
            #
            #
            # # Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$var.Phi.v[g] < 50
            # Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$LR.coef$var.alpha <- ifelse(1/hessc > 0, 1/hessc, 1e-6)
            old.alpha <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$LR.coef$mu.alpha


            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$LR.coef$mu.alpha <- (1- ww)*old.alpha + ww*op.var_delta$par
            # exp(op.var_delta$par + diag(Para.List[[names(data)[py]]]$st.sRF[[sub.rf.name]]$LR.coef$var.alpha)/2)
            assign("mu.alpha", op.var_delta$par, envir = .GlobalEnv)


            if(verbose.VB){
              # print(paste0("Expectation of the delta in the ",
              #            g, "th subdomain [", Round(op.var_delta$par, 5), "]\n"))
              cat(paste0("Expectation of the delta in the ",
                     g, "th subdomain [",
                     paste0(sprintf("%.4f", op.var_delta$pa), collapse = ", "),
                     "]\n"))
              # cat("............................s............................\n")
            }


            #  zeta ----
            op.zeta <- optim(c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[g])),
                             fn = non_f_range, #optim#pso::psoptim#nlminb
                              method = "L-BFGS-B",
                              lower   = c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$L[g])),
                              upper   = c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$H[g]))#, log(1e2)
                           )


            # if(is.null(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$a)){
            #   Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$a <- 5
            # }
            # if(is.null(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$b)){
            #   Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$b <- 0.1
            # }
            #
            #
            # op.zeta <- optim(c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$a),
            #                    log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$b)),
            #                  fn = non_f_range.kl, #optim#pso::psoptim#nlminb
            #                  method = "L-BFGS-B",
            #                  lower   = c(log(1e-3)),
            #                  upper   = c(log(1e2))#, log(1e2)
            # )


            # op.zeta
           # op.zeta <- nlminb(c(log(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[g])),
           #                    non_f_theta, #optim#pso::psoptim#nlminb
           #                    # method = "BFGS",
           #                    lower   = log(1e0),
           #                    upper   = log(50),
           #                    control = list(trace = verbose.VB, iter.max = 30)) verbose.VB, iter.max = 30))

            # print(op.zeta)
            hessc <- numDeriv::hessian(func = non_f_range, x = op.zeta$par)


            mu.phi <- (exp(op.zeta$par[1]))#/(exp(op.zeta$par[2]))


            old.phi <- Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[g]

            ww0 <-  1
            Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[g] <- (1- ww0)*old.phi + ww0*mu.phi


             rr <- 1#
             Hs.temp <- exp(-(rr*Hdist/Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[g]))

             gg <- Hs.temp/apply(Hs.temp, 1, sum)
             data[[names(data)[[py]]]]$Hs <- gg

             if(!is.null(test)){
               Hs.temp <-  exp(-(rr*test[[names(test)[[py]]]]$Hdist[, test[[names(test)[[py]]]]$Grid.infor$summary$g.index[[g]]]/Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[g]))

               gg <- Hs.temp/apply(Hs.temp, 1, sum)
               test[[names(test)[[py]]]]$Hs <- gg
             }

            if(verbose.VB){
              cat(paste0("Expectation of the Phi.v in the ",
                         g, "th subdomain [",
                         # Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$a, 3), "; ",
                         # Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$b, 3), "; ",
                         Round(Para.List[[names(data)[[py]]]]$st.sRF[[sub.rf.name]]$Phi.v$mu.Phi.v[1], 3), "]\n"))

            }
          }
        }
      }
    }
  }
  return(list(Para.List = Para.List, prior = prior, data = data, test = test))
}


