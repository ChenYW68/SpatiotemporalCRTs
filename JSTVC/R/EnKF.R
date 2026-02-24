EnKF <- function(data,
                       yts.minus.xts,
                       Para.List,
                       H,
                       Mv,
                       spTaper,
                       slope.G.mat    = FALSE,
                       ct             = 1,
                       Ne             = 100,
                       var.name       = "Intercept",
                       mu.sigma.sq    = 1,
                       inv.mu.sigma.sq = 1,
                       Py             = 0){
  # cat("***************************************************************** \n")
  # cat("                Start to execute EnKF! \n\n")
  # cat("*-----------------------------------------------------------\n")
  cat("* ---               \n")

  {
    spKnots.num <- data[[Py]]$Grid.infor$summary$nKnots
    ind.mesh    <- c(1:data[[Py]]$nKnots)
    Nt <- data[[Py]]$Nt
    n <- length(mu.sigma.sq)


    Sigma.R  <- spam_diag(n)*mu.sigma.sq
    Sigma.sR <- spam_diag(n)*inv.mu.sigma.sq

    schol.Sigma.sR <- spam::chol.spam(Sigma.sR)


    xf <- xp <- array(NA, dim   =  c(Nt + 1, spKnots.num, Ne),
                      dimnames = list(c(as.character(as.numeric(dimnames(data[[Py]]$Y_ts)[[1]][1]) - 1),
                                        dimnames(data[[Py]]$Y_ts)[[1]]),
                                      c(paste0("spKnots.", 1:spKnots.num)),
                                      as.character(paste0("Ens.", 1:Ne))))



    rho.time <- taper(0:Nt, ct + 2)
    names(rho.time) <- paste(0:Nt)

    for(g in 1:data[[Py]]$Grid.infor$summary$res){
      E <- matrix(rnorm(data[[Py]]$Grid.infor$level[[g]]$nKnots*Ne),
                  nrow = data[[Py]]$Grid.infor$level[[g]]$nKnots,
                  ncol = Ne)
      g.index <- data[[Py]]$Grid.infor$summary$g.index[[g]]

      {
        mu.alpha <- Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$LR.coef$mu.alpha

        mat.x <- as.matrix(data[[Py]]$Grid.infor$var.covariable)
        mu.var <- exp(mat.x %*% mu.alpha)
        rf.corr <-
          exp(-data[[Py]]$Grid.infor$level[[g]]$BAUs.Dist/Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$Phi.v$mu.Phi.v[1]
              )

         rf.corr <-  diag(as.vector(mu.var^(0.5))) %*% rf.corr %*% diag(as.vector(mu.var^(0.5)))
        Q <- spam::as.spam(Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$proc.tau.sq$mu.tau.sq[g]*solve(rf.corr))


        xf[1, g.index, ] <- spam::backsolve(as.spam(as.matrix(Matrix::chol(as.dgCMatrix.spam(Q)))), E)
      }
    }
  }
  ###############################################################

  if(Nt >= 1){
    pb <- progress::progress_bar$new(format = "EnKF: (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                     , total = Nt
                                     , clear = FALSE
                                     # , width= Nt
    )
    L_Q_upper <- list()
    for(g in 1:data[[Py]]$Grid.infor$summary$res){
      g.index <- data[[Py]]$Grid.infor$summary$g.index[[g]]
      {
        mu.alpha <- Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$LR.coef$mu.alpha

        mat.x <- as.matrix(data[[Py]]$Grid.infor$var.covariable)

        mu.var <- exp(mat.x %*% mu.alpha)
        rf.corr <-
          exp(-data[[Py]]$Grid.infor$level[[g]]$BAUs.Dist/Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$Phi.v$mu.Phi.v[1]
                    )
         rf.corr <-  diag(as.vector(mu.var^(0.5))) %*% rf.corr %*% diag(as.vector(mu.var^(0.5)))
        Q <- spam::as.spam(Para.List[[Py]]$st.sRF[[var.name]]$proc.tau.sq$mu.tau.sq[g]*solve(rf.corr))

        L_Q_upper[[g]] <- as.spam(as.matrix(Matrix::chol(as.dgCMatrix.spam(Q))))
      }

    }
    H.map <- H
    for (t in 2:(Nt + 1)){
      pb$tick()
      Sys.sleep(1 / (2*Nt))
      if(slope.G.mat){
        H.map <- spam::as.spam(H[t - 1,,])
      }

      ##########################################################################
      xp[max(1, t - ct - 3):(t - 1), , ] <- xf[max(1, t - ct - 3):(t - 1), , ]
      ##########################################################################
      ##########################################################################
      # if(data$Grid.infor$summary$res > 1){
      for(g in 1:data[[Py]]$Grid.infor$summary$res){
        E  <- matrix(rnorm(data[[Py]]$Grid.infor$level[[g]]$nKnots*Ne),
                     nrow = data[[Py]]$Grid.infor$level[[g]]$nKnots,
                     ncol = Ne)
        xi <- spam::backsolve(L_Q_upper[[g]], E)
        if(is.list(Mv)){
          xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ] <- xi +
            as.matrix(Mv[[g]]) %*% xf[t - 1, data[[Py]]$Grid.infor$summary$g.index[[g]], ]
        }else{
          xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ] <-
            as.matrix(Mv[g] * xf[t - 1, data[[Py]]$Grid.infor$summary$g.index[[g]], ] + xi)
        }
      }

      ##########################################################################
      e <- spam::backsolve(schol.Sigma.sR, matrix(rnorm(n*Ne), ncol = Ne))

      ##########################################################################
      yp          <- H.map %*% xp[t, 1:data[[Py]]$nKnots, ] + e

      rm(E, xi)

      {

        Phat <-  spam::as.spam(cov(t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[1]], ]),
                    t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[1]], ]))*spTaper$taper[[1]])

        if(data[[Py]]$Grid.infor$summary$res > 1){
          for(g in 2:data[[Py]]$Grid.infor$summary$res){
            A <- spam::as.spam(cov(t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ]),
                                   t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ]))*spTaper$taper[[g]])
            Phat <- spam::bdiag.spam(Phat, A)

          }
        }

        HP      <- H.map %*% Phat %*% spam::t.spam(H.map) + (Sigma.R)
        Chol.HP <- spam::chol.spam(HP)
        Sk      <- spam::t.spam(H.map) %*% spam::as.spam(spam::backsolve(Chol.HP, spam::forwardsolve(Chol.HP, yts.minus.xts[t - 1,] - yp)))
        Sk.mean <- spam::t.spam(H.map) %*% spam::as.spam(spam::backsolve(Chol.HP, spam::forwardsolve(Chol.HP, yts.minus.xts[t - 1,] - rowMeans(yp))))
        P.tH    <- Phat %*% t(H.map)
      }

      #  update step:
      for (j in max(1, t - ct):t){
        Phat <- Matrix::Diagonal(0)
        for(g in 1:data[[Py]]$Grid.infor$summary$res){
          A <- (cov(t(xp[j, data[[Py]]$Grid.infor$summary$g.index[[g]], ]),
                    t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ]))* as.numeric(rho.time[paste(t - j)]) *spTaper$taper[[g]])

          Phat <- Matrix::bdiag(Phat, spam::as.dgCMatrix.spam(A))
        }
        xf[j, , ] <- xp[j, , ] + spam::as.spam.dgCMatrix(Phat) %*% Sk
       }
    }
  }
  return(
    list(
      Vt.mean  = apply(xf, 1:2, mean),
      Vt.ens   = xf,
      rho.time = rho.time
    ))
}
