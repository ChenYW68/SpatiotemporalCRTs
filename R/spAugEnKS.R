Aug.spEnKS <- function(data,
                       yts.minus.xts,
                       Para.List,
                       H,
                       Mv,
                       spTaper,
                       Z.tsv          = NULL,
                       sZ.tsv         = NULL,
                       slope.G.mat    = FALSE,
                       ct             = 1,
                       Ne             = 100,
                       var.name       = "Intercept",
                       mu.sigma.sq    = 1,
                       Py             = 0){
  # cat("***************************************************************** \n")
  # cat("                Start to execute TaEnKs algorithm ! \n\n")
  # cat("*-----------------------------------------------------------\n")
  cat("* ---               \n")

  {
    spKnots.num <- data[[Py]]$Grid.infor$summary$nKnots
    ind.mesh    <- c(1:data[[Py]]$nKnots)

    n <- length(mu.sigma.sq)

    Sigma.R  <- spam_diag(n)*mu.sigma.sq
    Sigma.sR <- spam_diag(n)/mu.sigma.sq

    schol.Sigma.sR <- spam::chol.spam(Sigma.sR)
    Pz <- 0; gt <- NULL
    Nt <- data[[Py]]$Nt

    if(!is.null(Z.tsv)){
      Pz <- data[[Py]]$Pz
      for(i in 1:Pz){
        gt <- c(gt, paste0("Public.gt.", dimnames(data[[Py]]$Z_ts)[[1]][i]))
      }
    }

    if(!is.null(sZ.tsv)){
      for(py in 1:Py){
        sPz <- data[[py]]$sPz
        for(i in 1:sPz){
          gt <- c(gt, paste0(names(data)[py], ".gt.", dimnames(data[[py]]$Z_ts)[[1]][i]))
        }
      }
    }
    if(length(gt) > 0){
      #--- prior.ensemble
      xf <- xp <- array(NA, dim   =  c(data[[Py]]$Nt + 1, spKnots.num + length(gt), Ne),
                        dimnames = list(c(as.character(as.numeric(dimnames(data[[Py]]$Y_ts)[[1]][1]) - 1),
                                          dimnames(data[[Py]]$Y_ts)[[1]]),
                                        c(paste0("spKnots.", 1:spKnots.num), gt),
                                        as.character(paste0("Ens.", 1:Ne))))
    }else{
      xf <- xp <- array(NA, dim   =  c(data[[Py]]$Nt + 1, spKnots.num, Ne),
                        dimnames = list(c(as.character(as.numeric(dimnames(data[[Py]]$Y_ts)[[1]][1]) - 1),
                                          dimnames(data[[Py]]$Y_ts)[[1]]),
                                        c(paste0("spKnots.", 1:spKnots.num)),
                                        as.character(paste0("Ens.", 1:Ne))))
    }

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
          exp(-sqrt(data[[Py]]$Grid.infor$level[[g]]$BAUs.Dist^2/Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$Phi.v$mu.Phi.v[1]^2))

         rf.corr <-  diag(as.vector(mu.var^(0.5))) %*% rf.corr %*% diag(as.vector(mu.var^(0.5)))
        Q <- spam::as.spam(Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$proc.tau.sq$mu.tau.sq[g]*solve(rf.corr)) #
        xf[1, g.index, ] <- spam::backsolve(as.spam(as.matrix(Matrix::chol(as.dgCMatrix.spam(Q)))), E)
      }
    }
    if(!is.null(Z.tsv)){
      ind <- (data[[Py]]$Grid.infor$summary$nKnots + 1 ):(Pz + data[[Py]]$Grid.infor$summary$nKnots)
      E   <- matrix(rnorm(Pz*Ne), nrow = Pz, ncol = Ne)
      xf[1, ind, ] <- spam::backsolve(chol(solve(Para.List[[1]]$ini.alpha.tau.sq$mu.ini.alpha.tau.sq[1:Pz]*diag(Pz))), E)
    }

    if(!is.null(sZ.tsv)){
      for(py in 1:Py){
        if(py == 1){
          sPz.s <- data[[Py]]$Grid.infor$summary$nKnots + Pz + 1
          sPz   <- 0
          ind   <- (sPz.s):(dim(data[[py]]$sZ_ts)[[1]] + data[[Py]]$Grid.infor$summary$nKnots + Pz)
        }else{
          ind   <- (max(ind) + 1):(sPz + dim(data[[py]]$sZ_ts)[[1]] + data[[Py]]$Grid.infor$summary$nKnots + Pz)
        }
        sPz          <- sPz + dim(data[[py]]$sZ_ts)[[1]]
        E            <- matrix(rnorm(dim(data[[py]]$sZ_ts)[[1]]*Ne), nrow = dim(data[[py]]$sZ_ts)[[1]], ncol = Ne)
        xf[1, ind, ] <-
          spam::backsolve(chol(solve(Para.List[[names(data)[py]]]$self.ini.alpha.tau.sq$mu.ini.alpha.tau.sq[1:dim(data[[py]]$sZ_ts)[[1]]]*
                                       diag(dim(data[[py]]$sZ_ts)[[1]]))), E)
      }

    }
  }

  ###############################################################
  if(Nt >= 1){
    pb <- progress::progress_bar$new(format = "spAugEnKS: (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                     , total = Nt
                                     , clear = FALSE
    )
    L_Q_upper <- list()
    for(g in 1:data[[Py]]$Grid.infor$summary$res){
      g.index <- data[[Py]]$Grid.infor$summary$g.index[[g]]
      {
        mu.alpha <- Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$LR.coef$mu.alpha

        mat.x <- as.matrix(data[[Py]]$Grid.infor$var.covariable)
        mu.var <- exp(mat.x %*% mu.alpha)
        rf.corr <-
          exp(-sqrt(data[[Py]]$Grid.infor$level[[g]]$BAUs.Dist^2/Para.List[[names(data)[Py]]]$st.sRF[[var.name]]$Phi.v$mu.Phi.v[1]^2))

         rf.corr <-  diag(as.vector(mu.var^(0.5))) %*% rf.corr %*% diag(as.vector(mu.var^(0.5)))
        Q <- spam::as.spam(Para.List[[Py]]$st.sRF[[var.name]]$proc.tau.sq$mu.tau.sq[g]*solve(rf.corr)) #
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


      if(!is.null(Z.tsv)){
        ind          <- (data[[Py]]$Grid.infor$summary$nKnots + 1):(Pz + data[[Py]]$Grid.infor$summary$nKnots)
        E            <- matrix(rnorm(Pz*Ne), nrow = Pz, ncol = Ne)
        xp[t, ind, ] <- as.matrix(Para.List[[1]]$rho.alpha$mu.rho.alpha[1:Pz]*diag(Pz) %*% xf[t - 1, ind, ]) +
          spam::backsolve(chol(solve(diag(Pz)*Para.List[[1]]$alpha.tau.sq$mu.alpha.tau.sq[1:Pz, 1])), E)
      }

      if(!is.null(sZ.tsv)){
        for(py in 1:Py){
          if(py == 1){
            sPz.s <- data[[Py]]$Grid.infor$summary$nKnots + Pz + 1
            sPz   <- 0
            ind   <- (sPz.s):(dim(data[[py]]$sZ_ts)[[1]] + data[[Py]]$Grid.infor$summary$nKnots + Pz)
          }else{
            ind <- (max(ind) + 1):(sPz + dim(data[[py]]$sZ_ts)[[1]] + data[[Py]]$Grid.infor$summary$nKnots + Pz)
          }
          sPz          <- sPz + dim(data[[py]]$sZ_ts)[[1]]
          E            <- matrix(rnorm(dim(data[[py]]$sZ_ts)[[1]]*Ne), nrow = dim(data[[py]]$sZ_ts)[[1]], ncol = Ne)
          xp[t, ind, ] <- as.matrix(Para.List[[names(data)[py]]]$self.rho.alpha$mu.rho.alpha[1:dim(data[[py]]$sZ_ts)[[1]], 1]*
                                      diag(dim(data[[py]]$sZ_ts)[[1]]) %*% xf[t - 1, ind, ]) +
            spam::backsolve(chol(solve(diag(dim(data[[py]]$sZ_ts)[[1]])*
                                         Para.List[[names(data)[py]]]$self.alpha.tau.sq$mu.alpha.tau.sq[1:dim(data[[py]]$sZ_ts)[[1]], 1])), E)
        }
      }
      ##########################################################################
      e <- spam::backsolve(schol.Sigma.sR, matrix(rnorm(n*Ne), ncol = Ne))
      ##########################################################################
      yp          <- H.map %*% xp[t, 1:data[[Py]]$nKnots, ] + e
      if(!is.null(Z.tsv)){
        ind       <- (data[[Py]]$Grid.infor$summary$nKnots + 1):(Pz + data[[Py]]$Grid.infor$summary$nKnots)
        public.gt <- matrix(xp[t, ind, ], nrow = Pz)
        yp        <- yp + Z.tsv[[t - 1]] %*% public.gt
      }
      if(!is.null(sZ.tsv)){
        for(py in 1:Py){
          if(py == 1){
            sPz   <- 0
            sPz.s <- data[[Py]]$Grid.infor$summary$nKnots + Pz + 1
            ind   <- (sPz.s):(dim(data[[py]]$sZ_ts)[[1]] + data[[Py]]$Grid.infor$summary$nKnots + Pz)
          }else{
            ind   <- (max(ind) + 1):(sPz + dim(data[[py]]$sZ_ts)[[1]] + data[[Py]]$Grid.infor$summary$nKnots + Pz)
          }
          z.ind   <- (sPz + 1):(sPz + dim(data[[py]]$sZ_ts)[[1]])
          sPz     <- sPz + dim(data[[py]]$sZ_ts)[[1]]
          self.gt <- matrix(xp[t, ind, ], nrow = dim(data[[py]]$sZ_ts)[[1]])

          if(length(z.ind) == 1){
            yp <- yp + matrix(as.matrix(sZ.tsv[[t - 1]][, z.ind]), ncol = dim(data[[py]]$sZ_ts)[[1]]) %*% self.gt
          }else{
            yp <- yp + sZ.tsv[[t - 1]][, z.ind] %*% self.gt
          }
        }
      }
      rm(E, e,  xi)

      { 
        Phat <-  spam::as.spam.dgCMatrix(cov(t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[1]], ]),
                    t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[1]], ]))*spTaper$taper[[1]])
        if(data[[Py]]$Grid.infor$summary$res > 1){
          for(g in 2:data[[Py]]$Grid.infor$summary$res){
            A <- spam::as.spam(cov(t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ]),
                                   t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ]))*spTaper$taper[[g]])
            Phat <- spam::bdiag.spam(Phat, A)
          }
        }
        if(is.null(Z.tsv)&is.null(sZ.tsv)){
          HP      <- H.map %*% Phat %*% spam::t.spam(H.map) + (Sigma.R)
          Chol.HP <- spam::chol.spam(HP)
          Sk      <- spam::t.spam(H.map) %*% spam::as.spam(spam::backsolve(Chol.HP, spam::forwardsolve(Chol.HP, yts.minus.xts[t - 1,] - yp)))
        }else{
          if((!is.null(Z.tsv))&(is.null(sZ.tsv))){
            HP      <- cbind(H.map, Z.tsv[[t - 1]]) %*% spam::bdiag.spam(Phat, spam::as.spam(cov(t(matrix(xp[t, -ind.mesh, ], nrow = Pz))))) %*%
              spam::t.spam(cbind(H.map, Z.tsv[[t - 1]])) + (Sigma.R)
            Chol.HP <- spam::chol.spam(HP)
            Sk      <- spam::as.spam(spam::t.spam(cbind(H.map, Z.tsv[[t - 1]]))) %*%
              spam::as.spam(spam::backsolve(Chol.HP, spam::forwardsolve(Chol.HP, yts.minus.xts[t - 1,] - yp)))
          }
          if((is.null(Z.tsv))&(!is.null(sZ.tsv))){
            sPz <- ncol(sZ.tsv[[1]])
            HP  <- cbind(H.map, sZ.tsv[[t - 1]]) %*% spam::bdiag.spam(Phat, spam::as.spam(cov(t(matrix(xp[t, -ind.mesh, ], nrow = sPz))))) %*%
              spam::t.spam(cbind(H.map, sZ.tsv[[t - 1]])) + (Sigma.R)
            Chol.HP <- spam::chol.spam(HP)
            Sk      <- spam::as.spam(spam::t.spam(cbind(H.map, sZ.tsv[[t - 1]]))) %*%
              spam::as.spam(spam::backsolve(Chol.HP, spam::forwardsolve(Chol.HP, yts.minus.xts[t - 1,] - yp)))
          }
          if((!is.null(Z.tsv))&(!is.null(sZ.tsv))){
            sPz <- ncol(sZ.tsv[[1]])
            HP  <- cbind(H.map, Z.tsv[[t - 1]], sZ.tsv[[t - 1]]) %*%
              spam::bdiag.spam(Phat, spam::as.spam(cov(t(matrix(xp[t, -ind.mesh, ], nrow = Pz + sPz))))) %*%
              spam::t.spam(cbind(H.map, Z.tsv[[t - 1]], sZ.tsv[[t - 1]])) + (Sigma.R)
            Chol.HP <- spam::chol.spam(HP)
            Sk      <- spam::as.spam(spam::t.spam(cbind(H.map, Z.tsv[[t - 1]], sZ.tsv[[t - 1]]))) %*%
              spam::as.spam(spam::backsolve(Chol.HP, spam::forwardsolve(Chol.HP, yts.minus.xts[t - 1,] - yp)))
          }
        }
        rm(HP, Phat, Chol.HP)
      }

      #  update step:
      for (j in max(1, t - ct):t){
        Phat <- Matrix::Diagonal(0)
        for(g in 1:data[[Py]]$Grid.infor$summary$res){
          A <- (cov(t(xp[j, data[[Py]]$Grid.infor$summary$g.index[[g]], ]),
                    t(xp[t, data[[Py]]$Grid.infor$summary$g.index[[g]], ]))* as.numeric(rho.time[paste(t - j)]) *spTaper$taper[[g]])

          Phat <- Matrix::bdiag(Phat,A)
        }

        if((!is.null(Z.tsv))&(is.null(sZ.tsv))){
          Phat <- Matrix::bdiag(Phat, cov(t(matrix(xp[j, -ind.mesh, ], nrow = Pz)),
                                          t(matrix(xp[t, -ind.mesh, ], nrow = Pz))))

        }else if((is.null(Z.tsv))&(!is.null(sZ.tsv))){
          sPz  <- ncol(sZ.tsv[[1]])
          Phat <- Matrix::bdiag(Phat, cov(t(matrix(xp[j, -ind.mesh, ], nrow = sPz)),
                                          t(matrix(xp[t, -ind.mesh, ], nrow = sPz))))
        }else if((!is.null(Z.tsv))&(!is.null(sZ.tsv))){
          sPz  <- ncol(sZ.tsv[[1]])
          ind  <- (data[[Py]]$Grid.infor$summary$nKnots + 1):(Pz + data[[Py]]$Grid.infor$summary$nKnots)
          Phat <- Matrix::bdiag(Phat, cov(t(matrix(xp[j, ind, ], nrow = Pz)),
                                          t(matrix(xp[t, ind, ], nrow = Pz))),
                                cov(t(matrix(xp[j, -c(ind.mesh, ind), ], nrow = sPz)),
                                    t(matrix(xp[t, -c(ind.mesh, ind), ], nrow = sPz))))
        }
        xf[j, , ] <- xp[j, , ] + spam::as.spam.dgCMatrix(Phat) %*% Sk
      }
    }
  }

  Vt.mean <- apply(xf, 1:2, mean)
  rm(xp, L_Q_upper, Sk);
  return(
    list(
      Vt.mean  = Vt.mean,
      Vt.ens   = xf,
      rho.time = rho.time
    ))
}
