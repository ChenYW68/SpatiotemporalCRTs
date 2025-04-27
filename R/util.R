spatial.taper <- function(cs, data){
  spTaper <- list()
  # D <- vector()
  if(length(cs) < length(data)){
    cs <- rep(cs, length(data))
  }
  for(py in 1:length(data)){
    spTaper[[py]] <- list()
    for(g in 1:data[[py]]$Grid.infor$summary$res) {
      spTaper[[py]]$tuning[g] <- cs[py] 
      row.min <- min(data[[py]]$Grid.infor$level[[g]]$BAUs.Dist)
      spTaper[[py]]$tuning[g]  <- ifelse(spTaper[[py]]$tuning[g] < row.min, row.min, spTaper[[py]]$tuning[g])
      spTaper[[py]]$taper[[g]] <- as(fields::Wendland(data[[py]]$Grid.infor$level[[g]]$BAUs.Dist
                                                , theta = spTaper$tuning[g]
                                                , dimension = 1, k = 1),
                               "dgCMatrix")
    }
  }

  return(spTaper)
}
# transform spatial coordinates by WGS84
spCoords.transform <- function(loc, col = c("LON", "LAT"),
                               colname = c("LON_X", "LAT_Y"),
                               method = 1){
  nc <- ncol(loc)
  if(method == 1){
    d <- loc[, col]
    coordinates(d) <- col
    proj4string(d) <- CRS("+proj=longlat +datum=WGS84")# WGS 84
    CRS.new <- CRS("+proj=utm +zone=51 ellps=WGS84")
    sp.coord.trans <- spTransform(d, CRS.new)


    loc <- cbind(loc, sp.coord.trans@coords[,1], sp.coord.trans@coords[, 2])
  }else{
    LCC <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
+y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
    # transforming latitude and longitude using the Lambert conformal projection
    sp.coord <- SpatialPoints(coords = loc[, col],
                              proj4string = CRS("+proj=longlat +ellps=WGS84"))

    sp.coord.trans <- spTransform(sp.coord,  CRS(LCC))
    #
    #   #  %>% setnames(c("LON", "LAT"), c("LON_X", "LAT_Y"))
    sp.coord.trans <- as.data.frame(sp.coord.trans)
    #   colnames(sp.coord.trans) <- colname
    # loc = as.data.frame(loc)
    loc <- cbind(loc, sp.coord.trans[, 1], sp.coord.trans[, 2])

  }

  # grid.coords$LON_X = sp.coord.trans@coords[,1]
  # grid.coords$LAT_Y = d.ch1903@coords[,2]
  colnames(loc)[(nc + 1):(nc + 2)] <- colname
  return(loc)
}

# solve quantile
quant <- function(x){quantile(x, prob=c(0.025, 0.5, 0.975))}


# log likelihood of the HDCM in the first-level
loglik <- function(data, Para.List, yts.minus.xts, slope.Hv.Zg, slope.Au.Rg = 0)
{

  Nt <- data[[py]]$Nt
  n <- 0
  Py <- length(names(data))
  for(py in 1:Py){
    n <- n + data[[py]]$n
  }

  mu.sigma.sq <- vector()
  for (py in 1:Py){
    mu.sigma.sq[data[[py]]$index.y] <- Para.List[[py]]$obs.sigma.sq$mu.sigma.sq
  }

  inv.mu.sigma.sq <- spam::diag.spam(1/mu.sigma.sq)
  mu.sigma.sq <-  spam::diag.spam(mu.sigma.sq)

  Y_ts <- yts.minus.xts - slope.Hv.Zg - slope.Au.Rg

  # ylikelihood
  Loglik.y.1 <- spam::determinant.spam(mu.sigma.sq)$modulus*n*Nt


  Loglik.y.2 <- sapply(seq_len(Nt), function(t){
    t(Y_ts[t, ]) %*% inv.mu.sigma.sq %*% Y_ts[t, ]#/Para.List$obs.sigma.sq$mu.sigma.sq
  }, simplify = "matrix") %>% sum()

  Loglik <- (- Loglik.y.1 - Loglik.y.2)/2
  return(as.vector(Loglik))
}


# tapering function
taper <- function(distance, radius) Wendland(distance, radius, 1, 1)

#
chunk <- function(x,n)
{
  split(x, factor(sort(rank(x)%%n)))
}

Xts.beta.fun <- function(tsv.x, alpha, self = FALSE){
  if(!self){
    Nt <- length(tsv.x)
    if(length(as.vector(alpha)) > 1){
      temp <- sapply(seq_len(Nt), function(t)
        matrix(as.vector(alpha),
               nrow = 1) %*% tsv.x[[t]]
        , simplify = "matrix") #%>% t()
    }else{
      temp <- sapply(seq_len(Nt), function(t)
        tsv.x[[t]]*as.vector(alpha)
        , simplify = "matrix")

    }
    x_ts <- matrix(NA, nrow = Nt, ncol = ncol(tsv.x[[1]]))
    for(t in 1:Nt){
      x_ts[t, ] <- as.vector(temp[[t]]);
    }
  }else{
    Nt <- dim(tsv.x)[[2]]
    if(length(as.vector(alpha)) > 1){
      x_ts <- sapply(seq_len(Nt), function(t)
        matrix(as.vector(alpha),
               nrow = 1) %*% tsv.x[, t, ]
        , simplify = "matrix") %>% t()
    }else{
      x_ts <- sapply(seq_len(Nt), function(t)
        tsv.x[, t, ]*as.vector(alpha)
        , simplify = "matrix") %>% t()

    }
  }
  return(x_ts)
}



mat_vector_multi <- function(mat, vec)
{
  Multi = apply(t(vec), 1, FUN = '%*%', as.matrix(gpuR::t(mat)))
  return(Multi)
}
Empical_Cov_taper <- function(x, y, Taper, sym = F, nThreads = 8){
  non_zero_index.triu <- t(Matrix::which(Matrix::triu(Taper) != 0, arr.ind = T))

  n <- nrow(x)
  e <- ncol(x)
  m <- nrow(y)
  l <- ncol(non_zero_index.triu)
  mu_x <- rowMeans(x)
  mu_y <- rowMeans(y)

  storage.mode(n) <- "integer"
  storage.mode(m) <- "integer"
  storage.mode(e) <- "integer"
  storage.mode(l) <- "integer"
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  storage.mode(non_zero_index.triu) <- "integer"
  storage.mode(mu_x) <- "double"
  storage.mode(mu_y) <- "double"

  storage.mode(nThreads) <- "integer"


  cov.triu <- Empical_Cov_taper_c(x, y, mu_x, mu_y, non_zero_index.triu,
                                  n, m, e, l, nThreads)$R[, 1]

  # Emp_cov.1 <- spam(0, n, m)
  if(sym){
    # tril(Emp_cov) <- cov.triu
    index <- which(non_zero_index.triu[1, ] != non_zero_index.triu[2,])
    non_zero_index.triu_2 <- non_zero_index.triu[, index]
    cov.triu.2 <- cov.triu[index]
    Emp_cov <- Matrix::sparseMatrix(i = c(non_zero_index.triu[1, ],
                                          non_zero_index.triu_2[2, ]),
                                    j = c(non_zero_index.triu[2, ],
                                          non_zero_index.triu_2[1, ]),
                                    x = c(cov.triu, cov.triu.2))*Taper
  }else{
    non_zero_index.tril <- t(Matrix::which(Matrix::tril(Taper, -1) != 0, arr.ind = T))
    l <- ncol(non_zero_index.tril)
    storage.mode(non_zero_index.tril) <- "integer"
    storage.mode(l) <- "integer"

    cov.tril <- Empical_Cov_taper_c(x, y, mu_x, mu_y, non_zero_index.tril,
                                    n, m, e, l, nThreads)$R[, 1]
    Emp_cov <- Matrix::sparseMatrix(i = c(non_zero_index.triu[1, ],
                                          non_zero_index.tril[1, ]),
                                    j = c(non_zero_index.triu[2, ],
                                          non_zero_index.tril[2, ]),
                                    x = c(cov.triu, cov.tril))*Taper
  }
  return(Emp_cov)
}


tranform.list.to.matrix <- function(data, E_inv.sigma.sq = 1, E_inv.obs.delta = 1){


  name.y <- names(data)
  Py <- length(name.y)
  Nt <- data[[py]]$Nt
  n <- 0
  for(py in 1:Py){
    n <- n + data[[py]]$n
  }

  tsv.y <- matrix(NA, nrow = Nt, ncol = n)


  tsv.x <- tsv.sx <- tsv.z <- tsv.sz <- list()
  if(length(E_inv.obs.delta) < (Nt*Py)){
    E_inv.obs.delta <- matrix(E_inv.obs.delta, nrow = Nt, ncol = Py)
  }
  E_inv.sigma.sq <- E_inv.sigma.sq*E_inv.obs.delta

  for(t in 1:Nt) {
    y <- NULL
    for(py in 1:Py){
      y <- c(y, (data[[py]]$Y_ts[t, ]))
    }
    colnames(tsv.y) <- names(y)
    tsv.y[t, ] <- as.vector(y)



    # temp.sz <- Matrix::diag(0)

    temp.x <- temp.z <- temp.sz <- tem.X_ts <- tem.Z_ts <- tem.sZ_ts <- tem.sX_ts <-NULL
    for(py in 1:Py){
      if(!is.null(data[[py]]$X_ts)){
        if(data[[py]]$Px == 1){
          tem.X_ts <- matrix(data[[py]]$X_ts[, t, ], nrow = 1, ncol = data[[py]]$n)
        }else{
          tem.X_ts <- data[[py]]$X_ts[, t, ]
        }

      }
      temp.x <- cbind(temp.x, E_inv.sigma.sq[t, py] * tem.X_ts)
    }
    for(py in 1:Py){
      if(!is.null(data[[py]]$Z_ts)){
        if(data[[py]]$Pz == 1){
          tem.Z_ts <- matrix(data[[py]]$Z_ts[, t, ], nrow = data[[py]]$n, ncol = 1)
        }else{
          tem.Z_ts <- t(data[[py]]$Z_ts[, t, ])
        }

      }
      temp.z <- rbind(temp.z, tem.Z_ts)
    }


    for(py in 1:Py){
      if(!is.null(data[[py]]$sZ_ts)){
        if(data[[py]]$sPz == 1){
          tem.sZ_ts <- matrix(0, nrow = n, ncol = 1)
          tem.sZ_ts[data[[py]]$index.y, ] <- data[[py]]$sZ_ts[1, t, ]
        }else{
          tem.sZ_ts <- matrix(0, nrow = n, ncol = dim(data[[py]]$sZ_ts)[[1]])
          tem.sZ_ts[data[[py]]$index.y, ] <- t(data[[py]]$sZ_ts[, t, ])
        }
      }
      temp.sz <- cbind(temp.sz, tem.sZ_ts)
    }

    if(!is.null(temp.x)){
      tsv.x[[t]] <- spam::as.spam(as.matrix(temp.x))
    }
    if(!is.null(temp.z)){
      tsv.z[[t]] <- spam::as.spam(as.matrix(temp.z))
    }
    # tsv.sx[[t]] <- spam::as.spam(as.matrix(temp.sx))
    if(!is.null(temp.sz)){
      tsv.sz[[t]] <- spam::as.spam(as.matrix(temp.sz))
    }

  }
  if(length(tsv.x) == 0){
    tsv.x <- NULL
  }
  if(length(tsv.z) == 0){
    tsv.z <- NULL
  }
  if(length(tsv.sz) == 0){
    tsv.sz <- NULL
  }
  return(list(tsv.y = tsv.y, tsv.x = tsv.x, tsv.z = tsv.z,  tsv.sz = tsv.sz))
}


IniEnKF <- function (data, Para.List, test = NULL, IDE = TRUE, Hv.Zg = 0, Au.Rg = 0){
  #----------------------------------------------------------------------------
  #-- Stack Zts and H
  train.data.trans.tsv <- tranform.list.to.matrix(data)
  train.Z.tsv          <- test.Z.tsv <- train.sZ.tsv <- test.sZ.tsv <- NULL
  if(!is.null(test)){
    test.data.trans.tsv <- tranform.list.to.matrix(test)
  }
  if(!is.null(train.data.trans.tsv$tsv.z)){
    train.Z.tsv <- train.data.trans.tsv$tsv.z
    if(!is.null(test)){
      test.Z.tsv  <- test.data.trans.tsv$tsv.z
    }
  }
  if(!is.null(train.data.trans.tsv$tsv.sz)){
    train.sZ.tsv <- train.data.trans.tsv$tsv.sz
    if(!is.null(test)){
      test.sZ.tsv  <- test.data.trans.tsv$tsv.sz
    }
  }

  train.H.basis <- list()
  train.H.basis[[1]] <- data[[1]]$Hs
  if(length(data) > 1){
    for(py in 2:length(data)){
      train.H.basis[[py]] <- data[[py]]$Hs
    }
  }

  test.H.basis <- list()
  if(!is.null(test)){
    test.H.basis[[1]] <- test[[1]]$Hs
    if(length(data) > 1){
      for(py in 2:length(data)){
        test.H.basis[[py]] <- test[[py]]$Hs
      }
    }
  }

  # -- transform process
  Mv <- list()
  # term <- c("Intercept", names(Para.List[[1]]$st.RF))
  if(!is.null(data[[py]]$G_ts)){
    term <- c(dimnames(data[[py]]$G_ts)[[1]])
    for(Pg in 1:(length(term))){
      Mv[[term[Pg]]] <- Para.List[[1]]$st.RF[[term[Pg]]]$rho.v$mu.rho.v
      if(IDE){
        Mv[[term[Pg]]] <- list()
        for(g in 1:data$Grid.infor$summary$res){
          Mv[[term[Pg]]][[g]] <-  as(Para.List$st.RF[[term[Pg]]]$rho.v$mu.rho.v[g] *
                                       fields::Wendland(data$Grid.infor$level[[g]]$BAUs.Dist
                                                        , theta = Para.List$st.RF[[term[Pg]]]$Phi.v$mu.Phi.v[g]
                                                        , dimension = 1, k = 1), "dgCMatrix")


          lamda <- rARPACK::eigs_sym(Mv[[term[Pg]]][[g]], 1)$values
          if(lamda >= 1){
            if(is.complex(lamda[1]))
            {
              Lamda <- max(Mod(as.vector(lamda)))
              cat("\n .............................. \n")
              cat("Max-Eigen of the Mv:", Lamda)
              # print(t2 - t1)
              cat("\n .............................. \n")
              if(Lamda > 1) Mv[[term[Pg]]][[g]] <- Mv[[term[Pg]]][[g]]/(abs(Lamda) + 1E-3)
            }else{
              Lamda <- max(lamda)
              cat("\n .............................. \n")
              cat("Max-Eigen of the Mv:", Lamda)
              # print(t2 - t1)
              cat("\n .............................. \n")
              if(Lamda > 1) Mv[[term[Pg]]][[g]] <- Mv[[term[Pg]]][[g]]/(abs(Lamda) + 1E-3)
            }
          }
        }
      }
    }
  }

  for(py in 1:length(names(data))){
    if(!is.null(data[[py]]$sG_ts)){
      term <- dimnames(data[[py]]$sG_ts)[[1]]
      for(sPg in 1:(length(term))){
        Mv[[paste0(names(data)[py], ".", term[sPg])]] <- Para.List[[py]]$st.sRF[[paste0(names(data)[py], ".", term[sPg])]]$rho.v$mu.rho.v
      }
    }
  }


  yts.minus.xts <- train.data.trans.tsv$tsv.y - Hv.Zg - Au.Rg
  if(!is.null(data[[py]]$X_ts)){
    yts.minus.xts <- train.data.trans.tsv$tsv.y - Xts.beta.fun(
      tsv.x = train.data.trans.tsv$tsv.x,
      alpha = Para.List[[1]]$beta$mu.beta,
      self  = FALSE)

    for(py in 1:length(data)){
      if(!is.null(data[[py]]$sX_ts)){
        yts.minus.xts[, data[[py]]$index.y] <- yts.minus.xts[, data[[py]]$index.y] -
          Xts.beta.fun(tsv.x = data[[py]]$sX_ts,
                       alpha = Para.List[[py]]$sbeta$mu.beta,
                       self  = TRUE)
      }
    }
  }else{
    for(py in 1:length(data)){
      if(!is.null(data[[py]]$sX_ts)){
        yts.minus.xts[, data[[py]]$index.y] <- yts.minus.xts[, data[[py]]$index.y] -
          Xts.beta.fun(tsv.x = data[[py]]$sX_ts,
                       alpha = Para.List[[py]]$sbeta$mu.beta,
                       self  = TRUE)
      }
    }
  }


  return(list(yts.minus.xts = yts.minus.xts,
              train.Z.tsv = train.Z.tsv,
              test.Z.tsv = test.Z.tsv,
              train.sZ.tsv = train.sZ.tsv,
              test.sZ.tsv = test.sZ.tsv,
              train.H.basis = train.H.basis,
              test.H.basis = test.H.basis,
              # ini.alpha.cov = ini.alpha.cov,
              # alpha.cov = alpha.cov,
              # Mb = Mb,
              Mv = Mv
  ))
}

make.G.Mat <- function (data = NULL, test = NULL){
  data.G.Mat <- test.G.Mat <- Y.index <- Res.indic <- list()
  if(!is.null(data)){
    Nt     <- data[[py]]$Nt
    n      <- 0
    Py     <- length(names(data))
    nKnots <- vector()
    for(py in 1:Py){
      nKnots[py] <- data[[py]]$nKnots
      n <- n + data[[py]]$n
    }
    sPG  <- 0
    ID   <- NULL
    time <- dimnames(data[[py]]$Y_ts)[[1]]
    for(py in 1:Py){
      sPG <- sPG + data[[py]]$sPg
      ID  <- c(ID, paste0("Y", py, ".", dimnames(data[[py]]$Y_ts)[[2]]))
    }

    if(data[[py]]$Pg > 0){
      pb <- progress::progress_bar$new(format = "pub.G.Mat (train): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = data[[py]]$Pg
                                       , clear = FALSE)
      for(Pg in 1:data[[py]]$Pg){
        pb$tick()
        Sys.sleep(1 / (2*data[[py]]$Pg))
        Y.index[[dimnames(data[[py]]$G_ts)[[1]][Pg]]] <- 1:n
        Res.indic[[dimnames(data[[py]]$G_ts)[[1]][Pg]]] <- 0
        data.G.Mat[[dimnames(data[[py]]$G_ts)[[1]][Pg]]] <- array(0,
                                                                 dim =  c(Nt, n, nKnots),
                                                                 dimnames = list(time, ID,
                                                                                 paste0("Knots.", 1:nKnots)))
        for(py in 1:Py){
          for(t in 1:Nt){
            # G_ts <- NULL
            # for(py in 1:Py){
            # G_ts <- c(G_ts, as.vector(data[[py]]$G_ts[Pg, t, ]))
            G_ts <- c(as.vector(data[[py]]$G_ts[Pg, t, ]))
            data.G.Mat[[dimnames(data[[py]]$G_ts)[[1]][Pg]]][t, data[[py]]$index.y,] <- (G_ts*as.matrix(data[[py]]$Hs))
          }

        }
      }
    }
    if(sPG > 0){
      pb <- progress::progress_bar$new(format  = "self.G.Mat (train): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = sPG
                                       , clear = FALSE)
      for(py in 1:Py){
        if(data[[py]]$sPg > 0){
          for(sPg in 1:data[[py]]$sPg){
            Res.indic[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]] <- py
            Y.index[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]] <- data[[py]]$index.y
            pb$tick()
            Sys.sleep(1 / (2*sPG))
            data.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]] <- array(0, dim   =  c(Nt, n, nKnots[py]), dimnames = list(time, ID, paste0("Knots.", 1:nKnots[py])))
            for(t in 1:Nt){
              # G_ts                     <- rep(0, n)
              # G_ts[data[[py]]$index.y] <- as.vector(data[[py]]$sG_ts[sPg, t, ])
              G_ts                     <- as.vector(data[[py]]$sG_ts[sPg, t, ])
              data.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][sPg])]][t, data[[py]]$index.y,] <- (G_ts*as.matrix(data[[py]]$Hs))
            }
          }
        }
      }
    }
  }else{
    data.G.Mat <- Y.index <- Res.indic <- NULL
  }
  if(!is.null(test)){
    Nt     <- test[[1]]$Nt
    # nKnots <- test[[1]]$nKnots
    n.test <- 0
    Py     <- length(names(test))
    for(py in 1:Py){
      n.test <- n.test + test[[py]]$n
    }

    time.test <- ID.test <- NULL
    time.test <- c(dimnames(test[[1]]$Y_ts)[[1]])
    sPG       <- 0
    for(py in 1:Py){
      sPG <- sPG + test[[py]]$sPg
      ID.test <- c(ID.test, dimnames(test[[py]]$Y_ts)[[2]])
    }
    if(test[[1]]$Pg > 0){
      pb <- progress::progress_bar$new(format = "pub.G.Mat (test): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = test[[1]]$Pg
                                       , clear = FALSE)

      for(Pg in 1:test[[1]]$Pg){
        pb$tick()
        Sys.sleep(1 / (2*test[[1]]$Pg))
        test.G.Mat[[dimnames(test[[1]]$G_ts)[[1]][Pg]]] <- array(0,
                                                                 dim =  c(Nt, n.test, nKnots),
                                                                 dimnames = list(time.test, ID.test,
                                                                                 paste0("Knots.", 1:nKnots)))
        for(py in 1:Py){
          for(t in 1:Nt){
            #  G_ts <- NULL
            # G_ts <- c(G_ts, as.vector(test[[py]]$G_ts[Pg, t, ]))
            # for(py in 1:Py){
            G_ts <- c(as.vector(test[[py]]$G_ts[Pg, t, ]))
            test.G.Mat[[dimnames(test[[1]]$G_ts)[[1]][Pg]]][t, test[[py]]$index.y, ] <- G_ts*
              as.matrix(test[[py]]$Hs)
          }

        }
      }
    }
    if(sPG > 0){
      pb <- progress::progress_bar$new(format = "self.G.Mat (test): (:spin) |:bar| :current/:total (:percent in :elapsed)"
                                       , total = sPG
                                       , clear = FALSE)
      for(py in 1:Py){
        if(test[[py]]$sPg > 0){
          for(Pg in 1:test[[py]]$sPg){
            pb$tick()
            Sys.sleep(1 / (2*sPG))
            test.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][Pg])]] <- array(0,
                                                                                                     dim =  c(Nt, n.test, nKnots[py]),
                                                                                                     dimnames = list(time.test, ID.test,
                                                                                                                     paste0("Knots.", 1:nKnots[py])))
            for(t in 1:Nt){
              # G_ts                     <- rep(0, n.test)
              # G_ts[test[[py]]$index.y] <- as.vector(test[[py]]$sG_ts[Pg, t, ])
              G_ts <- as.vector(test[[py]]$sG_ts[Pg, t, ])
              test.G.Mat[[paste0(names(data)[py], ".", dimnames(data[[py]]$sG_ts)[[1]][Pg])]][t,test[[py]]$index.y,] <- G_ts*as.matrix(test[[py]]$Hs) #ini.IEnKF.Input$test.H.basis
            }
          }
        }

      }
    }

  }else{
    test.G.Mat <- NULL
  }
  return(list(data.G.Mat = data.G.Mat, test.G.Mat = test.G.Mat, Y.index = Y.index, Res.indic = Res.indic))
}


sfunc <- function(data, Ks, spTaper = NULL, cores = 1, py = 1){
  Pt_1_1 =  parallel::mcmapply(seq_len(data[[py]]$Nt + 1),
                               FUN = function(t){

                                # A <- spam::as.spam.dgCMatrix(spam::tcrossprod.spam(Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[1]], ])*spTaper$taper[[1]])

                                s11 <- spam::as.spam(spam::tcrossprod.spam(Ks$Vt.mean[t,
                                                                                      data[[py]]$Grid.infor$summary$g.index[[1]]])) + spam::as.spam.dgCMatrix(
                                                                                        cov(t(Ks$Vt.ens[t, data[[py]]$Grid.infor$summary$g.index[[1]], ]),
                                                                                            t(Ks$Vt.ens[t, data[[py]]$Grid.infor$summary$g.index[[1]], ]))*
                                                                                          Ks$rho.time[paste(0)] *spTaper$taper[[1]]);

                                 # s11 <- spam::as.spam(spam::tcrossprod.spam(Ks$Vt.mean[t,
                                 #                                                       data[[py]]$Grid.infor$summary$g.index[[1]]])) + spam::as.spam.dgCMatrix(
                                 #                                                         HDCM::Empical_Cov_taper(Ks$Vt.ens[t,
                                 #                                                                                           data[[py]]$Grid.infor$summary$g.index[[1]], ],
                                 #                                                                                 Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[1]], ],
                                 #                                                                                 Ks$rho.time[paste(0)] *spTaper$taper[[1]], sym = T,
                                 #                                                                                 nThreads = cores));
                                 if(data[[py]]$Grid.infor$summary$res > 1){
                                   for(g in 2:data[[py]]$Grid.infor$summary$res){
                                     s11 <- spam::bdiag.spam(s11,
                                                             spam::tcrossprod.spam(Ks$Vt.mean[t,
                                                                                              data[[py]]$Grid.infor$summary$g.index[[g]]]) + spam::as.spam.dgCMatrix(
                                                                                                HDCM::Empical_Cov_taper(Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                                                        Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                                                        Ks$rho.time[paste(0)] *spTaper$taper[[g]],
                                                                                                                        sym = T, nThreads = cores))
                                     );
                                   }
                                 };
                                 return(s11)},
                               SIMPLIFY = F, mc.cores = 1)


  St_1_0 =  parallel::mcmapply(seq_len(data[[py]]$Nt),
                               FUN = function(t){
                                 s10 <- spam::as.spam(spam::tcrossprod.spam(Ks$Vt.mean[t + 1,data[[py]]$Grid.infor$summary$g.index[[1]]],
                                                                            Ks$Vt.mean[t,data[[py]]$Grid.infor$summary$g.index[[1]]])) +
                                   spam::as.spam.dgCMatrix(
                                     HDCM::Empical_Cov_taper(Ks$Vt.ens[t + 1,data[[py]]$Grid.infor$summary$g.index[[1]], ],
                                                             Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[1]], ],
                                                             Ks$rho.time[paste(1)] *spTaper$taper[[1]],
                                                             sym = F,
                                                             nThreads = cores));
                                 if(data[[py]]$Grid.infor$summary$res > 1){
                                   for(g in 2:data[[py]]$Grid.infor$summary$res){
                                     s10 <- spam::bdiag.spam(s10,
                                                             spam::tcrossprod.spam(Ks$Vt.mean[t + 1,data[[py]]$Grid.infor$summary$g.index[[g]]],
                                                                                   Ks$Vt.mean[t,data[[py]]$Grid.infor$summary$g.index[[g]]]) +
                                                               spam::as.spam.dgCMatrix(
                                                                 HDCM::Empical_Cov_taper(Ks$Vt.ens[t + 1,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                         Ks$Vt.ens[t,data[[py]]$Grid.infor$summary$g.index[[g]], ],
                                                                                         Ks$rho.time[paste(1)] *spTaper$taper[[g]],
                                                                                         sym = F,
                                                                                         nThreads = cores))
                                     );
                                   }
                                 };
                                 return(s10)},
                               SIMPLIFY = F, mc.cores = 1)
  temp <- 0
  for(t in (2:data[[py]]$Nt)){
    temp <- temp + Pt_1_1[[t]]
  }

  S00 <- temp + Pt_1_1[[1]]
  S11 <- temp + Pt_1_1[[data[[py]]$Nt + 1]]
  S10 <- 0
  for(t in seq_len(data[[py]]$Nt)){
    S10 <- S10 + St_1_0[[t]]#temp10[,, t]
    # S01 <- S01 + temp01[[t]]#temp01[,, t]
  }
  rm(St_1_0, temp);gc()


  #beta.t
  beta.S11 <- beta.S10 <- beta.S00 <- beta.t00 <- list()
  # self.beta.S11 <- self.beta.S10 <- self.beta.S00 <- self.beta.t00 <- list()
  if(dim(Ks$Vt.ens)[[2]] > data[[py]]$nKnots){
    Pz <- 0
    if((!is.null(data[[py]]$Z_ts))){
      Pz <- dim(data[[py]]$Z_ts)[1]
      ind <-(data[[py]]$Grid.infor$summary$nKnots + 1):(data[[py]]$Pz + data[[py]]$Grid.infor$summary$nKnots)
      name <- "public.gt"
      beta.S11[[name]] <- beta.S00[[name]] <- matrix(NA, ncol = 1, nrow = data[[py]]$Pz)
      beta.S10[[name]] <- beta.t00[[name]] <- matrix(NA, ncol = 1, nrow = data[[py]]$Pz)
      rownames(beta.S11[[name]]) <-
        rownames(beta.S10[[name]]) <-
        rownames(beta.S00[[name]]) <-
        rownames(beta.t00[[name]]) <-
        colnames(Ks$Vt.mean[, ind])



      colnames(beta.S11[[name]]) <- "S11"
      colnames(beta.S10[[name]]) <- "S10"
      colnames(beta.S00[[name]]) <- "S00"
      colnames(beta.t00[[name]]) <- "beta.t00"

      for(pz in 1:data[[py]]$Pz){
        beta_1_1 <- parallel::mcmapply(seq_len(data[[py]]$Nt + 1),
                                       FUN = function(t){return(as.vector(tcrossprod(
                                         Ks$Vt.mean[t, ind[pz]]) +
                                           var(Ks$Vt.ens[t, ind[pz], ])))},
                                       SIMPLIFY = F, mc.cores = 1) %>% unlist()
        beta.S11[[name]][pz, ] <- sum(beta_1_1[-1])
        beta.S00[[name]][pz, ] <- sum(beta_1_1[-(data[[py]]$Nt + 1)])
        beta.t00[[name]][pz, ] <- sum(beta_1_1[1])
      }

      for (pz in 1:data[[py]]$Pz) {
        beta_1_0 =  parallel::mcmapply(seq_len(data[[py]]$Nt),
                                       FUN = function(t){return(as.vector(tcrossprod(
                                         Ks$Vt.mean[t + 1, ind[pz]] *
                                           Ks$Vt.mean[t, ind[pz]]) +
                                           cov(Ks$Vt.ens[t + 1, ind[pz], ],
                                               Ks$Vt.ens[t, ind[pz], ])))},
                                       SIMPLIFY = F, mc.cores = 1) %>% unlist()
        beta.S10[[name]][pz, ] <- sum(beta_1_0)
      }

    }

    for(py in 1:length(names(data))){
      if((!is.null(data[[py]]$sZ_ts))){

        if(py == 1){
          sPz.s <- data[[py]]$Grid.infor$summary$nKnots + Pz + 1
          sPz <- 0
          ind <- (sPz.s):(dim(data[[py]]$sZ_ts)[[1]] + data[[py]]$Grid.infor$summary$nKnots + Pz)

        }else{
          ind <- (max(ind) + 1):(sPz + dim(data[[py]]$sZ_ts)[[1]] + data[[py]]$Grid.infor$summary$nKnots + Pz)
        }
        sPz <- sPz + dim(data[[py]]$sZ_ts)[[1]]



        name <- paste0(names(data)[py], ".gt")
        beta.S11[[name]] <- beta.S00[[name]] <- matrix(NA, ncol = 1, nrow = dim(data[[py]]$sZ_ts)[[1]])
        beta.S10[[name]] <- beta.t00[[name]] <- matrix(NA, ncol = 1, nrow = dim(data[[py]]$sZ_ts)[[1]])
        rownames(beta.S11[[name]]) <-
          rownames(beta.S10[[name]]) <-
          rownames(beta.S00[[name]]) <-
          rownames(beta.t00[[name]]) <-
          colnames(Ks$Vt.mean[, ind])



        colnames(beta.S11[[name]]) <- "S11"
        colnames(beta.S10[[name]]) <- "S10"
        colnames(beta.S00[[name]]) <- "S00"
        colnames(beta.t00[[name]]) <- "beta.t00"

        for(pz in 1:data[[py]]$sPz){
          beta_1_1 <- parallel::mcmapply(seq_len(data[[py]]$Nt + 1),
                                         FUN = function(t){return(as.vector(tcrossprod(
                                           Ks$Vt.mean[t, ind[pz]]) +
                                             var(Ks$Vt.ens[t, ind[pz], ])))},
                                         SIMPLIFY = F, mc.cores = 1) %>% unlist()
          beta.S11[[name]][pz, ] <- sum(beta_1_1[-1])
          beta.S00[[name]][pz, ] <- sum(beta_1_1[-(data[[py]]$Nt + 1)])
          beta.t00[[name]][pz, ] <- sum(beta_1_1[1])
        }

        for (pz in 1:data[[py]]$sPz) {
          beta_1_0 =  parallel::mcmapply(seq_len(data[[py]]$Nt),
                                         FUN = function(t){return(as.vector(tcrossprod(
                                           Ks$Vt.mean[t + 1, ind[pz]] *
                                             Ks$Vt.mean[t, ind[pz]]) +
                                             cov(Ks$Vt.ens[t + 1, ind[pz], ],
                                                 Ks$Vt.ens[t, ind[pz], ])))},
                                         SIMPLIFY = F, mc.cores = 1) %>% unlist()
          beta.S10[[name]][pz, ] <- sum(beta_1_0)
        }

      }
    }
  }

  return(S = list(S00 = S00, S11 = S11, S10 = S10,
                  V.t00 = Pt_1_1[[1]],
                  beta.S11 = beta.S11, beta.S10 = beta.S10,
                  beta.S00 = beta.S00,
                  beta.t00 = beta.t00#, St.10 = St_1_0
  ))
}



non_f_var_alpha <- function(X){
  Adj.Mat <- exp(-(G/mu.Phi.v[1])) 
  mat.x <- as.matrix(var.covariable)

  mu.var <- exp(mat.x %*% c(X))
  Adj.Mat <- mu.tau.sq*solve(diag(as.vector(mu.var^(0.5))) %*% Adj.Mat %*% diag(as.vector(mu.var^(0.5))))

  log.det.D <- Matrix::determinant(Adj.Mat)$modulus/2 

  f     <- 
    (Nt + 1) *log.det.D - (sum(spam::diag.spam(Adj.Mat %*% S11))  +
                                      sum(spam::diag.spam(Adj.Mat %*% V.t00)) -
                                      2*Rho*sum(spam::diag.spam(Adj.Mat %*% S10)) +
                                      Rho.sq*sum(spam::diag.spam(Adj.Mat %*% S00)))/2 -
    0.5*matrix(X - LR.coef.prior$mu.delta, nrow = 1) %*% solve(LR.coef.prior$sigma.sq) %*% matrix(X - LR.coef.prior$mu.delta, ncol = 1)
  rm(Adj.Mat);
  return(-f)
}


non_f_range <- function(X){
  Hs.fun <- function(X, residuals, Hdist, sigma.sq, Vt.mean, index.y){
    phi <- exp(X[1])
    Hs.temp <- exp(-Hdist/phi)
    Hs <- as(Hs.temp/apply(Hs.temp, 1, sum), "sparseMatrix")
    Yt <-  as.vector(residuals[, index.y] - t(Hs %*% t(Vt.mean)))
    return(-Yt%*%Yt/(2*sigma.sq))
  }

  phi <- exp(X[1])#/exp(X[2])
  Adj.Mat <- exp(-G/phi)

  mat.x <- as.matrix(var.covariable)

  mu.var  <-  exp(mat.x %*% mu.alpha)
  Adj.Mat <- mu.tau.sq*solve(diag(as.vector(mu.var^(0.5))) %*% Adj.Mat %*% diag(as.vector(mu.var^(0.5))))
  log.det.D <- Matrix::determinant(Adj.Mat)$modulus/2

  f     <- Hs.fun(X, residuals, Hdist, sigma.sq, Vt.mean, index.y) -
    (Nt + 1) *log.det.D - (sum(spam::diag.spam(Adj.Mat %*% S11))  +
                                      sum(spam::diag.spam(Adj.Mat %*% V.t00)) -
                                      2*Rho*sum(spam::diag.spam(Adj.Mat %*% S10)) +
                                      Rho.sq*sum(spam::diag.spam(Adj.Mat %*% S00)))/2 +
    (prior.Phi.v.a - 1)*log(phi) - X[1]*prior.Phi.v.b

  rm(Adj.Mat);
  return(-f)
}


plot.lnormal <- function(mu, var, xlab = "x", high = 1, length.out = 1e2){
  x11 <- qlnorm(1e-8, mu, sqrt(var))
  x12 <- qlnorm(1 - 1e-8, mu, sqrt(var))
  x1 <- seq(x11, x12, by = (x12 - x11)/(length.out - 1))

  f1 <- function(x) dlnorm(x, mu, sqrt(var))
  plot(x1, f1(x1), cex = 2, lwd = 2,
       ylab = paste0("f(", xlab, ")"), type = "l",
       xlab = xlab, main = paste0("Density of ", xlab))
}





###################################################################
spT_validation <- function (z = NULL, zhat = NULL, sigma = NULL,
                            zhat.Ens = NULL, names = FALSE, CC = T){
  VMSE <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sum(u^2)/length(u), 4)
  }
  ## root mean square error
  RMSE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(z - zhat)
    u <- x[!is.na(x)]
    round(sqrt(sum(u^2)/length(u)), 4)
  }
  MAE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  MAPE <<- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- abs(c(zhat - z))/abs(z)
    u <- x[!is.na(x)]
    u <- u[!is.infinite(u)]
    round(sum(u)/length(u) * 100, 4)
  }
  ## normalised mean gross error
  NMGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    NMGE <- abs(c(x[[2]] - x[[1]]))/sum(x[[1]])
    round(mean(NMGE), 4)
  }

  BIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/length(u), 4)
  }
  rBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    round(sum(u)/(length(u) * mean(z, na.rm = TRUE)), 4)
  }
  ## normalised mean bias
  nBIAS <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    round(sum(x, na.rm = T)*100/sum(z, na.rm = T), 4)
  }

  rMSEP <- function(z, zhat) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    x <- c(zhat - z)
    u <- x[!is.na(x)]
    y <- c(mean(zhat, na.rm = TRUE) - z)
    v <- y[!is.na(y)]
    round(sum(u^2)/sum(v^2), 4)
  }
  ## correlation coefficient
  Coef <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    Coef <- suppressWarnings(cor(x[[2]], x[[1]])) ## when SD=0; will return NA
    round(Coef, 4)
  }

  ##  Coefficient of Efficiency
  COE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    COE <- 1 - sum(abs(x[[2]] - x[[1]])) / sum(abs(x[[1]] - mean(x[[1]])))
    round(COE, 4)
  }
  ## fraction within a factor of two
  FAC2 <<- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    ratio <- x[[2]]/x[[1]]
    len <- length(ratio)
    if (len > 0) {
      FAC2 <- length(which(ratio >= 0.5 & ratio <= 2)) / len
    } else {
      FAC2 <- NA
    }
    round(FAC2, 4)
  }
  ## mean gross error
  MGE <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame() %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)
    MGE <- round(mean(abs(x[[2]] - x[[1]])), 4)
  }

  ##  Index of Agreement
  IOA <- function(z, zhat) {
    z <- z %>% as.data.frame() %>% as.matrix() %>% as.vector()
    zhat <- zhat %>% as.data.frame()  %>% as.matrix() %>% as.vector()
    da <- data.frame(z = z, zhat = zhat)
    x <- na.omit(da)

    LHS <- sum(abs(x[[2]] - x[[1]]))
    RHS <- 2 * sum(abs(x[[1]] - mean(x[[1]])))

    if (LHS <= RHS) IOA <- 1 - LHS / RHS else IOA <- RHS / LHS - 1
    round(IOA, 4)
  }


  CRPS <<- function(z, zhat, sigma = NA){
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    if(is.na(sigma[1])){
      da <- data.frame(z = as.vector(z), zhat = as.vector(zhat))
      x <- na.omit(da)
      x$sigma <- sd(x$zhat)
    }else{
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat),
                       sigma = as.vector(sigma))
      x <- na.omit(da)
    }
    CRPS <- round(mean(verification::crps(x$z, cbind(x$zhat, x$sigma))$CRPS), 4)
  }

  # CRPS <<- function(x,mu,sigma)
  # {
  #   x.z = (x-mu)/sigma;
  #   sigma * (x.z*(2*pnorm(x.z)-1) + 2*dnorm(x.z) - 1/sqrt(pi))
  # }

  Interval.score <<- function(z, zhat, sigma = NA, alpha=0.05) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    if(is.na(sigma[1])){
      da <- data.frame(z = as.vector(z), zhat = as.vector(zhat))
      x <- na.omit(da)
      x$sigma <- sd(x$zhat)
    }else{
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat),
                       sigma = as.vector(sigma))
      x <- na.omit(da)
    }

    ## Compute the Interval score:
    hw <- -qnorm(alpha/2) * x$sigma
    scores <- 2 * hw + (2/alpha) * (((x$z - hw) - x$zhat) * (x$zhat < x$z - hw) +
                                      (x$zhat - (x$z + hw)) * (x$zhat > x$z + hw))
    # attr(scores, "score") <- mean(scores, na.rm=TRUE)
    scores <- mean(scores, na.rm=TRUE)
  }

  CRPS.ES <- function(z, zhat, zhat.Ens) {
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))

    Nt <- nrow(z)
    n <- ncol(z)
    # cat("Nt = ", Nt)
    if(!is.null(zhat.Ens))
    {
      # z <- as.matrix(z)
      test.index <- which(is.na(z), arr.ind = T)
      if(nrow(test.index) >0){ z[test.index] <- 0 }
      CRPS <- ES <- matrix(NA, nrow = nrow(z), ncol = n)
      for(t in 1:Nt)
      {
        CRPS[t, ] <- crps_sample(z[t, ], zhat.Ens[t,,])
        ES[t, ] <- es_sample(z[t, ], zhat.Ens[t,,])
      }
      CRPS[test.index] <-  ES[test.index] <- NA
      #crps
      CRPS <- round(mean(CRPS, na.rm = T), 4)
      ES <- round(mean(ES, na.rm = T), 4)
    }else{
      # z <- as.matrix(z)
      # zhat <- as.matrix(zhat)

      CRPS <- ES <- vector()
      for(t in 1:Nt)
      {
        test.index <- which(is.na(z[t, ]))%>% as.vector()
        if(length(test.index) == 0 ) {
          z0 <- z[t, ] %>% as.vector()
          zt <- zhat[t, ] %>% as.matrix()
        }else{
          z0 <- z[t, -test.index] %>% as.vector()
          zt <- zhat[t, -test.index] %>% as.matrix()
        }
        CRPS[t] <- mean(crps_sample(z0, zt))
        ES[t] <-  es_sample(z0, zt)
      }
      CRPS <- round(mean(CRPS, na.rm = T), 4)
      ES <- round(mean(ES, na.rm = T), 4)
    }
    return(list(CRPS = CRPS, ES = ES))
  }

  CRPS.samples <<- function(true, samples, weights=-1) {
    true <- as.vector(true)
    Ne <- dim(samples)[3]
    if(!is.null(Ne)){
      Sample <- samples
      samples <- NULL
      for(e in 1:Ne){
        samples <- cbind(samples, as.vector(Sample[,, e]))
      }
    }


    m <- length(true)

    if((is.vector(samples) & m > 1) || ncol(samples)==1) {

      expect1 = abs(samples-true)
      expect2 = 0

    }else{

      if(is.vector(samples)){ samples=matrix(samples,nrow=1,ncol=length(samples)) }

      n=ncol(samples)

      if(sum(weights==-1)>0) {
        weights=matrix(1/n,nrow=nrow(samples),ncol=n)
      }

      expect1 = rowSums(abs(samples-true)*weights)

      sortind=t(apply(samples,1,order))

      for(i in 1:m) {
        samples[i,]=samples[i,sortind[i,]]
        weights[i,]=weights[i,sortind[i,]]
      }

      cum.weights=t(apply(weights,1,cumsum))[,1:(n-1)]

      if(m==1) {
        expect2 = sum(cum.weights*(1-cum.weights)*(samples[2:n]-samples[1:(n-1)]))
      } else {
        expect2 = rowSums(cum.weights*(1-cum.weights)*(samples[,2:n]-samples[,1:(n-1)]))
      }

    }

    return(round(mean(expect1 - expect2), 4))

  }
  Ignorance_score <<- function(z, zhat, zhat.Ens){
    z <- as.matrix(as.data.frame(z))
    zhat <- as.matrix(as.data.frame(zhat))
    if(!is.null(zhat.Ens)){
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat),
                       zhat.Ens = as.vector(apply(zhat.Ens,
                                                  c(1, 2), sd)))
      x <- na.omit(da)
      IGN <- verification::crps(x$z,
                                cbind(x$zhat,
                                      x$zhat.Ens))$IGN
    }else{
      da <- data.frame(z = as.vector(z),
                       zhat = as.vector(zhat))
      x <- na.omit(da)
      IGN <- verification::crps(x$z,
                                cbind(x$zhat, sd(x$zhat)))$IGN
    }
    return(IGN)
  }
  if (names == TRUE) {
    cat("##\n Mean Squared Error (MSE) \n Root Mean Squared Error (RMSE) \n Mean Absolute Error (MAE) \n Mean Absolute Percentage Error (MAPE) \n Normalized Mean Error(NME: %) \n Bias (BIAS) \n Relative Bias (rBIAS) \n Normalized Mean Bias(NMB: %) \n Relative Mean Separation (rMSEP) \n Correlation coefficient (Coef.) \n The fraction of predictions within a factor of two of observations (FAC2) \n Continuous ranked probability score (CRPS) based on sample mean and standard deviation \n CRPS based on empirical distribution function \n Energy score (ES) based on empirical distribution function \n##\n")
  }
  if((!is.null(z))&(!is.null(zhat))){

    out <- NULL
    if(CC){
      out$MSE <- VMSE((z), (zhat))
    }

    out$RMSE <- RMSE((z), (zhat))
    if(CC){
      out$MAPE <- MAPE((z), (zhat))
      out$BIAS <- BIAS((z), (zhat))
      out$rBIAS <- rBIAS((z), (zhat))
    }

    out$Coef <- Coef((z), (zhat))
    out$FAC2 <- FAC2(z, zhat)
    out$CRPS <- CRPS(z, zhat, sigma = sigma)
    out$MAE <- MAE((z), (zhat))
    # out$INS <- Interval.score(z, zhat, sigma = sigma)
    if(!is.null(zhat.Ens)){
      out$CRPS.samples <- CRPS.samples(z, zhat.Ens)
    }else{
      out$CRPS.samples <- CRPS.samples(as.vector(z), as.vector(zhat))
    }

    # out$IGN <- Ignorance_score(z, zhat, zhat.Ens)

    if(CC){
      if(!is.null(zhat.Ens))
      {
        R <- CRPS.ES(z, zhat, zhat.Ens)
      }else{
        R <- CRPS.ES(z, zhat, zhat.Ens = NULL)
      }

      out$CRPS.sample <- R$CRPS
      out$ES.sample <- R$ES
      out$NMGE <- NMGE((z), (zhat))
      out$MGE <- MGE(z, zhat)
      out$IOA <- IOA(z, zhat)
      out$NMB <- nBIAS((z), (zhat))
      out$rMSEP <- rMSEP((z), (zhat))
      out$COE <- COE(z, zhat)
    }
    unlist(out)
  }else{
    return(0)
  }
}

Round <- function(x, n = 2)
{
  return(format(round(x, n), nsmall = n))
}

