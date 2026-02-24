###########################################################################
#                                   Create Y_ts/X_ts
###########################################################################
Construct.Fixed.effect.Data <- function(data = NULL,
                                include = list(YEAR = c(2015),
                                               month_day = c("01-01", "12-31")),
                                Y = "PM25",
                                X = NULL,
                                R = NULL,
                                Z = NULL,
                                G = NULL,
                                sX = NULL,
                                sR = NULL,
                                sZ = NULL,
                                sG = NULL,
                                date_time = "DATE_TIME",
                                siteid = "ID",
                                start.time = 0,
                                standard = F,
                                center = F,
                                start.index = 1,
                                initial.miss = NULL,
                                scaled.variable.x = NULL,
                                scaled.variable.r = NULL,
                                scaled.variable.z = NULL,
                                scaled.variable.g = NULL,
                                scaled.variable.sx = NULL,
                                scaled.variable.sr = NULL,
                                scaled.variable.sz = NULL,
                                scaled.variable.sg = NULL)
{
  if (is.null(data)) { stop("Must provide data.\n")}
  setDT(data)

  if(!is.null(include)){
    data_base <- data %>%
      dplyr::filter(
        # YEAR %in% year,
        between((as.Date(DATE_TIME)),
                (as.Date(paste0(include$YEAR[1], "-", include$month_day[1]))),
                (as.Date(paste0(include$YEAR[2], "-", include$month_day[2]))))#,
        # between(MONTH, 10, 12),
      ) %>% setorderv(cols = c(siteid, date_time)) %>% setDF()
  }else{
    data_base <- data %>% setDF()
  }


  # data_base_true <- data_base
  # setDT(data_base)
  Nt <- length(unique(data_base[, date_time]))
  SiteId <- unique(data_base[, siteid]) %>% sort()
  n <- length(SiteId)

  ######################################################################
  #                         create  data structure
  ######################################################################

  Y_ts <- matrix(NA, nrow = Nt, ncol = n)
  rownames(Y_ts) <- as.character(unique(data_base[, date_time]))
  colnames(Y_ts) <- as.character(SiteId)

  ini.miss <- NULL
  if(!is.null(initial.miss)){
    ini.miss <- matrix(NA, nrow = Nt, ncol = n)
    rownames(ini.miss) <- as.character(unique(data_base[, date_time]))
    colnames(ini.miss) <-as.character(SiteId)
  }

  X_ts <- R_ts <- Z_ts <- sZ_ts <- G_ts <- sX_ts <- sR_ts <- sG_ts <- NULL

  Px <- length(X)
  if(!is.null(X)){
    X_ts <- array(1, dim = c(Px, Nt, n)
                  , dimnames = list(as.character(c(X)),
                                    as.character(unique(data_base[, date_time])),
                                    as.character(SiteId)))
  }


  Pr <- 0
  if(!is.null(R)){
    Pr <- length(R)
    R_ts <- array(1, dim = c(Pr, Nt, n)
                  , dimnames = list(as.character(c(R)),
                                    as.character(unique(data_base[, date_time])),
                                    as.character(SiteId)))
  }


  Pz <- 0
  if(!is.null(Z)){
    Pz <- length(Z)
    Z_ts <- array(1, dim = c(Pz, Nt, n)
                  , dimnames = list(as.character(c(Z)),
                                    as.character(unique(data_base[, date_time])),
                                    as.character(SiteId)))
  }


  Pg <- 0
  # G_ts <- array(1, dim = c(Pg, Nt, n)
  #               , dimnames = list(as.character(c("Intercept", G)),
  #                                 as.character(unique(data_base[, date_time])),
  #                                 as.character(SiteId)))

  if(!is.null(G)){
    Pg <- length(G)
    G_ts <- array(1, dim = c(Pg, Nt, n)
                  , dimnames = list(as.character(c(G)),
                                    as.character(unique(data_base[, date_time])),
                                    as.character(SiteId)))
  }

  # Single covriates
  # sX_ts <- NULL
  # if(Py > 1){
    sPx <- 0
    if(!is.null(sX)){
      sPx <- length(sX)
      sX_ts <- array(1, dim = c(sPx, Nt, n)
                    , dimnames = list(sX,
                                      as.character(unique(data_base[, date_time])),
                                      as.character(SiteId)))
    }

    sRr <- 0
    if(!is.null(sR)){
      sPr <- length(sR)
      sR_ts <- array(1, dim = c(sPr, Nt, n)
                     , dimnames = list(sR,
                                       as.character(unique(data_base[, date_time])),
                                       as.character(SiteId)))
    }

    sPz <- 0
    if(!is.null(sZ)){
      sPz <- length(sZ)
      sZ_ts <- array(1, dim = c(sPz, Nt, n)
                     , dimnames = list(c(sZ),
                                       as.character(unique(data_base[, date_time])),
                                       as.character(SiteId)))
    }

   sPg <- 0
    if(!is.null(sG)){
      sPg <- length(sG)
      sG_ts <- array(1, dim = c(sPg, Nt, n)
                     , dimnames = list(sG,
                                       as.character(unique(data_base[, date_time])),
                                       as.character(SiteId)))
    }
  # }



  for(t in 1:Nt){
    Y_ts[t, ]  <- data_base[data_base$time.index == start.time + t,
                            colnames(data_base) == Y]
    if(!is.null(X)){
      for(k in 1:Px){
        X_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                   colnames(data_base) == X[k]]
      }
    }

    if(!is.null(initial.miss)){
      for(t in 1:Nt){
        ini.miss[t, ]  <- data_base[data_base$time.index == start.time + t,
                                colnames(data_base) == initial.miss]
      }
    }

    if(!is.null(R)){
      # if("Intercept" %in% R){star.ind <- 2}
        for(k in 1:Pr)
        {

          R_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                     colnames(data_base) == R[k]]
        }
    }


    if(!is.null(Z)){
      for(k in 1:Pz)
      {

        Z_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                   colnames(data_base) == Z[k]]
      }
    }


    if(!is.null(G)){
      for(k in 1:Pg)
      {
        G_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                   colnames(data_base) == G[k]]
      }
    }



    if(!is.null(sX)){
      for(k in 1:sPx)
      {
        sX_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                   colnames(data_base) == sX[k]]
      }
    }

    if(!is.null(sR)){
      for(k in 1:sPr)
      {
        sR_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                    colnames(data_base) == sR[k]]
      }
    }

    if(!is.null(sZ)){
      for(k in 1:sPz)
      {
        sZ_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                    colnames(data_base) == sZ[k]]
      }
    }

    if(!is.null(sG)){
      for(k in 1:sPg)
      {
        sG_ts[k, t, ]  <- data_base[data_base$time.index == start.time + t,
                                    colnames(data_base) == sG[k]]
      }
    }

  }


  if(isTRUE(standard)){
    if(!is.null(X)){
     if(is.null(scaled.variable.x)){
        for(k in 1:Px){
          temp <- X_ts[k,,]
         if(sum(temp[temp==1]) != (Nt*n)){
           scale.var <- scale(as.vector(X_ts[k,,]), center = center)
           if(center == TRUE){
             var.center <- as.numeric(attributes(scale.var)[2])
             var.scale <- as.numeric(attributes(scale.var)[3])
           }else{
             var.center <- 0
             var.scale <- as.numeric(attributes(scale.var)[2])
           }
           X_ts[k,,] <- matrix(scale.var[, 1],
                               nrow = Nt, ncol = ncol(X_ts[k,,]))
           scaled.variable.x <- rbind(scaled.variable.x, data.frame(Varible = X[k],
                                                                    center = var.center,
                                                                    scale = var.scale))
         }
        }
    }else{
      for(k in 1:Px){
        temp <- X_ts[k,,]
        if(sum(temp[temp==1]) != (Nt*n)){
          index <- which(scaled.variable.x$Varible == X[k])
          X_ts[k,,] <- (X_ts[k,,] - scaled.variable.x[index, 2])/scaled.variable.x[index, 3]
        }

      }
    }
}


    if(!is.null(R)){
      if(is.null(scaled.variable.r)){
        for(k in 1:Pr)#dim(X_ts)[1]
        {
          temp <- R_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(R_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            R_ts[k,,] <- matrix(scale.var[, 1],
                                nrow = Nt, ncol = ncol(R_ts[k,,]))
            scaled.variable.r <- rbind(scaled.variable.r, data.frame(Varible = R[k],
                                                                     center = var.center,
                                                                     scale = var.scale))
          }

        }
      }else{ for(k in 1:Pr){
        temp <- R_ts[k,,]
        if(sum(temp[temp==1]) != (Nt*n)){
          index <- which(scaled.variable.r$Varible == R[k])
          R_ts[k,,] <- (R_ts[k,,] - scaled.variable.r[index, 2])/scaled.variable.r[index, 3]
        }
      }
      }
    }


    if(!is.null(Z)){
      if(is.null(scaled.variable.z)){
        for(k in 1:Pz){
          temp <- Z_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(Z_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            Z_ts[k,,] <- matrix(scale.var[, 1],
                                nrow = Nt, ncol = ncol(Z_ts[k,,]))
            scaled.variable.z <- rbind(scaled.variable.z, data.frame(Varible = Z[k],
                                                                     center = var.center,
                                                                     scale = var.scale))
          }

        }
      }else{
        for(k in 1:Pz){
          temp <- Z_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            index <- which(scaled.variable.z$Varible == Z[k])
            Z_ts[k,,] <- (Z_ts[k,,] - scaled.variable.z[index, 2])/scaled.variable.z[index, 3]
          }
        }
      }
    }
    #G
    if(!is.null(G)){
      if(is.null(scaled.variable.g)){
        for(k in 1:Pg){
          temp <- G_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(G_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            G_ts[k,,] <- matrix(scale.var[, 1], nrow = Nt, ncol = ncol(G_ts[k,,]))
            scaled.variable.g <- rbind(scaled.variable.g, data.frame(Varible = G[k],
                                                                     center = var.center,
                                                                     scale = var.scale))
          }
        }
      }else{
        for(k in 1:Pg){
          temp <- G_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            index <- which(scaled.variable.g$Varible == G[k])
            G_ts[k,,] <- (G_ts[k,,] - scaled.variable.g[index, 2])/scaled.variable.g[index, 3]
          }
      }
      }
    }




#  single variates

    if(!is.null(sX)){
      if(is.null(scaled.variable.sx)){
        for(k in 1:sPx){
          temp <- sX_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(sX_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            sX_ts[k,,] <- matrix(scale.var[, 1],
                                 nrow = Nt, ncol = ncol(sX_ts[k,,]))
            scaled.variable.sx <- rbind(scaled.variable.sx, data.frame(Varible = sX[k],
                                                                       center = var.center,
                                                                       scale = var.scale))
          }
        }
      }else{
        for(k in 1:sPx){
        temp <- sX_ts[k,,]
        if(sum(temp[temp==1]) != (Nt*n)){
          index <- which(scaled.variable.sx$Varible == sX[k])
          sX_ts[k,,] <- (sX_ts[k,,] - scaled.variable.sx[index, 2])/scaled.variable.sx[index, 3]
           }
         }
      }
    }


    if(!is.null(sR)){
      if(is.null(scaled.variable.sr)){
        for(k in 1:sPr){
          temp <- sR_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(sR_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            sR_ts[k,,] <- matrix(scale.var[, 1],
                                 nrow = Nt, ncol = ncol(sR_ts[k,,]))
            scaled.variable.sr <- rbind(scaled.variable.sr, data.frame(Varible = sR[k],
                                                                       center = var.center,
                                                                       scale = var.scale))
          }
        }
      }else{
        for(k in 1:sPr){
          temp <- sR_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            index <- which(scaled.variable.sr$Varible == sR[k])
            sR_ts[k,,] <- (sR_ts[k,,] - scaled.variable.sr[index, 2])/scaled.variable.sr[index, 3]
          }
      }
      }
    }



    if(!is.null(sZ)){
      if(is.null(scaled.variable.sz)){
        for(k in 1:sPz) {
          temp <- sZ_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(sZ_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            sZ_ts[k,,] <- matrix(scale.var[, 1],
                                 nrow = Nt, ncol = ncol(sZ_ts[k,,]))
            scaled.variable.sz <- rbind(scaled.variable.sz, data.frame(Varible = sZ[k],
                                                                       center = var.center,
                                                                       scale = var.scale))
          }
        }
      }else{
        for(k in 1:sPz){
          temp <- sZ_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            index <- which(scaled.variable.sz$Varible == sZ[k])
            sZ_ts[k,,] <- (sZ_ts[k,,] - scaled.variable.sz[index, 2])/scaled.variable.sz[index, 3]
          }
      }
      }
    }


    if(!is.null(sG)){
      if(is.null(scaled.variable.sg)){
        for(k in 1:sPg) {
          temp <- sG_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            scale.var <- scale(as.vector(sG_ts[k,,]), center = center)
            if(center == TRUE){
              var.center <- as.numeric(attributes(scale.var)[2])
              var.scale <- as.numeric(attributes(scale.var)[3])
            }else{
              var.center <- 0
              var.scale <- as.numeric(attributes(scale.var)[2])
            }
            sG_ts[k,,] <- matrix(scale.var[, 1],
                                 nrow = Nt, ncol = ncol(sG_ts[k,,]))
            scaled.variable.sg <- rbind(scaled.variable.sg, data.frame(Varible = sG[k],
                                                                       center = var.center,
                                                                       scale = var.scale))
          }
        }
      }else{
        for(k in 1:sPg){
          temp <- sG_ts[k,,]
          if(sum(temp[temp==1]) != (Nt*n)){
            index <- which(scaled.variable.sg$Varible == sG[k])
            sG_ts[k,,] <- (sG_ts[k,,] - scaled.variable.sg[index, 2])/scaled.variable.sg[index, 3]
          }
        }
      }
    }
  }
  Y_X_Z_G_ts <- list(Y_ts = Y_ts,
                     ini.Miss_ts = ini.miss,
                     siteid = siteid,
                   X_ts = X_ts, R_ts = R_ts,
                   Z_ts = Z_ts, G_ts = G_ts,
                   scaled.variable.x = scaled.variable.x,
                   scaled.variable.r = scaled.variable.r,
                   scaled.variable.z = scaled.variable.z,
                   scaled.variable.g = scaled.variable.g,
                   sX_ts = sX_ts, sR_ts = sR_ts,
                   sZ_ts = sZ_ts, sG_ts = sG_ts,
                   scaled.variable.sx = scaled.variable.sx,
                   scaled.variable.sr = scaled.variable.sr,
                   scaled.variable.sz = scaled.variable.sz,
                   scaled.variable.sg = scaled.variable.sg)
  return(Y_X_Z_G_ts)
}
