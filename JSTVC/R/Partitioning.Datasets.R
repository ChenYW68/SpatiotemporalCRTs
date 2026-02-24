######################################################################
######################################################################
Partitioning.Dataset <- function(G.basic.data,
                                     Fixed.effect.Data,
                                     Object = "Object", num = 1,
                                     siteid = "ID"){
  # Hs = G.basic.data$Hs;


  Py <- length(Fixed.effect.Data)

  if(is.null(names(Fixed.effect.Data))){
    names(Fixed.effect.Data) <- c(paste0("Response.", 1:length(Fixed.effect.Data)))
  }
  train.r <- test.r <- list()
  for(py in 1:Py){
    Pred.coords <- G.basic.data[[py]]$Pred.coords

    Y_ts <- Fixed.effect.Data[[py]]$Y_ts;

    if(!is.null(Fixed.effect.Data[[py]]$ini.Miss_ts)){
      ini.Miss_ts <- Fixed.effect.Data[[py]]$ini.Miss_ts
    }

    if(!is.null(Fixed.effect.Data[[py]]$X_ts)){
      X_ts <- Fixed.effect.Data[[py]]$X_ts
    }
    if(!is.null(Fixed.effect.Data[[py]]$R_ts)){
      R_ts <- Fixed.effect.Data[[py]]$R_ts
    }

    if(!is.null(Fixed.effect.Data[[py]]$Z_ts)){
      Z_ts <- Fixed.effect.Data[[py]]$Z_ts
    }
    if(!is.null(Fixed.effect.Data[[py]]$G_ts)){
      G_ts <- Fixed.effect.Data[[py]]$G_ts
    }

    if(!is.null(Fixed.effect.Data[[py]]$sX_ts)){
      sX_ts <- Fixed.effect.Data[[py]]$sX_ts
    }
    if(!is.null(Fixed.effect.Data[[py]]$sR_ts)){
      sR_ts <- Fixed.effect.Data[[py]]$sR_ts
    }
    if(!is.null(Fixed.effect.Data[[py]]$sZ_ts)){
      sZ_ts <- Fixed.effect.Data[[py]]$sZ_ts
    }
    if(!is.null(Fixed.effect.Data[[py]]$sG_ts)){
      sG_ts <- Fixed.effect.Data[[py]]$sG_ts
    }

    setDF(Pred.coords)
    Obj.index <- which(colnames(Pred.coords) == Object)
    Obj       <- as.character(sort(unique(Pred.coords[, Obj.index])))[num]
    cat(paste0("\n To make space-time prediction on the [", Obj, "] group of ", names(Fixed.effect.Data)[py], " via JSTVC ...\n"))
    if(py == Py){
      cat(paste0("\n"))
    }
    ######################################################################
    ######################################################################
    setDF(Pred.coords)
    index <- which(colnames(Pred.coords) == siteid)
    test.ID <- unique(Pred.coords[as.character(Pred.coords[, Object]) == Obj, index]) %>% as.character()

    # test.Hs <- train.Hs <-  test.A <- train.A <- list()


    Hs <- G.basic.data[[py]]$Grid.infor$Hs
    # cat("...\n")
    test.Hs <- spam::as.spam.dgCMatrix(Hs[rownames(Hs) %in% test.ID, ])
    train.Hs <- spam::as.spam.dgCMatrix(Hs[rownames(Hs) %nin% test.ID, ])


    Hs.train.indx <- which(rownames(Hs) %nin% test.ID)
    Hs.test.indx <- which(rownames(Hs) %in% test.ID)


    test.A <- spam::as.spam.dgCMatrix(Hs[rownames(Hs) %in% test.ID, ])
    train.A <- spam::as.spam.dgCMatrix(Hs[rownames(Hs) %nin% test.ID, ])


    test.Hdist <- G.basic.data[[py]]$Grid.infor$summary$Hdist[rownames(G.basic.data[[py]]$Grid.infor$summary$Hdist) %in%
                                                          test.ID, ]
    train.Hdist <- G.basic.data[[py]]$Grid.infor$summary$Hdist[rownames(G.basic.data[[py]]$Grid.infor$summary$Hdist) %nin%
                                                          test.ID, ]

    test.site.Dist <- G.basic.data[[py]]$site.Dist[rownames(G.basic.data[[py]]$site.Dist) %in%
                                                     test.ID, colnames(G.basic.data[[py]]$site.Dist) %in%
                                                     test.ID]
    train.site.Dist <- G.basic.data[[py]]$site.Dist[rownames(G.basic.data[[py]]$site.Dist) %nin%
                                                      test.ID, colnames(G.basic.data[[py]]$site.Dist) %nin%
                                                      test.ID]


    Po.ID <- dimnames(Y_ts)[[2]]


    if(!is.null(Fixed.effect.Data[[py]]$ini.Miss_ts)){
        test.ini.Miss_ts <- ini.Miss_ts[,  Po.ID %in% test.ID]
        train.ini.Miss_ts <- ini.Miss_ts[, Po.ID %nin% test.ID]
    }else{
      test.ini.Miss_ts <- train.ini.Miss_ts <- NULL
    }




    if(!is.null(Fixed.effect.Data[[py]]$X_ts)){
      if(dim(Fixed.effect.Data[[py]]$X_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$X_ts)
        nx <- dim(Fixed.effect.Data[[py]]$X_ts)
        test.X_ts <-  array(X_ts[, , Po.ID %in% test.ID],
                            dim = c(1, nx[2], length(test.ID)),
                            dimnames = list(name[[1]], name[[2]], test.ID))
        train.X_ts <-  array(X_ts[, , Po.ID %nin% test.ID],
                             dim = c(1, nx[2], nx[3] - length(test.ID)),
                             dimnames = list(name[[1]], name[[2]],
                                             Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.X_ts <- X_ts[, , Po.ID %in% test.ID]
        train.X_ts <- X_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.X_ts <- train.X_ts <- NULL
    }




    if(dim(Y_ts)[1] > 1){
      test.Y_ts <- Y_ts[, Po.ID %in% test.ID]
      train.Y_ts <- Y_ts[, Po.ID %nin% test.ID]
    }else{
      test.Y_ts <- matrix(Y_ts[, Po.ID %in% test.ID],
                          nrow = dim(Y_ts)[1], ncol = length(test.ID))
      train.Y_ts <- matrix(Y_ts[, Po.ID %nin% test.ID],
                          nrow = dim(Y_ts)[1], ncol = length(test.ID))

      rownames(test.Y_ts) <- rownames(train.Y_ts) <- dimnames(Y_ts)[[1]]
      colnames(test.Y_ts) <- colnames(train.Y_ts) <- dimnames(Y_ts)[[2]]
    }

    if(!is.null(Fixed.effect.Data[[py]]$R_ts)){
      if(dim(Fixed.effect.Data[[py]]$R_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$R_ts)
        nr <- dim(Fixed.effect.Data[[py]]$R_ts)
        test.R_ts <-  array(R_ts[, , Po.ID %in% test.ID],
                            dim = c(1, nr[2], length(test.ID)),
                            dimnames = list(name[[1]], name[[2]], test.ID))

        train.R_ts <-  array(R_ts[, , Po.ID %nin% test.ID],
                             dim = c(1, nr[2], length(Po.ID) - length(test.ID)),
                             dimnames = list(name[[1]], name[[2]],
                                             Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.R_ts <- R_ts[, , Po.ID %in% test.ID]
        train.R_ts <- R_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.R_ts <- train.R_ts <- NULL
    }


    if(!is.null(Fixed.effect.Data[[py]]$Z_ts)){
      if(dim(Fixed.effect.Data[[py]]$Z_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$Z_ts)
        nz <- dim(Fixed.effect.Data[[py]]$Z_ts)
        test.Z_ts <-  array(Z_ts[, , Po.ID %in% test.ID],
                            dim = c(1, nz[2], length(test.ID)),
                            dimnames = list(name[[1]], name[[2]], test.ID))

        train.Z_ts <-  array(Z_ts[, , Po.ID %nin% test.ID],
                             dim = c(1, nz[2], length(Po.ID) - length(test.ID)),
                             dimnames = list(name[[1]], name[[2]],
                                             Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.Z_ts <- Z_ts[, , Po.ID %in% test.ID]
        train.Z_ts <- Z_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.Z_ts <- train.Z_ts <- NULL
    }

    if(!is.null(Fixed.effect.Data[[py]]$G_ts)){
      if(dim(Fixed.effect.Data[[py]]$G_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$G_ts)
        ng <- dim(Fixed.effect.Data[[py]]$G_ts)
        test.G_ts <-  array(G_ts[, , Po.ID %in% test.ID],
                            dim = c(1, ng[2], length(test.ID)),
                            dimnames = list(name[[1]], name[[2]], test.ID))
        train.G_ts <-  array(G_ts[, , Po.ID %nin% test.ID],
                             dim = c(1, ng[2], length(Po.ID) - length(test.ID)),
                             dimnames = list(name[[1]], name[[2]],
                                             Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.G_ts <- G_ts[, , Po.ID %in% test.ID]
        train.G_ts <- G_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.G_ts <- train.G_ts <- NULL
    }



    if(!is.null(Fixed.effect.Data[[py]]$sX_ts)){
      if(dim(Fixed.effect.Data[[py]]$sX_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$sX_ts)
        nsX <- dim(Fixed.effect.Data[[py]]$sX_ts)
        test.sX_ts <-  array(sX_ts[, , Po.ID %in% test.ID],
                             dim = c(1, nsX[2], length(test.ID)),
                             dimnames = list(name[[1]], name[[2]], test.ID))

        train.sX_ts <-  array(sX_ts[, , Po.ID %nin% test.ID],
                              dim = c(1, nsX[2], length(Po.ID) - length(test.ID)),
                              dimnames = list(name[[1]], name[[2]],
                                              Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.sX_ts <- sX_ts[, , Po.ID %in% test.ID]
        train.sX_ts <- sX_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.sX_ts <- train.sX_ts <- NULL
    }

    if(!is.null(Fixed.effect.Data[[py]]$sR_ts)){
      if(dim(Fixed.effect.Data[[py]]$sR_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$sR_ts)
        nsZ <- dim(Fixed.effect.Data[[py]]$sR_ts)
        test.sR_ts <-  array(sR_ts[, , Po.ID %in% test.ID],
                             dim = c(1, nsZ[2], length(test.ID)),
                             dimnames = list(name[[1]], name[[2]], test.ID))

        train.sR_ts <-  array(sR_ts[, , Po.ID %nin% test.ID],
                              dim = c(1, nsZ[2], length(Po.ID) - length(test.ID)),
                              dimnames = list(name[[1]], name[[2]],
                                              Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.sR_ts <- sR_ts[, , Po.ID %in% test.ID]
        train.sR_ts <- sR_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.sR_ts <- train.sR_ts <- NULL
    }



    if(!is.null(Fixed.effect.Data[[py]]$sZ_ts)){
      if(dim(Fixed.effect.Data[[py]]$sZ_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$sZ_ts)
        nsZ <- dim(Fixed.effect.Data[[py]]$sZ_ts)
        test.sZ_ts <-  array(sZ_ts[, , Po.ID %in% test.ID],
                             dim = c(1, nsZ[2], length(test.ID)),
                             dimnames = list(name[[1]], name[[2]], test.ID))

        train.sZ_ts <-  array(sZ_ts[, , Po.ID %nin% test.ID],
                              dim = c(1, nsZ[2], length(Po.ID) - length(test.ID)),
                              dimnames = list(name[[1]], name[[2]],
                                              Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.sZ_ts <- sZ_ts[, , Po.ID %in% test.ID]
        train.sZ_ts <- sZ_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.sZ_ts <- train.sZ_ts <- NULL
    }


    if(!is.null(Fixed.effect.Data[[py]]$sG_ts)){
      if(dim(Fixed.effect.Data[[py]]$sG_ts)[[1]] == 1){
        name <- dimnames(Fixed.effect.Data[[py]]$sG_ts)
        nsG <- dim(Fixed.effect.Data[[py]]$sG_ts)
        test.sG_ts <-  array(sG_ts[, , Po.ID %in% test.ID],
                             dim = c(1, nsG[2], length(test.ID)),
                             dimnames = list(name[[1]], name[[2]], test.ID))

        train.sG_ts <-  array(sG_ts[, , Po.ID %nin% test.ID],
                              dim = c(1, nsG[2], length(Po.ID) - length(test.ID)),
                              dimnames = list(name[[1]], name[[2]],
                                              Po.ID[Po.ID %nin% test.ID]))
      }else{
        test.sG_ts <- sG_ts[, , Po.ID %in% test.ID]
        train.sG_ts <- sG_ts[, , Po.ID %nin% test.ID]
      }
    }else{
      test.sG_ts <- train.sG_ts <- NULL
    }





    if(length(test.ID) == 1){
      if(!is.null(Fixed.effect.Data[[py]]$X_ts)){
      nx <- dim(train.X_ts)
      name <- dimnames(train.X_ts)
      test.X_ts <-  array(test.X_ts, dim = c(nx[1], nx[2], 1),
                          dimnames = list(name[[1]], name[[2]], test.ID))
      }

      if(!is.null(Fixed.effect.Data[[py]]$R_ts)){
        nr <- dim(train.R_ts)
        name <- dimnames(train.R_ts)
        test.R_ts <-  array(test.R_ts, dim = c(nr[1], nr[2], 1),
                            dimnames = list(name[[1]], name[[2]], test.ID))
      }


      if(!is.null(Fixed.effect.Data[[py]]$Z_ts)){
        nz <- dim(train.Z_ts)
        name <- dimnames(train.Z_ts)
        test.Z_ts <-  array(test.Z_ts, dim = c(nz[1], nz[2], 1),
                            dimnames = list(name[[1]], name[[2]], test.ID))
      }
      if(!is.null(Fixed.effect.Data[[py]]$G_ts)){
        ng <- dim(train.G_ts)
        name <- dimnames(train.G_ts)
        test.G_ts <-  array(test.G_ts, dim = c(ng[1], ng[2], 1),
                            dimnames = list(name[[1]], name[[2]], test.ID))
      }

      if(!is.null(Fixed.effect.Data[[py]]$sX_ts)){
        nsx <- dim(train.sX_ts)
        name <- dimnames(train.sX_ts)
        test.sX_ts <-  array(test.sX_ts, dim = c(nsx[1], nsx[2], 1),
                             dimnames = list(name[[1]], name[[2]], test.ID))
      }


      if(!is.null(Fixed.effect.Data[[py]]$sR_ts)){
        nsr <- dim(train.sR_ts)
        name <- dimnames(train.sR_ts)
        test.sR_ts <-  array(test.sR_ts, dim = c(nsr[1], nsr[2], 1),
                             dimnames = list(name[[1]], name[[2]], test.ID))
      }



      if(!is.null(Fixed.effect.Data[[py]]$sZ_ts)){
        nsz <- dim(train.sZ_ts)
        name <- dimnames(train.sZ_ts)
        test.sZ_ts <-  array(test.sZ_ts, dim = c(nsz[1], nsz[2], 1),
                             dimnames = list(name[[1]], name[[2]], test.ID))
      }
      if(!is.null(Fixed.effect.Data[[py]]$sG_ts)){
        nsg <- dim(train.sG_ts)
        name <- dimnames(train.sG_ts)
        test.sG_ts <-  array(test.sG_ts, dim = c(nsg[1], nsg[2], 1),
                             dimnames = list(name[[1]], name[[2]], test.ID))
      }
      test.Y_ts <- matrix(test.Y_ts, ncol = 1)
    }

    if(py == 1){
      train.index.y <-  1:dim(train.Y_ts)[2]
      test.index.y  <- 1:dim(test.Y_ts)[2]
    }else{
      train.index.y <-  (max(train.r[[py - 1]]$index.y) + 1):(dim(train.Y_ts)[2] + max(train.r[[py - 1]]$index.y))
      test.index.y  <-  (max(test.r[[py - 1]]$index.y) + 1):(dim(test.Y_ts)[2] + max(test.r[[py - 1]]$index.y))
    }

    # train.index.y <- list()
    # for(ind.py in 1:Py){
    #   train.index.y[[ind.py]] <-  ((ind.py - 1)*dim(train.Y_ts)[2] + 1):(ind.py*dim(train.Y_ts)[2])
    # }
    #
    # test.index.y <- list()
    # for(ind.py in 1:Py){
    #   test.index.y[[ind.py]] <-  ((ind.py - 1)*dim(test.Y_ts)[2] + 1):(ind.py*dim(test.Y_ts)[2])
    # }

    # which(colnames(Pred.coords) %in% Object)
    train.r[[py]] <- list(Y_ts                    = train.Y_ts,
                          ini.Miss_ts             = train.ini.Miss_ts,
                          X_ts                    = train.X_ts,
                          R_ts                    = train.R_ts,
                          Z_ts                    = train.Z_ts,
                          G_ts                    = train.G_ts,
                          sX_ts                   = train.sX_ts,
                          sZ_ts                   = train.sZ_ts,
                          sG_ts                   = train.sG_ts,
                          spCoordinates           = unique(Pred.coords[as.character(Pred.coords[, Object]) != Obj, ]),
                          Hs                      = train.Hs,
                          Hs.train.indx           = Hs.train.indx,
                          As                      = train.A,
                          Grid.infor              = G.basic.data[[py]]$Grid.infor,
                          Hdist                   = train.Hdist,
                          site.Dist               = train.site.Dist,
                          Nt                      = dim(train.Y_ts)[1],
                          n                       = dim(train.Y_ts)[2],
                          nPy                     = py,
                          Px                      = ifelse(!is.null(Fixed.effect.Data[[py]]$X_ts), dim(train.X_ts)[[1]], 0),
                          Pr                      = ifelse(!is.null(Fixed.effect.Data[[py]]$R_ts), dim(train.R_ts)[[1]], 0),
                          Pz                      = ifelse(!is.null(Fixed.effect.Data[[py]]$Z_ts), dim(train.Z_ts)[[1]], 0),
                          Pg                      = ifelse(!is.null(Fixed.effect.Data[[py]]$G_ts), dim(train.G_ts)[[1]], 0),
                          sPx                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sX_ts), dim(train.sX_ts)[[1]], 0),
                          sPr                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sR_ts), dim(train.sR_ts)[[1]], 0),
                          sPz                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sZ_ts), dim(train.sZ_ts)[[1]], 0),
                          sPg                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sG_ts), dim(train.sG_ts)[[1]], 0),
                          index.y                 = train.index.y,
                          nKnots                  = sum(G.basic.data[[py]]$Grid.infor$summary$Knots.count),
                          # sub.group               = G.basic.data[[py]]$sub.group,
                          mesh                    = G.basic.data[[py]]$mesh,
                          Max.Distance.Point.Grid = G.basic.data[[py]]$Max.Distance.Point.Grid,
                          Ch                      = G.basic.data[[py]]$Ch,
                          Sub.Varying.Ch          = G.basic.data[[py]]$Sub.Varying.Ch)


    test.r[[py]] <- list(Y_ts                    = test.Y_ts,
                         ini.Miss_ts             = test.ini.Miss_ts,
                         X_ts                    = test.X_ts,
                         R_ts                    = test.R_ts,
                         Z_ts                    = test.Z_ts,
                         G_ts                    = test.G_ts,
                         sX_ts                   = test.sX_ts,
                         sR_ts                   = test.sR_ts,
                         sZ_ts                   = test.sZ_ts,
                         sG_ts                   = test.sG_ts,
                         spCoordinates           = unique(Pred.coords[as.character(Pred.coords[, Object]) == Obj, ]),
                         Hs                      = test.Hs,
                         Hs.test.indx            = Hs.test.indx,
                         As                      = test.A,
                         Grid.infor              = G.basic.data[[py]]$Grid.infor,
                         Hdist                   = test.Hdist,
                         site.Dist               = test.site.Dist,
                         Nt                      = dim(test.Y_ts)[1],
                         n                       = dim(test.Y_ts)[2],
                         nPy                     = py,
                         Px                      = ifelse(!is.null(Fixed.effect.Data[[py]]$X_ts), dim(test.X_ts)[[1]], 0),
                         Pr                      = ifelse(!is.null(Fixed.effect.Data[[py]]$R_ts), dim(test.R_ts)[[1]], 0),
                         Pz                      = ifelse(!is.null(Fixed.effect.Data[[py]]$Z_ts), dim(test.Z_ts)[[1]], 0),
                         Pg                      = ifelse(!is.null(Fixed.effect.Data[[py]]$G_ts), dim(test.G_ts)[[1]], 0),
                         sPx                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sX_ts), dim(test.sX_ts)[[1]], 0),
                         sPr                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sR_ts), dim(test.sR_ts)[[1]], 0),
                         sPz                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sZ_ts), dim(test.sZ_ts)[[1]], 0),
                         sPg                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sG_ts), dim(test.sG_ts)[[1]], 0),
                         index.y                 = test.index.y,
                         nKnots                  = sum(G.basic.data[[py]]$Grid.infor$summary$Knots.count),
                         # sub.group               = G.basic.data[[py]]$sub.group,
                         # mesh = G.basic.data$mesh,
                         Max.Distance.Point.Grid = G.basic.data[[py]]$Max.Distance.Point.Grid,
                         Ch                      = G.basic.data[[py]]$Ch,
                         Sub.Varying.Ch          = G.basic.data[[py]]$Sub.Varying.Ch,
                         Object                  = Obj)
  }
  names(train.r) <- names(test.r) <- names(Fixed.effect.Data)
  return(list(train = train.r, test = test.r, Object = Obj))
}


