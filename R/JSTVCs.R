JSTVC <- function(Fixed.effect.Data,
                 G.basic.data,
                 center.y        = FALSE,
                 scale.y         = FALSE,
                 Object          = "ALL",
                 Obj.Seq         = 1,
                 Tab             = "JSTVC",
                 CV              = TRUE,
                 prior           = NULL,
                 Para.List       = NULL,
                 true.Para.List  = NULL,
                 Database        = list(DSN         = RODBC::odbcConnect("DSN_01", uid = "myname",
                                       pwd          = "mypwd",
                                       believeNRows = FALSE,
                                       case         = "toupper")),
                 verbose         = TRUE,
                 verbose.VB      = FALSE,
                 transf.Response = c("normal"),
                 save.Predict    = F,
                 Ne              = 100,
                 cs              = 0.5,
                 ct              = 1,
                 n.cores         = 1,
                 itMin           = 10,
                 itMax           = 50,
                 tol.real        = 0.001,
                 seed            = 1234,
                 plot            = TRUE,
                 nu              = c(1e-1, 1e-1),
                 var.select      = FALSE,
                 positive        = TRUE,
                 threshold       = c(0.5, 0.5))
{
  call <- match.call()


  if (is.null(Fixed.effect.Data)) {
    stop("Must provide the data for the response and covariates.\n")
  }
  if (is.null(G.basic.data)) {
    stop("Must provide G.basic.data data.\n")
  }

  res.num <- G.basic.data[[1]]$Grid.infor$summary$res
  Py <- length(Fixed.effect.Data)



  Px <- Pg <- Pz <- Pr <- sPx <- sPz <- sPg <- sPr <- vector()
  Px[1] <- Pz[1] <- Pr[1] <- Pg[1] <- 0
  for(py in 1:Py){
    Px[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$X_ts), 0,
                     dim(Fixed.effect.Data[[py]]$X_ts)[1])
    Pz[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$Z_ts), 0,
                     dim(Fixed.effect.Data[[py]]$Z_ts)[1])
    Pr[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$R_ts), 0,
                     dim(Fixed.effect.Data[[py]]$R_ts)[1])
    Pg[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$G_ts), 0,
                     dim(Fixed.effect.Data[[py]]$G_ts)[1])
  }
  Px <- max(Px); Pz <- max(Pz); Pr <- max(Pr); Pg <- max(Pg)

  sPx[1] <- sPz[1] <- sPg[1] <- sPr[1] <- 0
  for(py in 1:Py){
    sPx[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$sX_ts), 0,
                      dim(Fixed.effect.Data[[py]]$sX_ts)[1])
    sPz[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$sZ_ts), 0,
                      dim(Fixed.effect.Data[[py]]$sZ_ts)[1])
    sPg[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$sG_ts), 0,
                      dim(Fixed.effect.Data[[py]]$sG_ts)[1])
    sPr[py] <- ifelse(is.null(Fixed.effect.Data[[py]]$sR_ts), 0,
                      dim(Fixed.effect.Data[[py]]$sR_ts)[1])
  }

  Phi.v.min.max <- rep(NA, 2)
  Phi.v <- rep(NA, 3)

  if(Pr > 0){
    for(py in 1:Py){
      Name <- dimnames(Fixed.effect.Data[[py]]$R_ts)[[1]]
      for(p in 1:Pr){
        prior[[names(Fixed.effect.Data)[[py]]]]$s.RF[[Name[p]]] <- prior[[names(Fixed.effect.Data)[[py]]]]$s.RF[["Initial"]]
        Para.List[[names(Fixed.effect.Data)[[py]]]]$s.RF[[Name[p]]] <- Para.List[[names(Fixed.effect.Data)[[py]]]]$s.RF[["Initial"]]
      }
    }
  }

  for(py in 1:Py){
    if(sPr[py] > 0){
        Name <- dimnames(Fixed.effect.Data[[py]]$sR_ts)[[1]]
        for(p in 1:Pr){
          prior[[names(Fixed.effect.Data)[[py]]]]$s.RF[[paste0(names(Fixed.effect.Data)[py], ".", )]] <- prior[[names(Fixed.effect.Data)[[py]]]]$s.RF[["Initial"]]
          Para.List[[names(Fixed.effect.Data)[[py]]]]$s.RF[[paste0(names(Fixed.effect.Data)[py], ".", Name[p])]] <- Para.List[[names(Fixed.effect.Data)[[py]]]]$s.RF[["Initial"]]
        }
      }
  }


  if(Pg > 0){
    for(py in 1:Py){
      Name <- dimnames(Fixed.effect.Data[[py]]$G_ts)[[1]]
      for(p in 1:Pg){
        prior[[names(Fixed.effect.Data)[[py]]]]$st.RF[[Name[p]]] <- prior[[names(Fixed.effect.Data)[[py]]]]$st.RF[["Initial"]]
        Para.List[[names(Fixed.effect.Data)[[py]]]]$st.RF[[Name[p]]] <- Para.List[[names(Fixed.effect.Data)[[py]]]]$st.RF[["Initial"]]
      }
      Para.List[[names(Fixed.effect.Data)[[py]]]]$st.RF[["Initial"]] <- NULL
    }
  }

  for(py in 1:Py){
    if(sPg[py] > 0){
      Name <- dimnames(Fixed.effect.Data[[py]]$sG_ts)[[1]]
      for(p in 1:sPg[py]){
        prior[[names(Fixed.effect.Data)[[py]]]]$st.sRF[[paste0(names(Fixed.effect.Data)[py], ".", Name[p])]] <- prior[[names(Fixed.effect.Data)[[py]]]]$st.sRF[["Initial"]]
        Para.List[[names(Fixed.effect.Data)[[py]]]]$st.sRF[[paste0(names(Fixed.effect.Data)[py], ".", Name[p])]] <- Para.List[[names(Fixed.effect.Data)[[py]]]]$st.sRF[["Initial"]]
      }
      Para.List[[names(Fixed.effect.Data)[[py]]]]$st.sRF[["Initial"]] <- prior[[names(Fixed.effect.Data)[[py]]]]$st.sRF[["Initial"]] <- NULL
    }
  }

  for(py in 1:Py){
    if(is.null(Para.List[[names(Fixed.effect.Data)[[py]]]]$beta$mu.beta)){
      Para.List[[names(Fixed.effect.Data)[[py]]]]$beta$mu.beta <- matrix(NA, nrow = Px, ncol = 1)
      rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$beta$mu.beta) <- dimnames(Fixed.effect.Data[[py]]$X_ts)[[1]]
    }
    if(is.null(Para.List[[names(Fixed.effect.Data)[[py]]]]$sbeta$mu.beta)){
      Para.List[[names(Fixed.effect.Data)[[py]]]]$sbeta$mu.beta <- matrix(NA, nrow = sPx[py], ncol = 1)
      rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$sbeta$mu.beta) <- dimnames(Fixed.effect.Data[[py]]$sX_ts)[[1]]
    }


    rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$rho.alpha$mu.rho.alpha) <- rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$alpha.tau.sq$mu.alpha.tau.sq) <-
    rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$ini.alpha$mu.ini.alpha) <- rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$ini.alpha.tau.sq$mu.ini.alpha.tau.sq) <-
    c(dimnames(Fixed.effect.Data[[py]]$Z_ts)[[1]])

    rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$self.rho.alpha$mu.rho.alpha) <- rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$self.alpha.tau.sq$mu.alpha.tau.sq) <-
    rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$self.ini.alpha$mu.ini.alpha) <- rownames(Para.List[[names(Fixed.effect.Data)[[py]]]]$self.ini.alpha.tau.sq$mu.ini.alpha.tau.sq) <-
    c(dimnames(Fixed.effect.Data[[py]]$sZ_ts)[[1]])

  }


  if (!CV) {
    data <- list()
    for(py in 1:length(Fixed.effect.Data)){
      if(py == 1){
        index.y <-  1:dim(Fixed.effect.Data[[py]]$Y_ts)[2]
      }else{
        index.y <-  (max(data[[py - 1]]$index.y) + 1):(dim(Fixed.effect.Data[[py]]$Y_ts)[2] + max(data[[py - 1]]$index.y))
      }
      ini.Miss_ts <- NULL
      if(!is.null(Fixed.effect.Data[[py]]$ini.Miss_ts)){
        ini.Miss_ts <- Fixed.effect.Data[[py]]$ini.Miss_ts
      }

      data[[py]] <- list(
        Y_ts                    = Fixed.effect.Data[[py]]$Y_ts,
        ini.Miss_ts             = ini.Miss_ts,
        X_ts                    = Fixed.effect.Data[[py]]$X_ts,
        R_ts                    = Fixed.effect.Data[[py]]$R_ts,
        Z_ts                    = Fixed.effect.Data[[py]]$Z_ts,
        G_ts                    = Fixed.effect.Data[[py]]$G_ts,
        sX_ts                   =  Fixed.effect.Data[[py]]$sX_ts,
        sZ_ts                   =  Fixed.effect.Data[[py]]$sZ_ts,
        sG_ts                   =  Fixed.effect.Data[[py]]$sG_ts,
        spCoordinates           = unique(G.basic.data[[py]]$Pred.coords),
        Hs                      = spam::as.spam(as.matrix(G.basic.data[[py]]$Grid.infor$Hs)),
        As                      = spam::as.spam(as.matrix(G.basic.data[[py]]$Grid.infor$Hs)),
        Grid.infor              = G.basic.data[[py]]$Grid.infor,
        Hdist                   = G.basic.data[[py]]$Grid.infor$summary$Hdist,
        Nt                      = dim(Fixed.effect.Data[[py]]$Y_ts)[1],
        n                       = dim(Fixed.effect.Data[[py]]$Y_ts)[2],
        nPy                     = py,
        Px                      = ifelse(!is.null(Fixed.effect.Data[[py]]$X_ts),  dim(Fixed.effect.Data[[py]]$X_ts)[[1]], 0),
        Pr                      = ifelse(!is.null(Fixed.effect.Data[[py]]$R_ts), dim(Fixed.effect.Data[[py]]$R_ts)[[1]], 0),
        Pz                      = ifelse(!is.null(Fixed.effect.Data[[py]]$Z_ts), dim(Fixed.effect.Data[[py]]$Z_ts)[[1]], 0),
        Pg                      = ifelse(!is.null(Fixed.effect.Data[[py]]$G_ts), dim(Fixed.effect.Data[[py]]$G_ts)[[1]], 0),
        sPx                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sX_ts), dim(Fixed.effect.Data[[py]]$sX_ts)[[1]], 0),
        sPr                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sR_ts), dim(Fixed.effect.Data[[py]]$sR_ts)[[1]], 0),
        sPz                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sZ_ts), dim(Fixed.effect.Data[[py]]$sZ_ts)[[1]], 0),
        sPg                     = ifelse(!is.null(Fixed.effect.Data[[py]]$sG_ts), dim(Fixed.effect.Data[[py]]$sG_ts)[[1]], 0),
        index.y                 = index.y,
        nKnots                  = sum(G.basic.data[[py]]$Grid.infor$summary$Knots.count),
        mesh                    = G.basic.data[[py]]$mesh,
        Max.Distance.Point.Grid = G.basic.data[[py]]$Max.Distance.Point.Grid,
        Ch                      = G.basic.data[[py]]$Ch,
        Sub.Varying.Ch          = G.basic.data[[py]]$Sub.Varying.Ch)
    }
    names(data) <- names(Fixed.effect.Data)

    Tab_Name <- paste0(Tab, "_", if_else(month(Sys.Date()) >
                                           9, as.character(month(Sys.Date())), paste0("0",
                                            as.character(month(Sys.Date())))), "_", if_else(day(Sys.Date()) >
                                                                                                                                        9, as.character(day(Sys.Date())), paste0("0",
                                                                                                                                                                                 as.character(day(Sys.Date())))))
    if ((!is.null(Database)) & (!is.null(Database$DSN))) {
      sqlDrop(Database$DSN, paste0(Tab_Name), errors = F)
      Database$Table = Tab_Name
    }else {
      Database = list(DSN = NULL, Table = Tab_Name)
    }
      cat("Oject(dataset): All data are used to fit JSTVCs ...\n")
      CV.Re <- .VB.augEnKS(data            = data,
                             test            = NULL,
                             prior           = prior,
                             Para.List       = Para.List,
                             true.Para.List  = NULL,
                             center.y        = center.y,
                             scale.y         = scale.y,
                             Database        = Database,
                             Object          = "ALL",
                             save.Predict    = save.Predict,
                             verbose         = verbose,
                             verbose.VB      = verbose.VB,
                             transf.Response = transf.Response,
                             Ne              = Ne,
                             cs              = cs,
                             ct              = ct,
                             n.cores         = n.cores,
                             itMin           = itMin,
                             itMax           = itMax,
                             seed            = seed,
                             tol.real        = tol.real,
                             plot            = plot,
                             nu              = nu,
                             var.select      = var.select,
                             positive        = positive,
                             threshold               = threshold)
    CV.Re$Tuning.parameter$ch = G.basic.data$ch
    CV.Re$Call = call
    class(CV.Re) <- "JSTVCs"
  }else {
    CV.Re <- list()
    for (num in Obj.Seq) {
      cat("....\n")
      GSD <- Partitioning.Dataset(G.basic.data,
                                  Fixed.effect.Data = Fixed.effect.Data,
                                  Object = Object, num = num,
                                  siteid = Fixed.effect.Data[[1]]$siteid)
      Tab_Name <- paste0(Tab, "_", 
                         if_else(day(Sys.Date()) > 9, as.character(day(Sys.Date())),
                                 paste0("0", as.character(day(Sys.Date())))),
                         "_", GSD$Object)
      if ((!is.null(Database)) & (!is.null(Database$DSN))) {
        sqlDrop(Database$DSN, paste0(Tab_Name), errors = F)
        Database$Table = Tab_Name
      }else {
        Database = list(DSN = NULL, Table = Tab_Name)
      }

      R <- .VB.augEnKS(data            = GSD$train,
                         test            = GSD$test,
                         prior           = prior,
                         Para.List       = Para.List,
                         true.Para.List  = NULL,
                         center.y        = center.y,
                         scale.y         = scale.y,
                         Database        = Database,
                         Object          = GSD$Object,
                         verbose         = verbose,
                         verbose.VB      = verbose.VB,
                         transf.Response = transf.Response,
                         Ne              = Ne,
                         save.Predict    = save.Predict,
                         cs              = cs,
                         ct              = ct,
                         n.cores         = n.cores,
                         seed            = seed,
                         itMin           = itMin,
                         itMax           = itMax,
                         tol.real        = tol.real,
                         plot            = plot,
                         nu              = nu,
                         var.select      = var.select,
                         positive        = positive,
                         threshold               = threshold)
      R$Tuning.parameter$ch = G.basic.data$ch
      R$Tab_Name <- Tab_Name
      CV.Re[[num]] <- R
    }
  }
  return(CV.Re)
}