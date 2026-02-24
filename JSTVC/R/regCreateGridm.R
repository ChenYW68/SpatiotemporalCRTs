regCreateGridm <- function(sample.coords,
                           adjacent.matrix,
                           loc            = c("LON_X", "LAT_Y"),
                           pred.coords    = NULL,
                           BAUs.Dist      = NULL,
                           H.Grid_dist    = NULL,
                           R.sqrt         = 0,
                           site.id        = "ID",
                           sub.group      = "sub.group",
                           ch             = 0.1,
                           method         = c("PL", "Wendland"),
                           var.covariable = NULL)
{
  # sample.coords <- setorderv(sample.coords, c(site.id))
  # if(c("x") %nin% colnames(sample.coords)){
  #   setDF(sample.coords);setDF(obs.data)
  #   sample.coords <- spCoords.transform(sample.coords, method = distance.method);
  #   obs.data <- spCoords.transform(obs.data, method = distance.method);
  # }
  loc.ind <- which(colnames(sample.coords) %in% loc);
  setDF(sample.coords)
  if(is.null(pred.coords)){
    pred.coords <- sample.coords
    setcolorder(pred.coords, loc)
  }
  setDF(pred.coords)
  pred.loc.ind <- which(colnames(pred.coords) %in% loc);

  grid.coords <- sample.coords#[, loc.ind]
  nKnots <- nrow(grid.coords)

  # setnames(grid.coords, 1:2, c("LON_X", "LAT_Y"))
  ###############################################################

  xlim <- range(grid.coords$LON_X)
  ylim <- range(grid.coords$LAT_Y)


  n.clust <- 2^R.sqrt

  if(n.clust  <= 1){
    grid.coords$cluster <- 1

    part.centroids <- data.frame(Var1 = mean(xlim), Var2 = mean(ylim))

  }else{
    # R.sqrt <- 2
    # source("./HDCM/R/Irregblock.R")
    points 	<- data.frame(x = grid.coords$LON_X, y = grid.coords$LAT_Y)
    tree 	<- kdtree(points)
    treenew	<- tree[1:(n.clust - 1), ]
    blocks 	<- NULL
    # cat("......\n")
    # print(treenew)
    # print(head(as.matrix(grid.coords)))
    blocks 	<- kdtree_blocks(treenew, n.blocks = n.clust,
                             loc = as.matrix(grid.coords))

    # blocks
    grid.coords$cluster <- floor(blocks[, 1])

    part.centroids <- grid.coords %>% ddply(.(cluster) #, x, y
                                            , plyr::summarize
                                            , Var1 = mean(LON_X, na.rm = T)
                                            , Var2 = mean(LAT_Y, na.rm = T)
    )

  }



  ###############################################################
  Grid.infor <- list()
  Grid.infor$level <- list()
  Grid.infor$summary <- list()
  Grid.infor$summary$res <- n.clust
  group <- vector(); # g.index <- list()
  nKnots <- 0;
  Hdist <- vector()
  for(g in 1:n.clust){
    c1 <- which(grid.coords$cluster == g)

    Grid.infor$level[[g]] <- list()
    Grid.infor$level[[g]][["Adj.Mat"]] <- spam::as.spam(-adjacent.matrix[c1, c1])

    for(j in 1:length(c1)){
      # Grid.infor$level[[g]]$Adj.Mat[j, ] <- Grid.infor$level[[g]]$Adj.Mat[j, ]/-sum(Grid.infor$level[[g]]$Adj.Mat[j, ])
      Grid.infor$level[[g]][["Adj.Mat"]][j, j] <- ifelse(-sum(Grid.infor$level[[g]]$Adj.Mat[j, ]) != 0,
                                                         -sum(Grid.infor$level[[g]]$Adj.Mat[j, ]), 0)
    }

    Grid.infor$level[[g]][["latCoords"]] <- grid.coords[c1, ]

    Grid.infor$level[[g]]$BAUs.Dist <- BAUs.Dist
    Grid.infor$level[[g]]$var.dist <- rdist(scale(var.covariable))


    # X <- c(ch, 10)
    # mu.alpha <- c(5, rep(-1, ncol(var.covariable) - 1))
    # Adj.Mat <- exp(-sqrt(BAUs.Dist^2/exp(X[1])^2 + rdist(scale(var.covariable))^2/exp(X[2])^2))
    # mu.var <- exp(as.matrix(cbind(1, as.matrix(var.covariable[, -ncol(var.covariable)])))%*% mu.alpha)
    # # mu.var <- exp(mu.alpha[1] + var.covariable[, 1] * mu.alpha[-1])
    # C <- spam::as.spam(t(Adj.Mat * as.vector(mu.var^(0.5))) * as.vector(mu.var^(0.5)))
    # H <-   C%*% solve(C)

    Grid.infor$level[[g]]$Max.Dist <- max(Grid.infor$level[[g]]$BAUs.Dist)
    # Grid.infor$level[[g]]$BAUs.Dist <- Grid.infor$level[[g]]$BAUs.Dist#/Grid.infor$level[[g]]$Max.Dist




    # Grid.Coord.Cons <- rbind(Grid.Coord.Cons, temp)
    nKnots <- nKnots + length(c1)
    group[g] <- length(c1)

    if(g == 1){
      # Grid.infor$summary$Adj.Mat <- Matrix::bdiag(grid_coords[[g]]$Adj.Mat)
      Grid.infor$summary$g.index <- list()
      Grid.infor$summary$g.index[[g]] <- c1
      Grid.infor$summary$Hdist <- H.Grid_dist
      # Grid.infor$summary$BAUs.Dist <-  grid_coords[[g]]$BAUs.Dist
      Hdist[1] <- max(Grid.infor$summary$Hdist)
    }else{
      Grid.infor$summary$g.index[[g]] <- c1#(max(Grid.infor$summary$g.index[[g - 1]]) + 1):nKnots
      # Grid.infor$summary$Adj.Mat <- Matrix::bdiag(Grid.infor$summary$Adj.Mat,
      #                                              grid_coords[[g]]$Adj.Mat)
      hs <- H.Grid_dist
      Hdist[g] <- max(hs)
      Grid.infor$summary$Hdist <- cbind(Grid.infor$summary$Hdist, hs)
      # Grid.infor$summary$BAUs.Dist <- Matrix::bdiag(Grid.infor$summary$BAUs.Dist,
      #                                                grid_coords[[g]]$BAUs.Dist)
    }
    Grid.infor$level[[g]][["nKnots"]] <- length(c1)

  }
  Grid.infor$summary$nKnots <- nKnots



  ###########################################################################
  # if(is.null(thresh)){
  #   if(!is.null(ch)){
  #     thresh <- max(Grid.infor$summary$Hdist)*ch
  #   }
  # }

  threshold <- vector()

  ns <- nrow(pred.coords)

  ###########################################################################
  #                                 Create H
  ###########################################################################
  setDF(pred.coords)
  # if(method == "PL"){
  #   ###########################################################################
  #   #                                create  Hs
  #   ###########################################################################
  #
  #   Hs <- as(as.matrix(inla.spde.make.A(mesh = mesh,
  #                                       loc = as.matrix(pred.coords[, c("LON_X", "LAT_Y")]))),
  #            "sparseMatrix")
  #   rownames(Hs) <- as.character(pred.coords[, which(colnames(pred.coords) == site.id)])
  # }
  Grid.infor$var.covariable <- (var.covariable) #[, -ncol(var.covariable)]
  if(method == "Wenland"){
    Hs.temp <- matrix(0, nrow = ns, ncol = nKnots)
    nn <- 50
    for(g in Grid.infor$summary$res:1){
      threshold[g] <- ch#Hdist[g]*ch
      # row.min <- min(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]])
      #  row.min <- max(apply(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]], 1, min))
      # if(row.min == 0){
      #   row.min <- unique(sort(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]]))[2]
      # }
      #    # apply(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]], 1, min)
      #  threshold[g] <- ifelse(threshold[g] < row.min, 1.1*row.min, threshold[g])
      #  Hs.temp[, Grid.infor$summary$g.index[[g]]] <-
      #    fields::Wendland(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]],
      #             theta = threshold[g],
      #             dimension = 1, k = 1)
      row.min <- (apply(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]], 1, min))
      for(i in 1:length(row.min)){
         Hs.temp[i, Grid.infor$summary$g.index[[g]]] <-
          # fields::Wendland(Grid.infor$summary$Hdist[i, Grid.infor$summary$g.index[[g]]],
          #                  theta = ifelse(threshold[g] < row.min[i], row.min[i] + 1e-5, threshold[g]),
          #                  dimension = 1, k = 1)
        exp(-Grid.infor$summary$Hdist[i, Grid.infor$summary$g.index[[g]]]/ch)
      }
      #1. only setting with Hs.temp/apply(Hs.temp, 1, sum) not good
      #2. only sweep(gg, 2, colMeans(gg), "-"): larger variance
      #3. both: smallest variance
      gg <- Hs.temp/apply(Hs.temp, 1, sum)
      # gg <- sweep(gg, 2, colMeans(gg), "-")
      cat(paste0("Tapering spatial distance: ", round(threshold[g], 4), " (km). \n"))

    }
    id.x <- pred.coords[, which(colnames(pred.coords) == site.id)]

    # gg <- adjacent.matrix
    # rownames(gg) <- NULL
    # diag(gg) <- 1
    # gg <- gg[id.x, ]*Hs.temp
    # gg <- gg/apply(gg, 1, sum)

    # Hs.temp <- adjacent.matrix



    Grid.infor$Hs <- as(gg, "sparseMatrix")#as(Hs.temp, "sparseMatrix") #
    rownames(Grid.infor$Hs) <- as.character(id.x) #colnames(Grid.infor$Hs) <-
  }

  setDF(pred.coords)

  rownames(Grid.infor$summary$Hdist) <-
    as.character(pred.coords[, which(colnames(pred.coords) == site.id)])

  Grid.infor$summary$Knots.count <- group
  setDT(pred.coords)
  setDT(grid.coords)

  sub.group.ind <- which(colnames(sample.coords) %in% sub.group);
  sub.group <- NA
  if(length(sub.group.ind) > 0){
    setDF(Grid.Coords)
    sub.group <- Grid.Coords[, sub.group.ind]
    setDT(Grid.Coords)
  }

  # if(!is.null(var.covariable)){
  #   var.covariable
  # }

  # Grid.infor$var.dist <- rdist(scale(var.covariable))
  # Grid.infor$knots.nn.indx <- sample.nn.indx
  # Grid.infor$knots.dist <- rdist(coords[, 1:2])/1e3
  # Grid.infor$ord <- ord
  # Grid.infor$coords <-coords

  ###########################################################################
  Data_Str <- list(   Pred.coords               = pred.coords
                      , all.Grid.coords         = grid.coords
                      , Grid.infor              = Grid.infor
                      , sub.group               = sub.group
                      , Sub.Varying.Ch          = threshold
                      , Max.Distance.Point.Grid = Hdist
                      , Ch                      = ch)
  ##########################################################################
  return(Data_Str)
}
