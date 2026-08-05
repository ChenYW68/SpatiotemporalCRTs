#' Create grid and basis information
#'
#' @param sample.coords Sample coordinate data.
#' @param adjacent.matrix Adjacency matrix for sampled locations.
#' @param loc Coordinate column names.
#' @param pred.coords Prediction coordinate data.
#' @param BAUs.Dist Distance matrix among basic areal units.
#' @param H.Grid_dist Distance matrix between prediction sites and grid nodes.
#' @param R.sqrt Partitioning depth parameter.
#' @param site.id Site identifier column name.
#' @param sub.group Subgroup column name.
#' @param ch Spatial decay or bandwidth parameter.
#' @param var.covariable Covariates used for variance modeling.
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
                           var.covariable = NULL)
{
  loc.ind <- which(colnames(sample.coords) %in% loc);
  setDF(sample.coords)
  if(is.null(pred.coords)){
    pred.coords <- sample.coords
    setcolorder(pred.coords, loc)
  }
  setDF(pred.coords)
  pred.loc.ind <- which(colnames(pred.coords) %in% loc);

  grid.coords <- sample.coords
  nKnots <- nrow(grid.coords)

  ###############################################################

  xlim <- range(grid.coords$LON_X)
  ylim <- range(grid.coords$LAT_Y)


  n.clust <- 2^R.sqrt

  if(n.clust  <= 1){
    grid.coords$cluster <- 1

    part.centroids <- data.frame(Var1 = mean(xlim), Var2 = mean(ylim))

  }else{
    points 	<- data.frame(x = grid.coords$LON_X, y = grid.coords$LAT_Y)
    tree 	<- kdtree(points)
    treenew	<- tree[1:(n.clust - 1), ]
    blocks 	<- NULL
    # cat("......\n")
    blocks 	<- kdtree_blocks(treenew, n.blocks = n.clust,
                             loc = as.matrix(grid.coords))

    # blocks
    grid.coords$cluster <- floor(blocks[, 1])

    part.centroids <- grid.coords %>% ddply(.(cluster)
                                            , plyr::summarize
                                            , Var1 = mean(LON_X, na.rm = T)
                                            , Var2 = mean(LAT_Y, na.rm = T))

  }
  ###############################################################
  Grid.infor <- list()
  Grid.infor$level <- list()
  Grid.infor$summary <- list()
  Grid.infor$summary$res <- n.clust
  group <- vector();
  nKnots <- 0;
  Hdist <- vector()
  for(g in 1:n.clust){
    c1 <- which(grid.coords$cluster == g)

    Grid.infor$level[[g]] <- list()
    Grid.infor$level[[g]][["Adj.Mat"]] <- spam::as.spam(-adjacent.matrix[c1, c1]) # This object was generated here but was not used in JSTVC for this work.

    for(j in 1:length(c1)){
      Grid.infor$level[[g]][["Adj.Mat"]][j, j] <- ifelse(-sum(Grid.infor$level[[g]]$Adj.Mat[j, ]) != 0,
                                                         -sum(Grid.infor$level[[g]]$Adj.Mat[j, ]), 0)
    }

    Grid.infor$level[[g]][["latCoords"]] <- grid.coords[c1, ]

    Grid.infor$level[[g]]$BAUs.Dist <- BAUs.Dist
    Grid.infor$level[[g]]$var.dist <- rdist(scale(var.covariable))

    Grid.infor$level[[g]]$Max.Dist <- max(Grid.infor$level[[g]]$BAUs.Dist)
    nKnots <- nKnots + length(c1)
    group[g] <- length(c1)

    if(g == 1){
      Grid.infor$summary$g.index <- list()
      Grid.infor$summary$g.index[[g]] <- c1
      Grid.infor$summary$Hdist <- H.Grid_dist
      Hdist[1] <- max(Grid.infor$summary$Hdist)
    }else{
      Grid.infor$summary$g.index[[g]] <- c1
      hs <- H.Grid_dist
      Hdist[g] <- max(hs)
      Grid.infor$summary$Hdist <- cbind(Grid.infor$summary$Hdist, hs)
    }
    Grid.infor$level[[g]][["nKnots"]] <- length(c1)

  }
  Grid.infor$summary$nKnots <- nKnots
  ###########################################################################
  threshold <- vector()
  ns <- nrow(pred.coords)
  ###########################################################################
  #                  Initialize H, i.e., W in the paper
  ###########################################################################
  setDF(pred.coords)

  Grid.infor$var.covariable <- (var.covariable)
{
    Hs.temp <- matrix(0, nrow = ns, ncol = nKnots)
    nn <- 50
    for(g in Grid.infor$summary$res:1){
      threshold[g] <- ch
      row.min <- (apply(Grid.infor$summary$Hdist[, Grid.infor$summary$g.index[[g]]], 1, min))
      for(i in 1:length(row.min)){
        Hs.temp[i, Grid.infor$summary$g.index[[g]]] <-
          exp(-Grid.infor$summary$Hdist[i, Grid.infor$summary$g.index[[g]]]/ch)
      }
      gg <- Hs.temp/apply(Hs.temp, 1, sum)
      cat(paste0("Tapering spatial distance: ", round(threshold[g], 4), " (km). \n"))

    }
    id.x <- pred.coords[, which(colnames(pred.coords) == site.id)]
    Grid.infor$Hs <- as(gg, "sparseMatrix")
    rownames(Grid.infor$Hs) <- as.character(id.x)
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
