source("./LoadPackages/RDependPackages.R")
library(MASS)  # For multivariate normal simulation
library(fields)  # For visualization
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site
load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Site     <- rbind(Ken.Site, Tan.Site)
Score_Time_Series_Data <- rbind(Kenya_Score_Data, Tanzania_Score_Data)

ind      <- list()
ind[[1]] <- which(Site$flag =="Northwestern")
ind[[2]] <- which(Site$flag =="Northeastern")
ind[[3]] <- which(Site$flag =="Southern")
ind[[4]] <- which(Site$flag =="Western")
ind[[5]] <- which(Site$flag =="Eastern")

# Parameters

phi_t <- 0.2
phi_s <- 1e1*1e3
phi_x <- 1

nu     <- 0.5
delta  <- 2
a      <- 0.2
b      <- 0.2

delta0 <- 0.3
delta1 <- -0.6
delta2 <- -0.6

n_temporal    <- 5
temporal_grid <- seq(0, 1, len = n_temporal)

# Da <- Kenya_Score_Time_Series_Data[Kenya_Score_Time_Series_Data$flag %in% "Northwestern",]
# Ken.Site<- Site[Site$flag %in% "Northwestern",]
simData.DataBase <- list()
for(r in 1:5)
{
  distances <- Site[ind[[r]], ]$distances
  Temp <- Site[ind[[r]], c(1, 2)] %>% left_join(Score_Time_Series_Data, by = c("Village_ID", "Study_Arm"))#[Score_Time_Series_Data$Study_Arm %in% c(1, 2, 4), ]
  cov.ind <- which(colnames(Temp) %in% c("Village_ID",
                                         "log.log.mean.Intensity",
                                         "log.log.var.Intensity"))


  reg <- as.character(unique(Temp$flag, na.rm = T))[1]
  Var.variables <- unique(Temp[, cov.ind])


  # Var.variables$mean.Intensity <- exp(Var.variables$var.log.log.Prevalence)
  # Var.variables$var.Intensity <- exp(Var.variables$var.log.log.Prevalence)

  var.x <- exp(delta0 + delta1*Var.variables$log.log.mean.Intensity +
                 delta2*Var.variables$log.log.var.Intensity )


  plot(Var.variables$log.log.mean.Intensity, var.x)

  log_sigma <- function(z, delta0, delta1, delta2) {
    exp(delta0 + delta1 * z[1] + delta2 * z[2])
  }

  # Temporal decay function: Psi_t
  Psi_t <- function(tau, phi_t, a, b) {
   (abs(tau / phi_t)^(2) + 1)^(1) #b / a
    # exp(-tau / phi_t)
  }

  # Spatial dependence function: Psi_s
  Psi_s <- function(u, nu) {
    # gamma_factor <- #2^(1 - nu) / gamma(nu)
    # bessel_factor <- (sqrt(2 * nu) * u)^nu * besselK(sqrt(2 * nu) * u, nu)
    # gamma_factor * bessel_factor
    exp(-u)
  }

  # Full covariance function
  cov_function <- function(t, t_prime, s, s_prime, #x, x_prime,
                           phi_t, phi_s, phi_x,
                           nu, delta, a, b, z, z_prime,
                           delta0, delta1, delta2) {
    tau <- abs(t - t_prime)
    u   <- sqrt(sum((s - s_prime)^2)/phi_s^2) # + sum((x - x_prime)^2)/phi_x^2

    sigma_s       <- log_sigma(z, delta0, delta1, delta2)
    sigma_s_prime <- log_sigma(z_prime, delta0, delta1, delta2)
    Psi_t_val     <- Psi_t(tau, phi_t, a, b)
    Psi_s_val     <- Psi_s(u/(Psi_t_val)^(1e-2), nu)  #/sqrt(Psi_t_val)
    sqrt(sigma_s) * sqrt(sigma_s_prime) *  Psi_s_val *(Psi_t_val)^(-0.1)#*ifelse(sqrt(sum((s - s_prime)^2)) < 40, 1, 0)
  }



  n_spatial <- nrow(Site[ind[[r]], ])
  spatial_grid <- data.frame(Village_ID = Site[ind[[r]], ]$Village_ID,
                             x = Site[ind[[r]], ]$LON_X,
                             y = Site[ind[[r]], ]$LAT_Y,
                             z = Site[ind[[r]], ]$distances) %>%
    left_join(Var.variables, by = "Village_ID")





  # Combine spatial and temporal indices
  st.locations <- expand.grid(spatial_index = 1:nrow(spatial_grid), temporal_index = 1:n_temporal)
  setorderv(st.locations, c("spatial_index", "temporal_index"))
  nt <- nrow(st.locations)
  cov_matrix <- matrix(0, nt, nt)
  for (i in 1:nt) {
    for (j in 1:nt) {
      s       <- as.numeric(spatial_grid[st.locations$spatial_index[i], 2:4])
      s_prime <- as.numeric(spatial_grid[st.locations$spatial_index[j], 2:4])
      # x <- as.numeric(X[st.locations$spatial_index[i], ])
      # x_prime <- as.numeric(X[st.locations$spatial_index[j], ])
      z       <- as.numeric(spatial_grid[st.locations$spatial_index[i], 5:6])
      z_prime <- as.numeric(spatial_grid[st.locations$spatial_index[j], 5:6])
      t       <- temporal_grid[st.locations$temporal_index[i]]
      t_prime <- temporal_grid[st.locations$temporal_index[j]]
      cov_matrix[i, j] <- cov_function(t, t_prime, s, s_prime, #x, x_prime,
                                       phi_t, phi_s, phi_x, nu, delta, a, b,
                                       z,
                                       z_prime,
                                       delta0, delta1, delta2)
    }
  }
  simData.DataBase[[r]] <- list()
  simData.DataBase[[r]]$reg.Site <- data.table(
    Village_ID = as.vector(t(matrix(Site[ind[[r]], ]$Village_ID,
                             nrow = n_spatial,
                             ncol = n_temporal))),
    Study_Arm = as.vector(t(matrix(Site[ind[[r]], ]$Study_Arm,
                                    nrow = n_spatial,
                                    ncol = n_temporal))),
    LON = as.vector(t(matrix(Site[ind[[r]], ]$LON,
                             nrow = n_spatial,
                             ncol = n_temporal))),
    LAT = as.vector(t(matrix(Site[ind[[r]], ]$LAT,
                             nrow = n_spatial,
                             ncol = n_temporal))),
    LON_X = as.vector(t(matrix(Site[ind[[r]], ]$LON_X,
                             nrow = n_spatial,
                             ncol = n_temporal))),
    LAT_Y = as.vector(t(matrix(Site[ind[[r]], ]$LAT_Y,
                             nrow = n_spatial,
                             ncol = n_temporal))),
    flag = as.vector(t(matrix(Site[ind[[r]], ]$flag,
                                    nrow = n_spatial,
                                    ncol = n_temporal))),
    time.scale = rep(seq(0, 1, length = n_temporal), times = n_spatial),
    time = rep((2010 + 1):(n_temporal + 2010), times = n_spatial))
   simData.DataBase[[r]]$reg.cov_matrix <- cov_matrix
}
sim.Cov.Data <- simData.DataBase

for(r in 1:5){
  diag(sim.Cov.Data[[r]]$reg.cov_matrix) <- diag(sim.Cov.Data[[r]]$reg.cov_matrix) #+ 1e-2
  D <- chol(sim.Cov.Data[[r]]$reg.cov_matrix)

  cat(range(diag(sim.Cov.Data[[r]]$reg.cov_matrix)), "\n")

  w <- t(matrix(rnorm(nrow(sim.Cov.Data[[r]]$reg.cov_matrix)), ncol= nrow(sim.Cov.Data[[r]]$reg.cov_matrix))%*% D)
  # simulated_process_matrix <- matrix(w, nrow = n_temporal, ncol = n_spatial)
  #length(unique(sim.Cov.Data[[r]]$reg.Site$Village_ID))
  sim.Cov.Data[[r]]$reg.Site$W_ts <- w[, 1]
  if(r == 1){
    simData.DataBase <- sim.Cov.Data[[r]]$reg.Site
  }else{
    simData.DataBase <- rbind(simData.DataBase, sim.Cov.Data[[r]]$reg.Site)
  }
}
save(sim.Cov.Data, file =  "./data/sim.Cov.Data.rds")
load("./data/sim.Cov.Data.rds")




range(diag(cov_matrix))
range(simData.DataBase$W_ts)
cov_matrix[1:5, 1:5]
# Simulate spatiotemporal process
# set.seed(123)
# mean_values <- rep(0, n)  # Zero mean

# Reshape the result into a spatiotemporal grid



# plot(Ken.Site$LON, Ken.Site$LAT, )

#
# simData.DataBase <- data.table(
#   LON = as.vector(t(matrix(Ken.Site$LON,
#                            nrow = n_spatial,
#                            ncol = n_temporal))),
#   LAT = as.vector(t(matrix(Ken.Site$LAT,
#                            nrow = n_spatial,
#                            ncol = n_temporal))),
#   time.scale = rep(seq(0, 1, length = n_temporal), times = n_spatial),
#   time = rep((2010 + 1):(n_temporal + 2010), times = n_spatial),
#   W_ts = w[, 1])

da <- simData.DataBase[simData.DataBase$LAT > -1,]
m <- range(da$W_ts)
p <- ggplot(da, aes(x = LON, y = LAT, col = W_ts), size = 10) +
  geom_tile() +
  geom_point() + #range(simData.DataBase$W_ts)
  facet_wrap(vars(time), nrow = 1) +
  # scale_color_gradientn(colours = rainbow(10),
  #                       limits = c(floor(m[1]), ceiling(m[2])),
  #                       breaks = round(seq(floor(m[1]), ceiling(m[2]), length = 4), 0)
  # )+
  scale_color_gradient2(low = "blue", mid = "green", high = "red", midpoint = 0,
                        limits = c((m[1]), (m[2])),
                        breaks = round(seq((m[1]), (m[2]), length = 4), 2)) +
  # theme_minimal() +
  theme_light() +
  labs(title = 'Spatiotemporal processes over times',
       x = 'Longitude',
       y = 'Latitude',
       fill = 'Value',
       color = "Effects   ") + # xlab("Longitude") + ylab("Latitude") +
  # labs(color = "Arms") +
  theme(axis.text = element_text(size = 23, colour = "black")
        # ,axis.text.x = element_text(hjust = 0.25, size = 35, colour = "black")
        , axis.title   = element_text(size = 28, colour = "black")
        , legend.title = element_text(size = 25, colour = "black")
        , legend.text  = element_text(size = 25, colour = "black")
        , strip.text   = element_text(size = 25, colour = "black")
        # , legend.title = element_blank()
        , strip.background = element_rect(colour = "grey80", fill = "grey80")
        , legend.background = element_rect(colour = 'transparent', fill = 'transparent')
        , legend.key.width = unit(5,"line")
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , legend.position  =  c("top")##
        # , plot.margin = margin(2, 2, 2, 2, "pt")
        # , legend.margin = margin(0,unit="cm")
        # , legend.margin = margin(10, 10, 10, 10, "pt")
  )
p
range(simData.DataBase$W_ts)
