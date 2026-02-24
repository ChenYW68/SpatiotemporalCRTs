source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
Ken.Site <- Site
Ken.G    <- G.mat

load("./data/Tanzania_Score_Data_r.RData")
Tan.Site <- Site
Tan.G    <- G.mat
Site     <- rbind(Ken.Site, Tan.Site)
#-----------------------------------------
region_flags <- c("Northwestern", "Northeastern", "Southern", "Western", "Eastern")
Ken_indices  <- lapply(region_flags[1:3], function(f) which(Ken.Site$flag == f))
Tan_indices  <- lapply(region_flags[4:5], function(f) which(Tan.Site$flag == f))
Score_Data       <- rbind(Kenya_Score_Data, Tanzania_Score_Data)

Score_Data$Prevalence <- Score_Data$Prevalence*100



Temp0 <- Score_Data %>%
  group_by(Study_Arm, Year) %>%
  dplyr::summarise(
    mu = mean(Prevalence, na.rm = TRUE),
    SD = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(Year, Study_Arm, mean_sd) %>%
  pivot_wider(names_from = Year, values_from = mean_sd)
Temp0



# -------------------------
# Table S1: Prevalence
# -------------------------
Temp1 <- Score_Data %>%
  group_by(Study_Arm, flag) %>%
  dplyr::summarise(
    mu = mean(Prevalence, na.rm = TRUE),
    SD = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(flag, Study_Arm, mean_sd) %>%
  pivot_wider(names_from = flag, values_from = mean_sd)

# -------------------------
# Add "Overall" column for each Study_Arm (mean(SD) across all flags)
# -------------------------
overall_arm <- Score_Data %>%
  group_by(Study_Arm) %>%
  dplyr::summarise(
    mu = mean(Prevalence, na.rm = TRUE),
    SD = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(Study_Arm, mean_sd) %>%
  rename(Overall = mean_sd)

Temp1 <- Temp1 %>%
  left_join(overall_arm, by = "Study_Arm")

# -------------------------
# Add "Overall" row for each flag (mean(SD) across all Study_Arm)
# -------------------------
overall_flag <- Score_Data %>%
  group_by(flag) %>%
  dplyr::summarise(
    mu = mean(Prevalence, na.rm = TRUE),
    SD = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(flag, mean_sd) %>%
  pivot_wider(names_from = flag, values_from = mean_sd)

all <- paste0(round(mean(Score_Data$Prevalence, na.rm = TRUE), 1),
              " (", round(sd(Score_Data$Prevalence, na.rm = TRUE), 1), ")")

Temp1 <- rbind(Temp1, cbind(Study_Arm = "Overall", overall_flag, Overall = all))

# -------------------------
# Table S2: Intensity
# -------------------------
Temp2 <- Score_Data %>%
  group_by(Study_Arm, flag) %>%
  dplyr::summarise(
    mu = mean(exp(exp(log.log.mean.Intensity)), na.rm = TRUE),
    SD = sd(exp(exp(log.log.mean.Intensity)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(flag, Study_Arm, mean_sd) %>%
  pivot_wider(names_from = flag, values_from = mean_sd)

# Add "Overall" column
overall_arm2 <- Score_Data %>%
  group_by(Study_Arm) %>%
  dplyr::summarise(
    mu = mean(exp(exp(log.log.mean.Intensity)), na.rm = TRUE),
    SD = sd(exp(exp(log.log.mean.Intensity)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(Study_Arm, mean_sd) %>%
  rename(Overall = mean_sd)

Temp2 <- Temp2 %>%
  left_join(overall_arm2, by = "Study_Arm")

# Add "Overall" row
overall_flag2 <- Score_Data %>%
  group_by(flag) %>%
  dplyr::summarise(
    mu = mean(exp(exp(log.log.mean.Intensity)), na.rm = TRUE),
    SD = sd(exp(exp(log.log.mean.Intensity)), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = paste0(round(mu, 1), " (", round(SD, 1), ")")) %>%
  dplyr::select(flag, mean_sd) %>%
  pivot_wider(names_from = flag, values_from = mean_sd)

all <- paste0(round(mean(exp(exp(Score_Data$log.log.mean.Intensity)), na.rm = TRUE), 1),
              " (", round(sd(exp(exp(Score_Data$log.log.mean.Intensity)), na.rm = TRUE), 1), ")")
Temp2 <- rbind(Temp2, cbind(Study_Arm = "Overall", overall_flag2, Overall = all))


writexl::write_xlsx(as.data.frame(Temp1)[, c(1, 4, 3, 5, 6, 2, 7)], path = "./result/case/Table_S1.xlsx")
writexl::write_xlsx(as.data.frame(Temp2)[, c(1, 4, 3, 5, 6, 2, 7)], path = "./result/case/Table_S2.xlsx")

tab <- round(cor(Score_Data[, c(9:16, 23:28)]), 2)
writexl::write_xlsx(as.data.frame(tab), path = "./Result/Table_S3.xlsx")
















