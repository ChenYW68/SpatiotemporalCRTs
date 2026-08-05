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

Baseline_Data_2011 <- Score_Data %>%
  filter(Year == 2011) %>%
  mutate(Study_Arm = as.character(Study_Arm))

Baseline_2011 <- Baseline_Data_2011 %>%
  group_by(Study_Arm) %>%
  dplyr::summarise(
    Villages = dplyr::n_distinct(Village_ID),
    mu = mean(Prevalence, na.rm = TRUE),
    SD = sd(Prevalence, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(`Baseline prevalence, mean (SD) %` = paste0(round(mu, 3), " (", round(SD, 3), ")")) %>%
  dplyr::select(
    Study_Arm,
    Villages,
    `Baseline prevalence, mean (SD) %`
  ) %>%
  arrange(Study_Arm)

Overall_2011 <- data.frame(
  Study_Arm = "Overall",
  Villages = dplyr::n_distinct(Baseline_Data_2011$Village_ID),
  check.names = FALSE
)

Overall_2011$`Baseline prevalence, mean (SD) %` <- paste0(
  round(mean(Baseline_Data_2011$Prevalence, na.rm = TRUE), 3),
  " (",
  round(sd(Baseline_Data_2011$Prevalence, na.rm = TRUE), 3),
  ")"
)

Baseline_2011 <- dplyr::bind_rows(Baseline_2011, Overall_2011)

Arm_Pairs <- utils::combn(sort(unique(Baseline_Data_2011$Study_Arm)), 2, simplify = FALSE)

Baseline_t_tests_2011 <- lapply(Arm_Pairs, function(arms) {
  arm_1_data <- Baseline_Data_2011 %>%
    filter(Study_Arm == arms[1]) %>%
    pull(Prevalence)
  arm_2_data <- Baseline_Data_2011 %>%
    filter(Study_Arm == arms[2]) %>%
    pull(Prevalence)

  test_res <- t.test(arm_1_data, arm_2_data)

  data.frame(
    Cluster_1 = arms[1],
    Cluster_2 = arms[2],
    Mean_Cluster_1 = round(mean(arm_1_data, na.rm = TRUE), 3),
    Mean_Cluster_2 = round(mean(arm_2_data, na.rm = TRUE), 3),
    # Mean_Difference = round(unname(diff(test_res$estimate)), 1),
    # CI_Lower = round(test_res$conf.int[1], 1),
    # CI_Upper = round(test_res$conf.int[2], 1),
    P_Value = signif(test_res$p.value, 3),
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

Villages_by_arm_year <- Score_Data %>%
  mutate(Study_Arm = as.character(Study_Arm)) %>%
  group_by(Study_Arm, Year) %>%
  dplyr::summarise(
    Villages = dplyr::n_distinct(Village_ID),
    .groups = "drop"
  ) %>%
  arrange(Study_Arm, Year) %>%
  pivot_wider(
    names_from = Year,
    values_from = Villages
  )


writexl::write_xlsx(
  list(
    Pairwise_t_tests_2011 = as.data.frame(Baseline_t_tests_2011),
    Baseline_2011_by_arm = as.data.frame(Baseline_2011),
    Villages_by_arm_year = as.data.frame(Villages_by_arm_year)
  ),
  path = "./result/summary/Table_S1_Baseline.xlsx"
)

Baseline_2011
Baseline_t_tests_2011
Villages_by_arm_year
