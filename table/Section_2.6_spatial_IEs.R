source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
load("./data/Tanzania_Score_Data_r.RData")
Score_Data  <- rbind(Kenya_Score_Data, Tanzania_Score_Data)
# "Southern", "Northwestern", "Northeastern"
#  "Western", "Eastern"

arm <- 3
region <- unique(Score_Data$flag)
Arm <- Score_Data %>%
  filter(Study_Arm %in% c(arm)
         , flag %in% c(region[1:3])
         #, !is.na(Prevalence)
         )



temp1 <- Arm %>% group_by(Village_ID) %>%
  filter(IEs.CWT.sp.Neigh.500.30 > 0) %>%
  dplyr::summarise(
    CWT.sp = mean(IEs.CWT.sp.Neigh.500.50, na.rm = TRUE),  # good practice to handle NAs
    SBT.sp = mean(IEs.SBT.sp.Neigh.500.50, na.rm = TRUE),
    .groups = "drop"
  )


# id.upp <- temp1[temp1$CWT.sp > mean(temp1$CWT.sp), ]
# id.low <- temp1[temp1$CWT.sp <= mean(temp1$CWT.sp), ]

id.upp <- temp1[temp1$SBT.sp > mean(temp1$SBT.sp), ]
id.low <- temp1[temp1$SBT.sp <= mean(temp1$SBT.sp), ]



Arm1.1 <- Arm %>%
  filter(Village_ID %in% c(id.upp$Village_ID)
         , !is.na(Prevalence))

Arm1.2 <- Arm %>%
  filter(Village_ID %in% c(id.low$Village_ID)
         , !is.na(Prevalence))


relative_change.1 <- Arm1.1 %>%
  filter(Year %in% c(2011, 2015)) %>%       # keep only the years of interest
  dplyr::select(Village_ID, Year, Prevalence) %>%  # keep relevant columns
  pivot_wider(names_from = Year, values_from = Prevalence, names_prefix = "Year_") %>%
  mutate(Relative_Change = (Year_2015 - Year_2011)) %>% as.data.frame()

relative_change.2 <- Arm1.2 %>%
  filter(Year %in% c(2011, 2015)) %>%       # keep only the years of interest
  dplyr::select(Village_ID, Year, Prevalence) %>%  # keep relevant columns
  pivot_wider(names_from = Year, values_from = Prevalence, names_prefix = "Year_") %>%
  mutate(Relative_Change = (Year_2015 - Year_2011))%>% as.data.frame()


mean(relative_change.1$Relative_Change)
mean(relative_change.2$Relative_Change)


# t.test(relative_change.1$Relative_Change,
#        relative_change.2$Relative_Change)

wilcox.test(relative_change.1$Relative_Change,
            relative_change.2$Relative_Change, paired = FALSE)





