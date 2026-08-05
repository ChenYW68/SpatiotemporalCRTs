#*****You may encounter the following error: "Error 1 occurred while building the shared library."
#*****Simply rerun the code to address this problem.

#------------------------------------------------------------------------
#                      ----  Part A: simulation ----
------------------------------------------------------------------------
#**********************************************************************#
#*
# Section 7.2 Spatiotemporal random fields with Gneiting covariance ----
#*
#**********************************************************************#
#-- Section 7.2.1 Sensitivity for ensemble size: Figure S6
source("./R/Simulation/RandomField/rand_Tuning_Ensemble_Size_VB_EnKF_298.R")
#----------------------------------------------------------------------
#-- Sections 7.2.3 and 7.2.4 Gneiting spatiotemporal covariance
source("./R/Simulation/RandomField/rand_Run_parrallel_JSTVC_250.R")
source("./R/Simulation/RandomField/rand_Run_parrallel_JSTVC_298.R")
source("./R/Simulation/RandomField/rand_Run_parrallel_JSTVC_IE_250.R")
source("./R/Simulation/RandomField/rand_Run_parrallel_JSTVC_IE_298.R")
source("./R/Simulation/RandomField/rand_Run_parrallel_Others_250.R")
source("./R/Simulation/RandomField/rand_Run_parrallel_Others_298.R")
#----------------------------------------------------------------------
rm(list = ls())
#-- Section 7.5 for generating results to generate Figures S7 and S9
#------------------------------------------------------------------------
# running VB using MCMC = FALSE; MCMC using MCMC = TRUE
#------------------------------------------------------------------------
# n = 250
MCMC     <- FALSE
source("./R/Simulation/RandomField/rand_single_rand_Run_vb-mcmc_n_250.R")
MCMC     <- TRUE
source("./R/Simulation/RandomField/rand_single_rand_Run_vb-mcmc_n_250.R")

# n = 298
MCMC     <- FALSE
source("./R/Simulation/RandomField/rand_single_rand_Run_vb-mcmc_n_298.R")
MCMC     <- TRUE
source("./R/Simulation/RandomField/rand_single_rand_Run_vb-mcmc_n_298.R")
#**********************************************************************#
#End of Section 7.2


#**********************************************************************#
#*
#                     Section 7.3 Low-rank basis functions ----
#*
#**********************************************************************#
source("./R/Simulation/SmoothedFun/smoothed_Run_parrallel_JSTVC_250.R")
source("./R/Simulation/SmoothedFun/smoothed_Run_parrallel_JSTVC_298.R")
source("./R/Simulation/SmoothedFun/smoothed_Run_parrallel_JSTVC_IE_250.R")
source("./R/Simulation/SmoothedFun/smoothed_Run_parrallel_JSTVC_IE_298.R")
source("./R/Simulation/SmoothedFun/smoothed_Run_parrallel_others_250.R")
source("./R/Simulation/SmoothedFun/smoothed_Run_parrallel_others_298.R")
#**********************************************************************#
rm(list = ls())
#-- Section 7.5 for generating results to generate Figures S8 and S10
#------------------------------------------------------------------------
# running VB using MCMC = FALSE; MCMC using MCMC = TRUE
#------------------------------------------------------------------------
# n = 250
MCMC     <- FALSE
source("./R/Simulation/SmoothedFun/smoothed_single_Run_vb-mcmc_n_250.R")
MCMC     <- TRUE
source("./R/Simulation/SmoothedFun/smoothed_single_Run_vb-mcmc_n_250.R")

# n = 298
MCMC     <- FALSE
source("./R/Simulation/SmoothedFun/smoothed_single_Run_vb-mcmc_n_298.R")
MCMC     <- TRUE
source("./R/Simulation/SmoothedFun/smoothed_single_Run_vb-mcmc_n_298.R")
#**********************************************************************#
#End of Section 7.3


#**********************************************************************#
#*
#            Section 7.4  Misspecified decay function structures ----
#*
#**********************************************************************#
#-- Section 7.4 Misspecified decay function structures
source("./R/Simulation/RandomField_xMisspecified_Decay/rand_Run_parrallel_JSTVC_250.R")
source("./R/Simulation/RandomField_xMisspecified_Decay/rand_Run_parrallel_JSTVC_298.R")
source("./R/Simulation/RandomField_xMisspecified_Decay/rand_Run_parrallel_JSTVC_IE_250.R")
source("./R/Simulation/RandomField_xMisspecified_Decay/rand_Run_parrallel_JSTVC_IE_298.R")
source("./R/Simulation/RandomField_xMisspecified_Decay/rand_Run_parrallel_Others_250.R")
source("./R/Simulation/RandomField_xMisspecified_Decay/rand_Run_parrallel_Others_298.R")
#**********************************************************************#
#End of Section 7.4

#**********************************************************************#
#*
#                           Part A Summary  ----
#*
#*Tables are saved to "./result/summary"
#*Figures are saved to "./figure"
#**********************************************************************#
#Generate Tables S6–S13
source("./table/Table_S6_S13.R")

#Generate figures: Fig S6-S12
source("./plot/simulation/FigS6_tuning_ensemble_size.R")
source("./plot/simulation/FigS7_Gneiting_trace_mcmc_vb.R")
source("./plot/simulation/FigS8_smoothed_trace_mcmc_vb.R")
source("./plot/simulation/FigS9_Gneiting_surface_mcmc_vb.R")
source("./plot/simulation/FigS10_smoothed_surface_mcmc_vb.R")
source("./plot/simulation/FigS11_Gneiting_monitor_convergence.R")
source("./plot/simulation/FigS12_smoothed_monitor_convergence.R")
#**********************************************************************#
#End of all simulation sections



#------------------------------------------------------------------------
#                         Part B: Real data analysis ----
#
#*Tables are saved to "./result/summary"
#*Figures are saved to "./figure"
#------------------------------------------------------------------------

#**********************************************************************#
#*
#                      Part B.1 Data exploration ----
#*
#**********************************************************************#
#Generate figures
source("./plot/case/EDA/Fig2_Heterogneity.R")
source("./plot/case/EDA/FigS1_Arms_Map_cluster.R")
source("./plot/case/EDA/FigS2_Boxplot.R")
source("./plot/case/EDA/FigS3_Hist_transformation.R")
source("./plot/case/EDA/FigS4_Variance_Relationships.R")
#Generating tables for EDA
source("./table/Table_S1_EDA_baseline.R")
source("./table/Table_S2-S4_EDA_real_data.R")
#**********************************************************************#
#End of Part B.1



#**********************************************************************#
#*
#                     Part B.2 Cross-validation ----
#*
#**********************************************************************#
#CV: Cross-validation is computationally intensive
source("./R/Case/Section_4.1_JSTVC_CV_Predictions.R")

# Table S14
source("./table/Table_S14_Real_data.R")# based on results from CV
#**********************************************************************#
#End of Part B.2



#**********************************************************************#
#*
#                      Part B.3 Fitting ----
#*
#**********************************************************************#
#Fitting across different scenarios
source("./R/Case/Section_4.2_(i)_Standard_JSTVC.R")
source("./R/Case/Section_4.2_(ii)_JSTVC_ie.R")
source("./R/Case/Section_4.2_(iii)-JSTVC_mixed_DE_IE.R")
source("./R/Case/Section_4.3_(i)_Exponential_decay_JSTVC.R")
source("./R/Case/Section_4.3_(ii)_Untransformation_JSTVC.R")
source("./R/Case/Section_4.3_(iv)-JSTVC_spIE_vb.R")
source("./R/Case/Section_4.3_(v)-1_Standard_JSTVC_Kenya_vb.R")
source("./R/Case/Section_4.3_(v)-3_Standard_JSTVC_Tanzania_vb.R")
source("./R/Case/Section_4.3_(vi)-3_Delete_Outliers_JSTVC.R")
source("./R/Case/Section_9.3.6-STVC_non_decay_time_IEs.R")

#The following tasks may take a few days due to the large number of MCMC iterations:
source("./R/Case/Section_4.3_(iii)_MCMC_JSTVC.R")
source("./R/Case/Section_4.3_(v)-2_Standard_JSTVC_Kenya_mcmc.R")
source("./R/Case/Section_4.3_(v)-4_Standard_JSTVC_Tanzania_mcmc.R")
#**********************************************************************#
#End of Part B.3



#**********************************************************************#
#*
#                      Part B.4 Summary ----
#*
#**********************************************************************#
#Generating tables/figures
source("./plot/case/CV/Fig4_FigS27_Table_S15_CV_Fitting.R")


#-Assessing DE, IE, and ATE
source("./plot/case/Ranking/0_all_ATE_orginal_scale.R")

#-Rankings
source("./plot/case/Ranking/1_Fig5_Rankings.R")
source("./plot/case/Ranking/2_FigS13_Rankings_three_competing_methods.R")
source("./plot/case/Ranking/3_FigS14_Rankings_JSTVC_variants.R")

source("./plot/case/Fitting/FigS15_RealData_Trace_mcmc_vb_both.R")

source("./plot/case/Ranking/4_FigS16_Rankings_JSTVC_spIEs.R")
source("./plot/case/Ranking/5_FigS17_Rankings_two_Ken.R")
source("./plot/case/Ranking/6_FigS18_Rankings_two_Tan.R")

source("./plot/case/Fitting/FigS19_RealData_Trace_mcmc_vb_Kenya.R")
source("./plot/case/Fitting/FigS20_RealData_Trace_mcmc_vb_Tanzania.R")

# case = 1 for generating Fig 6 and Fig S22
# case = 2 for generating Fig S24
# case = 3 for generating Fig S25
# case = 4 for generating Fig S23
case <- 1
source("./plot/case/Fitting/Fig6_FigS22_S23_S24_S25_DI_samples.R")
case <- 2
source("./plot/case/Fitting/Fig6_FigS22_S23_S24_S25_DI_samples.R")
case <- 3
source("./plot/case/Fitting/Fig6_FigS22_S23_S24_S25_DI_samples.R")
case <- 4
source("./plot/case/Fitting/Fig6_FigS22_S23_S24_S25_DI_samples.R")


source("./plot/case/Ranking/7_FigS21_Rankings_non_decay.R")

source("./plot/case/Fitting/Fig7_Spatiotemporal_patterns.R")
source("./plot/case/Fitting/FigS26_monitor_convergence_realData.R")
#**********************************************************************#
#End of Part B.4







