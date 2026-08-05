load(file = "./result/case/SCORE_CV_models_RMSE_CRPS_no_ie.RData")
JSTVC_IE <- CV.error
load(file = "./result/case/SCORE_CV_models_RMSE_CRPS_ie.RData")
da <- rbind(data.frame(Criteria = "RMSE",
                       Arm.AR1 = round(mean(CV.error$Arm_AR1.RMSE), 4)*100,
                       Sub.AR1 = round(mean(CV.error$Sub_AR1.RMSE), 4)*100,
                       JSTVC_xi = round(mean(CV.error$JSTVC_xi.RMSE), 4)*100,
                       JSTVC_IE = round(mean(JSTVC_IE$STVC_ie.RMSE), 4)*100,
                       JSTVC = round(mean(CV.error$STVC.RMSE), 4)*100),
            data.frame(Criteria = "CRPS",
                       Arm.AR1 = round(mean(CV.error$Arm_AR1.CRPS), 4)*100,
                       Sub.AR1 = round(mean(CV.error$Sub_AR1.CRPS), 4)*100,
                       JSTVC_xi = round(mean(CV.error$JSTVC_xi.CRPS), 4)*100,
                       JSTVC_IE = round(mean(JSTVC_IE$STVC_ie.CRPS), 4)*100,
                       JSTVC = round(mean(CV.error$STVC.CRPS), 4)*100))


writexl::write_xlsx(
  list(CV = da),
  path = "./result/summary/Table_S14_CV.xlsx"
)
