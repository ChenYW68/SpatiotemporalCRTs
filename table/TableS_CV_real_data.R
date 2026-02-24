load(file = "./result/case/SCORE_CV_models_RMSE_CRPS_ie.RData")
da <- rbind(data.frame(Criteria = "RMSE",
                       Arm.AR1 = mean(CV.error$Arm_AR1.RMSE),
                       Sub.AR1 = mean(CV.error$Sub_AR1.RMSE),
                       JSTVC_xi = mean(CV.error$JSTVC_xi.RMSE),
                       JSTVC = mean(CV.error$STVC.RMSE)),
            data.frame(Criteria = "CRPS",
                       Arm.AR1 = mean(CV.error$Arm_AR1.CRPS),
                       Sub.AR1 = mean(CV.error$Sub_AR1.CRPS),
                       JSTVC_xi = mean(CV.error$JSTVC_xi.CRPS),
                       JSTVC = mean(CV.error$STVC.CRPS)))
da
