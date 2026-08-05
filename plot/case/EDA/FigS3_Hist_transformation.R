rm(list = ls())
#-----------------------------------------
source("./LoadPackages/RDependPackages.R")
load("./data/Kenya_Score_Data_r.RData")
load("./data/Tanzania_Score_Data_r.RData")
da1 <- da2 <- rbind(Kenya_Score_Data, Tanzania_Score_Data)

da1$Trans        <- "Original scale"
da2$Prevalence   <- log(-log(1 - da2$Prevalence))
da2$Trans        <- "Double logarithmic scale"
Score_Data       <- rbind(da1[!is.na(da1$Prevalence),], da2[!is.na(da2$Prevalence),])
Score_Data$Trans <- factor(Score_Data$Trans, levels = c("Original scale", "Double logarithmic scale"))
p <- ggplot() +
  geom_histogram(data  = Score_Data %>% filter(Trans == "Original scale"),
                 aes(x = Prevalence, y = ..density../sum(..density..)),
                 binwidth = 0.03,
                 fill     = "skyblue",
                 color    = "black",
                 alpha    = 0.3) +
  geom_histogram(data = Score_Data %>% filter(Trans == "Double logarithmic scale"),
                 aes(x = Prevalence, y = ..density../sum(..density..)),
                 binwidth = 0.25,
                 fill     = "skyblue",
                 color    = "black",
                 alpha    = 0.3) +
  facet_wrap(~Trans, ncol = 2, scales = "free") +
  scale_size_binned() +
  theme_bw()          +
  labs(x = "Outcome: Prevalence", y = "Density") +
  theme(axis.text        = element_text(size = 12, colour = "black"),
        axis.title       = element_text(size = 12, colour = "black"),
        strip.text       = element_text(size = 12, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "none")
ggsave(plot = p, filename = "./figure/FigS3_Hist.pdf", width = 10, height = 5)
