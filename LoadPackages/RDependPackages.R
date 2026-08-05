# options(rgdal_show_exportToProj4_warnings="none")
# # 1 loading packages ---------------------------------------
packages <- c("data.table",
              "ggplot2",
              "parallel",
              "latex2exp",
              "lubridate",
              "dplyr",
              "INLA",
              "Hmisc",
              "MASS",
              "tidyr",
              "RColorBrewer",
              "progress",
              "fields",
              "Rcpp",
              "writexl",
              "readxl",
              "ggmap",
              "mapproj",
              "spam",
              "scoringutils",
              "rARPACK",
              "mvnfast",
              "MCMCpack"
              )
# ,'MASS'
# 2  library
for(i in 1:length(packages))
{
  if(!lapply(packages[i], require,
             character.only = TRUE)[[1]])
  {
    install.packages(packages[i])
    # library(packages[i])
    lapply(packages[i], require,
           character.only = TRUE)
  }else{lapply(packages[i], require,
               character.only = TRUE)}
}
# x=lapply(packages, require, character.only = TRUE)
# rm(list=ls())
rm(i, packages)


