#!/usr/bin/env Rscript
#################################################################
# Function: Configurate the R packages for Flex Meta-Storms
# Call: Rscript RM_Config.R
# Authors: Mingqian Zhang,Xiaoquan Su
# Updated at Oct 10, 2022
# Bioinformatics Group, College of Computer Science & Technology, Qingdao University
#################################################################

## install necessary libraries
p <- c("optparse","RColorBrewer")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="https://mirrors.aliyun.com/CRAN/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))

cat("**R Packages Configuration Complete**\n")

#  check environment variables
Env <-Sys.getenv("FlexMetaStorms")
if(nchar(Env)<1){
  cat('Please set the environment variable \"FlexMetaStorms\" to the directory\n')
 }

