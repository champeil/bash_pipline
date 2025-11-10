# this script is for load necessary packages
# author: laojp
# time: 2023.11.28
# position: SYSUCC bioinformatic
# version: 1.0
# usage: process_load_packages(c("packages"))

process_load_packages <- function(packages){
  for (i in packages) {
    if (!require(i, character.only = TRUE)) {
      stop(paste("package: ", i, " is not exits"))
    }
  }
}