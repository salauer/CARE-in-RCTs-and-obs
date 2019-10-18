## install packages required to run code

if(!require(pacman)){
    install.packages("pacman", repos='https://cloud.r-project.org')
    library(pacman)
}

pload(dplyr, ggplot2, grid, gridExtra, readr, xtable, doMC, here)
