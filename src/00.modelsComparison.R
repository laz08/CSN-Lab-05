# Load and install necessary packages
requiredPackages <- c("igraph", "ggplot2", "Matrix", "ggthemes", "gridExtra")

for (pac in requiredPackages) {
    if(!require(pac,  character.only=TRUE)){
        install.packages(pac, repos="http://cran.rstudio.com")
        library(pac,  character.only=TRUE)
    } 
}
rm(pac)
rm(requiredPackages)


# Set WD and load data
wd = getwd()
if(grepl("nora", wd)) {
    setwd("~/Documents/18-19/CSN/LABS/05/src")
} else {
    # Set Piero Working directory
    setwd("")
}
rm(wd)

#################
####  Flags  ####
#################

APPLY_PROFILING = FALSE

LOAD_EXISTING_RUN_BA = TRUE
LOAD_EXISTING_RUN_BA_RAND = TRUE
LOAD_EXISTING_RUN_BA_NO_GROWTH = TRUE



TESTING = TRUE

#################
### Variables ###
#################

t.max = 10000
timestamps <- c(1, 10, 100, 1000) # 4 vertices to track

#################
### Sourcing  ###
#################
source("commonFunctions.R")
source("01.BA_original.R")
source("02.BA_growth_rand_att.R")
source("03.BA_no_growth.R")

if(TESTING){
    theoretic_ki = c()
    m.0 = 2
    for (t in seq(t.max+1)) {
        theoretic_ki = append(theoretic_ki, (m.0*((t.max/t)**0.5)))
    }
}
