# Load and install necessary packages
requiredPackages <- c("igraph", "ggplot2", "Matrix")

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
    setwd("~/Documents/18-19/CSN/LABS/05")
} else {
    # Set Piero Working directory
    setwd("")
}
rm(wd)


#################
### Sourcing  ###
#################
source("commonFunctions.R")

APPLY_PROFILING = FALSE
#################
### Functions ###
#################

getTwoRandomNodes <- function(sum.kj, m.0, nodeList, n, k){
    
    n = n - 1
    
    if(sum.kj == 0){
        # First iteration, disconnected graph
        pi.k = rep(1, n)
    } else {
        # Compute probs.
        ## Adding offset so everyone has a chance to be selected
        pi.k = (k[1:n]/sum.kj) + 0.1
    }
    selectedNodes <- sample(nodeList[1:n], m.0, prob = pi.k)
    return(selectedNodes)
}




generateBarabasiAlbertModel <- function(ts, n.0, m.0, v.track) {
    
    n.max = ts + n.0
    nodeList = seq(n.max)   ## Vector representing idx from 1..n.max
    k = rep(0, n.max)       ## Degrees of each node. Start on 0
    sum.kj = 0              ## Total sum of degrees
    
    adj.mat = Matrix(0, ncol = n.max, nrow = n.max, sparse=TRUE) 
    
    n = n.0 + 1
    
    v.track.1 <- c(0)
    v.track.2 <- c(0)
    v.track.3 <- c(0)
    v.track.4 <- c(0)
    
    for(t in seq(ts)){
        
        shouldGetRandNodes <- TRUE
        
        while(shouldGetRandNodes){
            selectedNodes <- getTwoRandomNodes(sum.kj, m.0, nodeList, n, k)
            
            n1 = selectedNodes[1]
            n2 = selectedNodes[2]
        
            ## Not allowing multiedges and assuring the stub connects to two diff. nodes.
            if(n1 != n2 & adj.mat[n, n1] == 0 & adj.mat[n, n2] == 0){
                shouldGetRandNodes = FALSE
            }
        }
        
        selectedNodes <- getTwoRandomNodes(sum.kj, m.0, nodeList, n, k)
        
        n1 = selectedNodes[1]
        n2 = selectedNodes[2]
        
        
        # Add edges on adj.mat
        adj.mat[n, n1] <- 1
        adj.mat[n1, n] <- 1
        
        adj.mat[n, n2] <- 1
        adj.mat[n2, n] <- 1
        
        ## Add num. of edges added on each it to the total sum
        sum.kj <- sum.kj + m.0
        
        ## Increase degree of each vertice
        k[n1] <- k[n1] + 1
        k[n2] <- k[n2] + 1
        
        
        # Keep track of our 4 vertices
        v.track.1 <- append(v.track.1, k[v.track[1]])
        v.track.2 <- append(v.track.2, k[v.track[2]])
        v.track.3 <- append(v.track.3, k[v.track[3]])
        v.track.4 <- append(v.track.4, k[v.track[4]])
        
        
        # Increase n
        n = n + 1
        
        # Printing to keep track of the execution :)
        if(t%%100 == 0){
            cat("t: ", t,"\n")
        }
    }
    
    # Save 4 vertices tracking
    saveNodesDegreeOnFile(v.track[1], v.track.1, "BA_v_")
    saveNodesDegreeOnFile(v.track[2], v.track.2, "BA_v_")
    saveNodesDegreeOnFile(v.track[3], v.track.3, "BA_v_")
    saveNodesDegreeOnFile(v.track[4], v.track.4, "BA_v_")
    
    return(k)
}



#################
## Running code #
#################
#################


runBarabasiAlbertModelConstruction <- function(t.max) {
    
  n.0 <- 3
  m.0 <- 2
  timestamps <- c(1, 10, 100, 1000) # 4 vertices to track
  
  start = Sys.time()
  k = generateBarabasiAlbertModel(t.max, n.0, m.0, timestamps)
  end = Sys.time()
  
  elapsedTime = end - start
  cat("Elasped time: ", elapsedTime, "\n")
}


t.max = 10000

if(APPLY_PROFILING){
    Rprof(tmp <- tempfile())
    runBarabasiAlbertModelConstruction(t.max)
    Rprof()
    summaryRprof(tmp)
} else {
    runBarabasiAlbertModelConstruction(t.max)
}


ks = read.csv2("BA_v_00001.csv")
ks = as.vector(ks)
plot(ks, ylab = "k", xlab = "Num of nodes", main="Degrees distribution")



