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
    setwd("~/git/csn-labs/03")
} else {
    # Set Pietro Working directory
    setwd("")
}
rm(wd)

#################
### Functions ###
#################

getTwoRandomNodes <- function(m, m.0, nodeList, n){
    
    n = n - 1
    sum.kj = sum(m)
    pi.k = c()
    
    if(sum.kj == 0){
        # First iteration, disconnected graph
        pi.k = rep(1, n)
    } else {
        # Compute probs
        for(i in seq(n)) {
            ki = sum(m[i, ])
            pi.k = append(pi.k, ki/sum.kj)
        }
    }
    selectedNodes <- sample(nodeList[1:n], m.0, prob = (pi.k + 0.1))
    return(selectedNodes)
}


getNodesDegree <- function(adj.mat) {
    
    ks = c()
    
    for (n in seq(nrow(adj.mat))) {
        ki = sum(adj.mat[n, ])
        ks = append(ks, ki)
    }
    return(ks)
}



generateBarabasiAlbertModel <- function(ts, n.0, m.0, timeseriesTimestamps) {
    
    nodeList = seq(ts + n.0)
    matSize = ts + n.0
    adj.mat = Matrix(0, ncol = matSize, nrow = matSize, sparse=TRUE) 
    
    n = n.0 + 1
    
    for(t in seq(ts)){
      selectedNodes <- getTwoRandomNodes(adj.mat, m.0, nodeList, n)
      
      n1 = selectedNodes[1]
      n2 = selectedNodes[2]
      
      adj.mat[n, n1] <- 1
      adj.mat[n1, n] <- 1
      
      adj.mat[n, n2] <- 1
      adj.mat[n2, n] <- 1
      
      n = n + 1
      
      # Printing to keep track of the execution :)
      if(t%%100 == 0){
          cat("t: ", t,"\n")
      }
    }
    return(adj.mat)
}



#################
## Running code #
#################
#################

t.max <- 10000

t.max = 10
n.0 <- 3
m.0 <- 2
timeseriesTimestamps <- c(1, 10, 100, 1000)

adj.mat = generateBarabasiAlbertModel(t.max, n.0, m.0, timeseriesTimestamps)
getNodesDegree(adj.mat)
