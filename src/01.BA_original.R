PREFIX = "../data/BA_v_"

#################
### Variables ###
#################
filenamesBA <- c("../data/BA_v_00001.csv", 
                          "../data/BA_v_00010.csv", 
                          "../data/BA_v_00100.csv",
                          "../data/BA_v_01000.csv")

#################
### Functions ###
#################

getTwoRandomNodesPreferentialAtt <- function(sum.kj, m.0, nodeList, n, k){
    
    n = n - 1
    
    if(sum.kj == 0){
        # First iteration, disconnected graph
        pi.k = rep(1, n)
    } else {
        # Compute probs.
        pi.k = (k[1:n]/sum.kj)
    }
    selectedNodes <- sample(nodeList[1:n], m.0, prob = pi.k)
    return(selectedNodes)
}



# TODO: This is hardcoded for only  m.0!! Change!!
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
            selectedNodes <- getTwoRandomNodesPreferentialAtt(sum.kj, m.0, nodeList, n, k)
            
            n1 = selectedNodes[1]
            n2 = selectedNodes[2]
        
            ## Not allowing multiedges and assuring the stub connects to two diff. nodes.
            if(n1 == n2){
                shouldGetRandNodes = TRUE
            } else{
                shouldGetRandNodes = FALSE
            }
        }
       
        # Add edges on adj.mat
        adj.mat[n, n1] <- 1
        adj.mat[n1, n] <- 1
        
        adj.mat[n, n2] <- 1
        adj.mat[n2, n] <- 1
        
        ## Add num. of edges + 1 added on each it to the total sum
        ## Which is the increase of degree for each iteration
        sum.kj <- sum.kj + (m.0 + 1)
        
        ## Increase degree of each vertice
        k[n1] <- k[n1] + 1
        k[n2] <- k[n2] + 1
        k[n] <- k[n] + 2
        
        # Keep track of our 4 vertices
        v.track.1 <- append(v.track.1, k[v.track[1] + n.0])
        v.track.2 <- append(v.track.2, k[v.track[2] + n.0])
        v.track.3 <- append(v.track.3, k[v.track[3] + n.0])
        v.track.4 <- append(v.track.4, k[v.track[4] + n.0])
        
        
        # Increase n
        n = n + 1
        
        # Printing to keep track of the execution :)
        if(t%%100 == 0){
            cat("t: ", t,"\n")
        }
    }
    
    # Save 4 vertices tracking
    saveNodesDegreeOnFile(v.track[1], v.track.1, PREFIX)
    saveNodesDegreeOnFile(v.track[2], v.track.2, PREFIX)
    saveNodesDegreeOnFile(v.track[3], v.track.3, PREFIX)
    saveNodesDegreeOnFile(v.track[4], v.track.4, PREFIX)
    
    return(k)
}



#################
## Running code #
#################

runBA <- function(t.max) {
    
  n.0 <- 3
  m.0 <- 2
  
  start = Sys.time()
  k = generateBarabasiAlbertModel(t.max, n.0, m.0, timestamps)
  end = Sys.time()
  
  elapsedTime = end - start
  cat("Elasped time: ", elapsedTime, "\n")
  return(k)
}



if(!LOAD_EXISTING_RUN_BA){
    
    if(APPLY_PROFILING){
        Rprof(tmp <- tempfile())
        runBA(t.max)
        Rprof()
        summaryRprof(tmp)
    } else {
        final.k = runBA(t.max)
    }
}

table.BA <- loadModelExecutions(filenamesBA)
tail(table.BA)
