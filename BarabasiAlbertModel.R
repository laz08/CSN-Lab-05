library("Matrix")

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




ts = 10000 #TODO: 100000       # Max ts
n.0 <- 3
m.0 = 2

nodeList = seq(ts)
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
    
    if(t%%100 == 0){
        cat("t: ", t,"\n")
    }
}
adj.mat



getNodesDegree <- function(adj.mat) {
  
  ks = c()
  for (n in seq(nrow(adj.mat))) {
      
      ki = sum(adj.mat[n, ])
      ks = append(ks, ki)
  }

  return(ks)
}

getNodesDegree(adj.mat)
