#' Calculate the branching times and leaf heights of a tree
#' Branching events labelled by a 1, leaves (sampling/death) by a 0
#' 
#' (based on getx() function by Tanja Stadler)
#'
#' @param tree The tree to calculate the times on
#' @param reverse If true return time from the present, otherwise time from the root (default=FALSE)
#' @param includeRoot If true include the root in the time estimation (only has an effect if reverseTime is false) (default=TRUE)
#' 
#' @export
getTreeIntervals <- function(tree, reverse=FALSE, includeRoot=TRUE) {
    
  nodes <- sort(unique(c(tree$edge)))
  ttype <- times <- numeric(length(nodes))
  
  ttype[tree$edge[1,1]] <- 1
  timespoly <- c()
  
  # Internal nodes (branching events)
  for (j in (tree$edge[1,1]+1):length(nodes)) {
    ttype[j] <- 1
    temp     <- which(tree$edge[,2] == j)
    ancestor <- tree$edge[temp,1]
    times[j] <- times[ancestor] + tree$edge.length[temp]	
    
    # Check for polytomies at this edge
    polytomy <- length(which(tree$edge[,1] == j))
    if (polytomy > 2) {
        temptimes <- rep(times[j],polytomy-2)
        timespoly <- c(timespoly,temptimes)
    }
  }
  
  # External nodes (leaves)
  for (j in 1:(tree$edge[1,1]-1)) {
    temp     <- which(tree$edge[,2] == j)
    ancestor <- tree$edge[temp,1]
    times[j] <- times[ancestor] + tree$edge.length[temp]	
  }
  
  times <- c(times,timespoly)
  ttype <- c(ttype,rep(1,length(timespoly)))
  
  if (reverse == TRUE) {
      maxt  <- max(times)
      times <- maxt - times
  } else 
  if (includeRoot == TRUE) {
    times <- times + tree$root.edge
  }
  
  return (cbind(times,ttype))
}

#' Calculate the branching times of a tree
#'  
#' (based on getx() function by Tanja Stadler)
#'
#' @param tree The tree to calculate the times on
#' @param reverse If true return time from the present, otherwise time from the root (default=FALSE)
#' @param includeRoot If true include the root in the time estimation (only has an effect if reverseTime is false) (default=TRUE)
#' 
#' @export
getBranchingTimes <- function(tree, reverse=FALSE, includeRoot=TRUE) {
  times <- getTreeIntervals(tree, reverse, includeRoot)
  return (times[which(times[,2] == 1),1])  
}

#' Calculate the leaf heights of a tree
#' 
#' (based on getx() function by Tanja Stadler)
#'
#' @param tree The tree to calculate the times on
#' @param reverse If true return time from the present, otherwise time from the root (default=FALSE)
#' @param includeRoot If true include the root in the time estimation (only has an effect if reverseTime is false) (default=TRUE)
#' 
#' @export
getLeafTimes <- function(tree, reverse=FALSE, includeRoot=TRUE) {
  times <- getTreeIntervals(tree, reverse, includeRoot)
  return (times[which(times[,2] == 0),1])  
}
