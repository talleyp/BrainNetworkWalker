r_run_random_walk <- function(P, start, numSteps, numRuns, movement, filename)
{
    N <- nrow(P)
    walker <- matrix(data = NaN, nrow= numRuns, ncol = numSteps)
    for(run in 1:numRuns)
    {
        walker[run,1] <- start
        for(step in 2:numSteps)
        {
            loc <- walker[run, step - 1]
            p <- P[loc,]
            p <- cumsum(p)
            aux <- p - movement[run, step]
            walker[run, step] = which(aux >= 0)[1]
        }
    }
    
    write.csv(walker, file = filename)
    
    lt_mfpt <- matrix(data = NaN, nrow = numSteps, ncol = N)
    fraction_arrival <- matrix(data = NaN, nrow = numSteps, ncol = N)
    numArr <- matrix(data = 0, nrow = numSteps, ncol = N)
    
    for(target in 1:N){
        # Addition of clause such that MFPT and FracArr and narr do not include start condition for the walker matrix
        if(target == start){
            #print("Before")
            twalker = walker
            walker[,1] = 0
        }
        for(time in 1:numSteps)
        {
            walker_limited <- matrix(walker[,1:time], ncol = time)
            arrived <- walker_limited[which(walker_limited == target, arr.ind = TRUE)[,1],]
            fraction_arrival[time, target] <- length(unique(which(walker_limited == target, 
                                                                  arr.ind = TRUE)[,1])) / numRuns
            
            if(length(arrived) > 0)
            {
                if(is.null(ncol(arrived)))
                {
                    lt_mfpt[time, target] <- min(which(arrived == target, arr.ind = TRUE))
                    if(time == 1 && target == start){
                        numArr[time, target] <- 1
                    }
                    else{
                        numArr[time, target] <- length(which(arrived == target, arr.ind = TRUE))
                    }
                }
                else
                {
                    lt_mfpt[time, target] <- mean(apply(arrived, 1, 
                                                        function(x) min(which(x == target, arr.ind = TRUE))), na.RM=T)
                    numArr[time, target] <- mean(apply(arrived, 1, 
                                                  function(x) length(which(x == target, arr.ind = TRUE))), na.RM=T)
                }
            }
        }
        # Add clause for mean return time to set it such that the mfpt and fracArr cannot include the initial step
        # This returns the walker matrix to the original
        if(target == start){
            #print("After")
            walker = twalker
        }
    }
    
    return(list(fraction_arrival=fraction_arrival, lt_mfpt= lt_mfpt , numArr = numArr))
}