r_run_random_walk <- function(P, start, target, numSteps, numRuns)
{
    walker <- matrix(data = NaN, nrow= numRuns, ncol = numSteps)
    arrCount <- rep(0, numRuns)
    for(run in 1:numRuns)
    {
        walker[run,1] <- start
        step <- 2
        for(step in 1:numSteps)
        {
            loc <- walker[run, step - 1]
            p <- P[loc,]
            p <- cumsum(p)
            movement <- runif(n = 1, min = 0, max = 1)
            aux <- p - movement
            walker[run, step] = which(aux >= 0)[1]
            step = step + 1
        }
    }
    
    
    lt_mfpt <- matrix(data = NaN, nrow = nrow(P), ncol = numSteps)
    fraction_arrival <- matrix(data = NaN, nrow = nrow(P), ncol = numSteps)
    numArr <- matrix(data = 0, nrow = numSteps, ncol = numSteps)
    for(target in 1:nrow(P)){
        for(time in 1:numSteps)
        {
            walker_limited <- matrix(walker[,1:time], ncol = time)
            arrived <- walker_limited[which(walker_limited == target, arr.ind = TRUE)[,1],]
            fraction_arrival[target, time] <- length(unique(which(walker_limited == target, 
                                                                  arr.ind = TRUE)[,1])) / numRuns
            if(length(arrived) > 0)
            {
                if(is.null(ncol(arrived)))
                {
                    lt_mfpt[target, time] <- min(which(arrived == target, arr.ind = TRUE))
                    numArr[target, time] <- length(which(arrived == target, arr.ind = TRUE))
                }
                else
                {
                    lt_mfpt[target, time] <- mean(apply(arrived, 1, 
                                                        function(x) min(which(x == target, arr.ind = TRUE))))
                    numArr[target, time] <- mean(apply(arrived, 1, 
                                                  function(x) length(which(x == target, arr.ind = TRUE))))
                }
            }
        }
    }
    
    return(list(fraction_arrival=fraction_arrival, lt_mfpt= lt_mfpt , numArr = numArr))
}