calc_mfpt <-  function(walker, N, numSteps){  
    lt_mfpt <- matrix(data = NaN, nrow = numSteps, ncol = N)
    for(target in 1:N){
        for(time in 1:numSteps)
        {
            walker_limited <- as.matrix(walker[,1:time], ncol = time)
            arrived <- walker_limited[which(walker_limited == target, arr.ind = TRUE)[,1],]

            if(length(arrived) > 0)
            {
                if(is.null(ncol(arrived)))
                {
                    lt_mfpt[time, target] <- min(which(arrived == target, arr.ind = TRUE))
                }
                else
                {
                    lt_mfpt[time, target] <- min(apply(arrived, 1, 
                                                        function(x) min(which(x == target, arr.ind = TRUE))), na.RM=T)
                }
            }
        }
    }
    return(lt_mfpt)
}