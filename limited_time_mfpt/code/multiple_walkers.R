multiple_walkers <- function(csvpath, maxwalkers, numSteps){
    # csvpath: path to the csv files you want to read
    # maxwalkers: the largest value of particles you want to send at one time
    library(gtools)
    library(doParallel)
    library(foreach)
    cl <- makeCluster(2)
    registerDoParallel(cl)
    
    Files <- list.files(csvpath, full.names = T)
    Files <- Files[mixedorder(Files)]
    N <- length(Files)

    mult_walker_mfpt = array(data=NaN, dim = c(maxwalkers, numSteps, N, N))
    # Read in file name
    list <- foreach(start = 1:N) %dopar% {
        source('/run/media/renge/Beta/Documents/Research/limited_time_mfpt/code/calc_mfpt.R')
        
        data = read.csv(Files[start])
        data = data[,-1]
        data[,1] = 0 # Allows Return time to be calculated
        nSim = 500 #nrow(data) is too large
        # Set number of walkers
        for(i in 1:maxwalkers){
            # make vector of walker bin sizes
            partseq = seq(1,nSim,i)
            num = length(partseq) - 1
            mfptholder = array(data=NaN, dim = c(num, numSteps, N))
            # send to calc_mfpt one batch at a time
            for(j in 1:num){
                # set to temp array where column is number of walkers, row is target, page is steps
                mfptholder[j, , ] = calc_mfpt(data[partseq[j]:(partseq[j+i]-1), ], N, numSteps)
            }
            # take mean of column 
            # set to mult_walker_mfpt[i, , start, ]
            mult_walker_mfpt[i, , start, ] = colMeans(aperm(mfptholder, c(1, 2, 3)), na.rm=T)
        }
        return(start)
    }
    return(mult_walker_mfpt)
}