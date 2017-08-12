spl_drift <- function(splobj, mfpt){
    rows = dim(mfpt)[2]
    cols = dim(mfpt)[3]
    steps = dim(mfpt)[1]
    drift <- array(data=NaN, dim = c(steps, cols, rows))
    for(i in 1:steps){
        drift[i,,] <- splobj / mfpt[i,,]
    }
    
    return(drift)
}