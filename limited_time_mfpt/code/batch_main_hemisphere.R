library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(modMax)

rm(list=ls())
graphics.off()


path.lmftp <- paste(getwd(), "limited_time_mfpt", sep = "/")
path.data <- paste(path.lmftp,"data", "subject", sep = "/")
path.code <- paste(path.lmftp, "code", sep = "/")
path.result <- paste(path.lmftp, "results", sep = "/")
path.csv <- paste(path.result, "csv", "100s_10000r_2", sep="/")
path.image <- paste(path.result, "images", sep="/")
source(paste(path.code, "run_random_walk.R", sep = "/"))
source(paste(path.code, "cor_matrices.R", sep = "/"))
source(paste(path.code, "heatPlot.R", sep = "/"))
source(paste(path.code, "makesym.R", sep = "/"))
source(paste(path.code, "zscore.R", sep = "/"))
source(paste(path.code, "makedelta.R", sep = "/"))
source(paste(path.code, "test_functions.R", sep = "/"))
source(paste(path.code, "multiple_walkers.R",sep="/"))

FC <- as.matrix(read.csv(paste(path.data, "FC.csv", sep = "/"), header=FALSE))
SCden <- as.matrix(read.csv(paste(path.data, "SCden.csv", sep = "/"), header=FALSE))

N <- ncol(FC)
FC <- FC[(1:(N/2)), (1:(N/2))]
SCden <- SCden[(1:(N/2)), (1:(N/2))]
N <-ncol(FC)
Pden <- SCden / rowSums(SCden)


# Setup run parameters
numSteps <- 100
numRuns <- 10000

# Create a 3d array where each slice is a number of steps 
# Columns are source, and rows are target
fracArr_SCden <- array(data=NaN, dim = c(numSteps, N, N))
mfpt_SCden <- array(data=NaN, dim = c(numSteps, N, N))
nArr_SCden <- array(data=NaN, dim = c(numSteps, N, N))

randMat <- matrix(runif(numRuns * numSteps, min = 0, max = 1), nrow= numRuns, ncol = numSteps)


for(start in 1:N){
    file.csv <- paste(path.csv,paste(paste("walker", start, sep="_"), "csv", sep="."), sep = "/")
    ## Returns the fraction arrival, and the limited time mean first passage time
    mydata <- r_run_random_walk(Pden,start,numSteps,numRuns,randMat,file.csv)
    ## Splits them into 2 columns respectively
    ST_fraction_arrived <- mydata$fraction_arrival
    ST_lt_mfpt <- mydata$lt_mfpt
    ST_nArr <- mydata$numArr
    
    fracArr_SCden[ , start, ] <- ST_fraction_arrived
    mfpt_SCden[ , start, ] <- ST_lt_mfpt
    nArr_SCden[ , start, ] <- ST_nArr
}
Z_mfpt <- zscore(mfpt_SCden)
sym_fracArr_SCden <- makesym(fracArr_SCden, numSteps)
sym_mfpt_SCden <- makesym(mfpt_SCden, numSteps)
sym_nArr_SCden <- makesym(nArr_SCden, numSteps)
sym_z <- makesym(Z_mfpt, numSteps)

delta <- makedelta(mfpt_SCden)

# Gets the correlation between FC and the simulations
corMat <- cor_matrices(numSteps, FC, sym_fracArr_SCden, sym_mfpt_SCden, sym_nArr_SCden, sym_z, delta)

# Test functions between fraction arrival and mfpt
corFun <- test_functions(numSteps, FC, sym_fracArr_SCden, sym_mfpt_SCden)

# Grab index of maximum correlation to send to heatplot 
mfptDen <- which(abs(corMat$cor_FC_mfpt) == max(abs(corMat$cor_FC_mfpt), na.rm=TRUE), arr.ind = TRUE)
fracADen <- which(abs(corMat$cor_FC_FA) == max(abs(corMat$cor_FC_FA), na.rm = TRUE), arr.ind = TRUE)
nDen <- which(abs(corMat$cor_n_scden) == max(abs(corMat$cor_n_scden), na.rm=TRUE), arr.ind = TRUE)
nDen <- which(abs(corMat$cor_n_scden) == max(abs(corMat$cor_n_scden), na.rm=TRUE), arr.ind = TRUE)

multCor<- which(abs(corFun$mult) == max(abs(corFun$mult), na.rm=TRUE), arr.ind = TRUE)



heatStep <- fracADen
heatPlot(FC, "FC", file = paste(path.image, "FC.png", sep = "/"))
heatPlot(mfpt_SCden[heatStep, , ], title = "MFPT using SCDen", 
         file = paste(path.image, "mfptDen.png", sep = "/"))
heatPlot(fracArr_SCden[heatStep, , ], "FracA using SCDen",
         file = paste(path.image, "fracADen.png", sep = "/"))
heatPlot(nArr_SCden[heatStep, , ], "nArr using SCden",
         file = paste(path.image, "nDen.png", sep = "/"))
heatPlot(Z_mfpt[heatStep, , ], "Z using SCden",
         file = paste(path.image, "ZDen.png", sep = "/"))


tmat <- sym_fracArr_SCden[fracADen, ,]
tmat <- tmat[upper.tri(tmat)]
tmat_in <- tmat[!is.infinite(tmat)]
upFC <- FC[upper.tri(FC)]
upFC_in <- upFC[!is.infinite(tmat)]
cor(upFC_in, tmat_in)

## Scatter for the log of fracA at highest correlation
mat <- sym_fracArr_SCden[fracADen,,]
upFA = mat[upper.tri(mat)]
upFA = log(upFA)
remin <- !is.infinite(upFA)
l_fa_sub = which(upFA>-3, arr.ind=T)
fa_l = upFA[l_fa_sub]
fa_l = fa_l[remin]
fc_long = FC[upper.tri(FC)]
fc_l = fc_long[l_fa_sub]
fc_l = fc_l[remin]
line = lm(fc_l~fa_l)
corfa <- cor(fa_l,fc_l, use='pairwise.complete.obs')
qplot(fa_l, fc_l) + geom_smooth(method = "lm") + labs(x = 'Fraction Arrival', y = 'FC',
                                                      title = paste(corfa, "Correlation, with log"))
ggsave(paste(path.image,"cor_fa_log.png", sep = "/"))
#plot(fa_l, y = fc_l, xlab = 'Log of Fraction Arrival', ylab = 'FC')
#abline(line, col = 'blue')



## Scatter for fracA of highest correlation
mat <- sym_fracArr_SCden[fracADen,,]
upFA = mat[upper.tri(mat)]
#upFA = log(upFA)
remin <- !is.infinite(upFA)
l_fa_sub = which(upFA>.05, arr.ind=T)
fa_l = upFA[l_fa_sub]
fa_l = fa_l[remin]
fc_long = FC[upper.tri(FC)]
fc_l = fc_long[l_fa_sub]
fc_l = fc_l[remin]
corfa <- cor(fa_l,fc_l, use='pairwise.complete.obs')
qplot(fa_l, fc_l) + geom_smooth(method = "lm") + labs(x = 'Fraction Arrival', y = 'FC',
                                                      title = paste(corfa, "Correlation, without log"))
ggsave(paste(path.image,"cor_fa_no_log.png",sep = "/"))


## Scatter for mult of highest correlation
mat <- sym_fracArr_SCden[multCor,,] * sym_mfpt_SCden[multCor,,]
upFA = mat[upper.tri(mat)]
#upFA = log(upFA)
remin <- !is.infinite(upFA)
l_fa_sub = which(upFA>.05, arr.ind=T)
fa_l = upFA[l_fa_sub]
fa_l = fa_l[remin]
fc_long = FC[upper.tri(FC)]
fc_l = fc_long[l_fa_sub]
fc_l = fc_l[remin]
corfa <- cor(fa_l,fc_l, use='pairwise.complete.obs')
qplot(fa_l, fc_l) + geom_smooth(method = "lm") + labs(x = 'Fraction Arrival * MFPT', y = 'FC',
                                                      title = paste(corfa, "Correlation, mfpt*fracArr"))
ggsave(paste(path.image,"cor_mult",sep = "/"))


## Calculate SPL and Driftness
graphSC <- graph.adjacency(SCden, mode='undirected',weighted=T)
SPLmat <- shortest.paths(a, algorithm = 'dijkstra')

drift <- spl_drift(SPLmat, sym_mfpt_SCden)


# multiple MFPT
minwalkers = 5
maxwalkers = 10
nSim = 500
library(gtools)
library(doParallel)
library(foreach)
library(abind)
cl <- makeCluster(3)
registerDoParallel(cl)

Files <- list.files(path.csv, full.names = T)
Files <- Files[mixedorder(Files)]
N <- length(Files)

mult_walker_mfpt = NULL#array(data=NaN, dim = c(N,N,maxwalkers, numSteps))
# Read in file name
acomb <- function(...) abind(..., along=4)
ptm <- proc.time()
mult_walker_mfpt = foreach(start = 1:N, .combine='acomb', .multicombine=TRUE) %dopar% {
    source('/run/media/renge/Beta/Documents/Research/limited_time_mfpt/code/calc_mfpt.R')
    start_mult_walker_mfpt= array(data=NA, dim = c(N,length(minwalkers:maxwalkers), numSteps))
    print("read data next")
    data = read.csv(Files[start])
    data = data[,-1]
    data[,1] = 0 # Allows Return time to be calculated
    #nrow(data) is too large
    # Set number of walkers
    for(i in minwalkers:maxwalkers){
        # make vector of walker bin sizes
        partseq = seq(1,nSim,i)
        num = length(partseq) 
        mfptholder = array(data=NA, dim = c(num, numSteps, N))
        # send to calc_mfpt one batch at a time
        for(j in 1:num){
            # set to temp array where column is number of walkers, row is target, page is steps
            mfptholder[j, , ] = calc_mfpt(data[partseq[j]:(partseq[j]+i-1),], N, numSteps)
            
        }
        # take mean of column 
        # set to mult_walker_mfpt[i, , start, ]
        start_mult_walker_mfpt[, i-minwalkers+1, ] = colMeans(aperm(mfptholder, c(1, 2, 3)), na.rm=T)
    }
    return(start_mult_walker_mfpt)
}
stopCluster(cl)
proc.time()-ptm

switch_order = aperm(mult_walker_mfpt, c(4,1,2,3))

mult_cor = matrix(data=NA, nrow=numSteps, ncol=length(minwalkers:maxwalkers))
for(i in 1:(length(minwalkers:maxwalkers))){
    for(j in 1:numSteps){
        mult_cor[j,i] = cor(c(FC), c(switch_order[ , , i, j]),use = "pairwise.complete.obs")
    }
}