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
path.csv <- paste(path.result, "csv", sep="/")
path.image <- paste(path.result, "images", sep="/")
source(paste(path.code, "run_random_walk.R", sep = "/"))
source(paste(path.code, "cor_matrices.R", sep = "/"))
source(paste(path.code, "heatPlot.R", sep = "/"))
source(paste(path.code, "makesym.R", sep = "/"))


FC <- as.matrix(read.csv(paste(path.data, "FC.csv", sep = "/"), header=FALSE))
SCden <- as.matrix(read.csv(paste(path.data, "SCden.csv", sep = "/"), header=FALSE))
SCfa <- as.matrix(read.csv(paste(path.data, "SCfa.csv", sep = "/"), header=FALSE))

N <- ncol(FC)
Pden <- SCden / rowSums(SCden)
Pfa <- SCfa / rowSums(SCfa)


# Setup run parameters
numSteps <- 100
numRuns <- 1000

# Create a 3d array where each slice is a number of steps 
# Columns are source, and rows are target
fracArr_SCden <- array(data=NaN, dim = c(numSteps, N, N))
mfpt_SCden <- array(data=NaN, dim = c(numSteps, N, N))
nArr_SCden <- array(data=NaN, dim = c(numSteps, N, N))
fracArr_SCfa <- array(data=NaN, dim = c(numSteps, N, N))
mfpt_SCfa <- array(data=NaN, dim = c(numSteps, N, N))
nArr_SCfa <- array(data=NaN, dim = c(numSteps, N, N))

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
for(start in 1:N){
    file.csv <- paste(path.csv,paste(paste("fa","walker", start, sep="_"), "csv", sep="."), sep = "/")
    ## Returns the fraction arrival, and the limited time mean first passage time
    mydata <- r_run_random_walk(Pfa,start,numSteps,numRuns,randMat,file.csv)
    ## Splits them into 2 columns respectively
    ST_fraction_arrived <- mydata$fraction_arrival
    ST_lt_mfpt <- mydata$lt_mfpt
    ST_nArr <- mydata$numArr
    
    fracArr_SCfa[ , start, ] <- ST_fraction_arrived
    mfpt_SCfa[ , start, ] <- ST_lt_mfpt
    nArr_SCfa[ , start, ] <- ST_nArr
}
sym_fracArr_SCden <- makesym(fracArr_SCden, numSteps)
sym_fracArr_SCfa <- makesym(fracArr_SCfa, numSteps)
sym_mfpt_SCden <- makesym(mfpt_SCden, numSteps)
sym_mfpt_SCfa <- makesym(mfpt_SCfa, numSteps)
sym_nArr_SCden <- makesym(nArr_SCden, numSteps)
sym_nArr_SCfa <- makesym(nArr_SCfa, numSteps)

# Gets the correlation between FC and the simulations
corMat <- cor_matrices(numSteps, FC, sym_fracArr_SCden, sym_mfpt_SCden, 
                       sym_fracArr_SCfa, sym_mfpt_SCfa,
                       sym_nArr_SCfa, sym_nArr_SCden)

# Grab index of maximum correlation to send to heatplot 
mfptDen <- which(abs(corMat$cor_FC_mfpt) == max(abs(corMat$cor_FC_mfpt), na.rm=TRUE), arr.ind = TRUE)
fracADen <- which(abs(corMat$cor_FC_FA) == max(abs(corMat$cor_FC_FA), na.rm = TRUE), arr.ind = TRUE)
mfptFa <- which(abs(corMat$cor_FC_mfpt_scfa) == max(abs(corMat$cor_FC_mfpt_scfa), na.rm=TRUE), arr.ind = TRUE)
fracAFa <- which(abs(corMat$cor_FC_FA_scfa) == max(abs(corMat$cor_FC_FA_scfa), na.rm=TRUE), arr.ind = TRUE)
nFa <- which(abs(corMat$cor_n_scfa) == max(abs(corMat$cor_n_scfa), na.rm=TRUE), arr.ind = TRUE)
nDen <- which(abs(corMat$cor_n_scden) == max(abs(corMat$cor_n_scden), na.rm=TRUE), arr.ind = TRUE)

heatPlot(FC, "FC", file = paste(path.image, "FC.png", sep = "/"))
heatPlot(mfpt_SCden[mfptDen, , ], title = "MFPT using SCDen", 
         file = paste(path.image, "mfptDen.png", sep = "/"))
heatPlot(fracArr_SCden[fracADen, , ], "FracA using SCDen",
         file = paste(path.image, "fracADen.png", sep = "/"))
heatPlot(mfpt_SCfa[mfptFa, , ], "MFPT using SCFa",
         file = paste(path.image, "mfptFA.png", sep = "/"))
heatPlot(fracArr_SCfa[fracAFa, , ], "FracA using SCfa",
         file = paste(path.image, "fracAFA.png", sep = "/"))
heatPlot(nArr_SCfa[nFa, , ], "nArr using SCfa",
         file = paste(path.image, "nFA.png", sep = "/"))
heatPlot(nArr_SCden[nDen, , ], "nArr using SCden",
         file = paste(path.image, "nDen.png", sep = "/"))



greedy(sym_mfpt_SCfa[50,,])
