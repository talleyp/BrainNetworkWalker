# TO DO: Correlation of matrices with FC

library(igraph)

path.row = paste(path.image, "rowMatrices", sep="/")
path.col = paste(path.image, "colMatrices", sep="/")
path.mix = paste(path.image, "mixMatrices", sep="/")

mfpt100 = mfpt_SCden[100,,]
MixMaxCor = matrix(data=0, nrow=139, ncol=139)
for(i in 1:138){
    rowi = mfpt100[i,]
    for(j in (i+1):139){
        rowj = mfpt100[j,]
        corij = cor(rowi, rowj, use = "pairwise.complete.obs")
        MixMaxCor[i,j] = corij
    }
}

# Column pairwise correlation
for(i in 1:138){
    coli = mfpt100[,i]
    for(j in (i+1):139){
        colj = mfpt100[,j]
        corij = cor(coli, colj, use = "pairwise.complete.obs")
        MixMaxCor[j,i] = corij
    }
}

symMix <- (MixMaxCor[ , ] + t(MixMaxCor[,])) / 2
diag(symMix[,]) <- 0
# sympos = sym
# for(i in 1:139){
#     for(j in 1:139){
#         if(sympos[i,j]<0){
#             sympos[i,j] = 0
#         }
#     }
# }
gMixed = graph_from_adjacency_matrix(sym, mode = "undirected", weighted = TRUE)
clMixed = spinglass.community(gMixed)
MixCl = NULL
for(i in 1:length(cl$csize)){
    MixCl = c(MixCl, cl[i])
}
MixCl = unlist(MixCl, use.names = F)

heatPlot(sym[MixCl,MixCl], "spinglass","spinglass.png")

cor_mix <- rep(NULL, numSteps)
cor_row<- rep(NULL, numSteps)
cor_col <- rep(NULL, numSteps)

for(step in 1:100){
    MixCor = matrix(data=0, nrow=139, ncol=139)
    RowCor = matrix(data=0, nrow=139, ncol=139)
    ColCor = matrix(data=0, nrow=139, ncol=139)
    
    #row pairwise correlation
    for(i in 1:138){
        rowi = mfpt_SCden[step, i,]
        for(j in (i+1):139){
            rowj = mfpt_SCden[step, j,]
            corij = cor(rowi, rowj, use = "pairwise.complete.obs")
            RowCor[i,j] = corij
            RowCor[j,i] = corij
            MixCor[i,j] = corij
        }
    }
    
    # Column pairwise correlation
    for(i in 1:138){
        coli =  mfpt_SCden[step, , i]
        for(j in (i+1):139){
            colj = mfpt_SCden[step, , j]
            corij = cor(coli, colj, use = "pairwise.complete.obs")
            ColCor[i,j] = corij
            ColCor[j,i] = corij
            MixCor[j,i] = corij
        }
    }
    symMix <- (MixCor[ , ] + t(MixCor[,])) / 2
    diag(symMix[,]) <- 0
    titleRow = paste(step, "row")
    fileRow = paste(path.row, paste(step,"_row.png",sep=""), sep="/")
    titleCol = paste(step, "col")
    fileCol = paste(path.col, paste(step,"_col.png",sep=""), sep="/")
    titleMix = paste(step, "mix")
    fileMix = paste(path.mix, paste(step,"_mix.png",sep=""), sep="/")
    
    cor_mix[step] <- cor(c(FC), c(symMix),use = "pairwise.complete.obs")
    cor_row[step] <- cor(c(FC), c(RowCor),use = "pairwise.complete.obs")
    cor_col[step] <- cor(c(FC), c(ColCor),use = "pairwise.complete.obs")
    
    heatPlot(RowCor[MixCl,MixCl], titleRow, fileRow)
    heatPlot(ColCor[MixCl,MixCl], titleCol, fileCol)
    heatPlot(symMix[MixCl,MixCl], titleMix,fileMix)
}

stepVec <- seq(1, numSteps)
cor_matrix = as.data.frame(cbind(stepVec, cor_mix, cor_col, cor_row))
ggplot(data = cor_matrix, aes(x= stepVec, y = value, color = variable)) + 
    geom_point(aes(y = cor_mix, col = "cor_mix")) + 
    geom_line(aes(y = cor_mix, col = "cor_mix")) + 
    geom_point(aes(y = cor_col, col = "cor_col")) +  
    geom_line(aes(y = cor_col, col = "cor_col")) + 
    geom_point(aes(y = cor_row, col = "cor_row")) + 
    geom_line(aes(y = cor_row, col = "cor_row")) + 
    labs(x="Step", y="Correlation")
ggsave(paste(getwd(),"/limited_time_mfpt","/results","/pairwise_cor.png", sep = '/'))
