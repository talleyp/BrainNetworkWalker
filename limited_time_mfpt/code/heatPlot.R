heatPlot <- function(data, title, file){
    longData <- melt(data)
    myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")
    zp1 <- ggplot(longData,
                  aes(x = Var2, y = Var1, fill = value))
    zp1 <- zp1 + geom_tile()
    zp1 <- zp1 + scale_fill_gradientn(colours = myPalette(100))
    zp1 <- zp1 + scale_x_discrete(expand = c(0, 0))
    zp1 <- zp1 + scale_y_discrete(expand = c(0, 0))
    zp1 <- zp1 + coord_equal()
    zp1 <- zp1 + theme_bw()
    zp1 <- zp1 + ggtitle(title) + labs(x = "target", y = "source")
    ggsave(file)
    hist(data[,1])
    print(zp1)
}