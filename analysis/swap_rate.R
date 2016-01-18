x2 <- read.delim("live_data_output/072112_02_t11_50_neighbor_swaps_membonly.txt")
x2 <- x2[x2$time < 850+950,]
z2 <- t(data.frame(error.rate = x2$swap.events / x2$cells.analyzed, ok.rate = 1 - x2$swap.events / x2$cells.analyzed))


x <- read.delim("live_data_output/072112_02_t11_50_neighbor_swaps.txt")
x <- x[x$time < 850+950,]
z <- t(data.frame(error.rate = x$swap.events / x$cells.analyzed, ok.rate = 1 - x$swap.events / x$cells.analyzed))

pdf("swap_rate_new.pdf", 5, 5)
op <- par(mar=(c(4.3, 4.4, 1, 1.25) + 0.1), ps=18) #,lwd=3)

## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
barplot(x2$swap.events / x2$cells.analyzed, col="gray",
        names.arg = x2$time-850, ylab="Estimate of tracking error rate", ylim=c(0,0.5), legend=F, xlab="Time (seconds)")#, lwd=3)

## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
##barplot(x$swap.events / x$cells.analyzed, col="red",
##        names.arg = x$time-850, ylab="Estimate of tracking error rate", ylim=c(0,0.5), legend=F, xlab="Time (seconds)", lwd=3, add=T)

barplot(x$swap.events / x$cells.analyzed, col="red",
        names.arg = x$time-850, ylim=c(0,0.5), legend=F, add=T) #,lwd=3)

par(op)
dev.off()

## pdf("swap_rate2.pdf", 7, 4)
## op <- par(mar=(c(4.3, 4.4, 1, 1.25) + 0.1), lwd=3, ps=22)
## barplot(z, col=c("red", "green"),
##         names.arg = x$time, ylab="Error Rate", ylim=c(0,1), legend=F, xlab="Time (seconds)", lwd=3)
## dev.off()
