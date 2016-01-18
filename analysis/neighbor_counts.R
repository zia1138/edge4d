#x <- read.delim("072112_02_t11_50_neighbor_counts_noborder.txt")
#hist(x$neighbor.cnt, breaks=8)

pdf("neighbor_counts.pdf")
op <- par(mar=(c(5, 4.4, 1, 1.25) + 0.1), lwd=2, ps=22)
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
y <- read.delim("live_data_output/072112_02_t11_50_neighbor_counts_noborder.txt")
hist(y$neighbor.cnt, col="gray", xlab="number of neighbors", main="", ylab="reconstructed cell shapes", breaks=9)
par(op)
dev.off()

## hist(x$neighbor.cnt, breaks=8, add=T, col="red")

## y <- data.frame(time=x$time, neighbor.cnt=x$neighbor.cnt)



## hist(x$neighbor.cnt[x$time <= 1050])

## boxplot(neighbor.cnts ~ time, data=y, method="jitter")
