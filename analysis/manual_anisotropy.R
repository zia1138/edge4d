data <- read.delim(gzfile("live_data_output/Manual Aspect Ratio.txt.gz"))

N <- dim(data)[1]

even <- seq(1,N,2)
odd <- seq(2,N,2)

data.even <- data[even,]
data.odd <- data[odd,]

aspect <- data.frame(time = data.even$time, len.A = data.even$length, len.B = data.odd$length)

aspect$ratio <- 0
for(i in 1:dim(aspect)[1]) {
  aspect$ratio[i] <- max(aspect$len.A[i], aspect$len.B[i]) / min(aspect$len.A[i], aspect$len.B[i])
}

wilcox.test(aspect$ratio[aspect$time == 850], aspect$ratio[aspect$time == 0])
boxplot(aspect$ratio[aspect$time == 0], aspect$ratio[aspect$time == 850])

pdf("aspect_ratio_manual.pdf", 5,5)
op <- par(mar=(c(3, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
boxplot(aspect$ratio[aspect$time == 0], aspect$ratio[aspect$time == 850], names=c("0 sec", "800 sec"), lwd=2, col=c("gray", "orange"), ylab="cell aspect ratio", outline=F, ylim=c(1,2.0))
par(op)
dev.off()
