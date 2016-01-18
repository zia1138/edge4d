membrane.only <- read.delim("live_data_output/072112_02_t11_50_trajectory_lengths_membraneonly.txt")
combined <- read.delim("live_data_output/072112_02_t11_50_trajectory_lengths.txt")

pdf("traj_lengths.pdf", 5, 5)

op <- par(mar=(c(3, 4.4, 1, 1.25) + 0.1), lwd=2, ps=22)
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.

boxplot(membrane.only$traj.length,
        combined$traj.length, names=c("membranes", "combined"), lwd=2, col=c("darkgreen", "gray"),
        ylab="number of time steps tracked", notch=T)


par(op)
dev.off()


