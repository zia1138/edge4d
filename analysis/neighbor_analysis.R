data <- read.delim("live_data_output/072112_02_t11_50_neighbor_analysis.txt")
d <- data[data$traj.len > 10 ,]

pdf("rank_neighbors.pdf")

op <- par(mar=(c(5, 4.4, 1, 1.25) + 0.1), lwd=2, ps=20)
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.

boxplot(d$rank.contact.1 / d$traj.len,
        d$rank.contact.2 / d$traj.len,
        d$rank.contact.3 / d$traj.len,
        d$rank.contact.4 / d$traj.len,
        d$rank.contact.5 / d$traj.len,
        d$rank.contact.6 / d$traj.len,
        d$rank.contact.7 / d$traj.len,
        d$rank.contact.8 / d$traj.len,
        d$rank.contact.9 / d$traj.len,
        d$rank.contact.10 / d$traj.len, outline=T, cex=0.5,
        names=1:10,
        col=c(rep("orange",6), rep("gray", 4)),
        xlab="rank of neighbor based on duration of contact",
        ylab="percent of trajectory neighbor in contact" )
par(op)

dev.off()


pdf("traj_90pcent_contact.pdf")

op <- par(mar=(c(5, 4.4, 1, 1.25) + 0.1), lwd=2, ps=20)
z <- c(dim(d)[1],
       sum(d$rank.contact.1 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.2 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.3 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.4 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.5 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.6 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.7 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.8 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.9 / d$traj.len >= 0.9, na.rm=T),
       sum(d$rank.contact.10 / d$traj.len >= 0.9, na.rm=T)
       )
  
barplot(z, names=c(0:10), col=c("green", rep("orange", 6), rep("gray", 4)), ylab="number of cell trajectories", xlab="rank of neighbor based on duration of contact")

dev.off()
