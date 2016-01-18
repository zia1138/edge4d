library(plyr)

data1.name <- "072112_02_t11_50"
data1.file <- "live_data_output/072112_02_t11_50_all_traj.txt"
time0.data1 <- 850
PF.basal.cells.data1 <- c(0, 32, 38, 41, 47, 70, 79, 112, 117, 120, 144, 124, 152, 156, 183, 195, 202, 203)
PF.1.neighbors.data1 <- c(25, 34, 57, 74, 75, 87, 89, 98, 108, 109, 114, 116, 124, 146, 148, 177, 181, 198, 206, 221, 238)
PF.2.neighbors.data1 <- c(2, 10, 20, 30, 60, 64, 67, 86, 94, 100, 125, 136, 155, 162, 175, 182, 194, 215, 222, 315)
PF.3.neighbors.data1 <- c(16, 24, 33, 52, 56, 78, 93, 97, 107, 111, 119, 129, 132, 135, 171, 180, 186, 192, 200, 218, 235, 259, 268, 269, 270, 323, 325, 326, 335, 337, 338, 356)

AF.basal.cells.data1 <- c(137,110,35,40,36,3,140,141,142,178,26)
AF.1.neighbors.data1 <- c(6, 9, 23, 27, 42, 48, 65, 71, 73, 95, 105, 113, 118, 123, 130, 173, 190, 193, 196, 213)
AF.2.neighbors.data1 <- c(4, 18, 58, 68, 72, 82, 84, 85, 88, 106, 131, 133, 145, 150, 153, 157, 159, 207, 214, 216, 226, 236, 247, 271, 279, 317)
AF.3.neighbors.data1 <- c(15, 17, 22, 55, 59, 121, 128, 134, 151, 167, 197, 204, 205, 207, 208, 216, 224, 232, 254, 264, 271, 279, 289)

# time points 800 and 150 are the cleanest
vol.AB.analysis <- function(name, file, time0, basal.cells, late.time = 800) {
  traj <- read.delim(file) 
  vol.time0 <- c()
  vol.late <- c()
  AB.len.time0 <- c()
  AB.len.late <- c()
  time.late <- time0 + late.time
  for(id in basal.cells) {
    cur.traj <- traj[traj$trajectory.id == id,]
    cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]

    vol.time0 <- c(vol.time0, median(cur.traj$volume[cur.traj$time == time0], na.rm=T) )
    vol.late <- c(vol.late, median(cur.traj$volume[cur.traj$time == time.late], na.rm=T) )

    AB.len.time0 <- c(AB.len.time0, median(cur.traj$AB.length[cur.traj$time == time0], na.rm=T) )
    AB.len.late <- c(AB.len.late, median(cur.traj$AB.length[cur.traj$time == time.late], na.rm=T) )
  }
  ret <- data.frame(basal.cells, vol.time0, vol.late, AB.len.time0, AB.len.late)
  ret$name <- name
  ret
}

Z3 <- vol.AB.analysis(data1.name, data1.file, time0.data1, c(PF.basal.cells.data1, AF.basal.cells.data1))

pdf("volume_conserved.pdf", 5, 5)
op <- par(mar=(c(3, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
boxplot(Z3$vol.time0,Z3$vol.late, ylim=c(500,2000), ylab=expression(volume~microns^3), names=c("0 sec", "800 sec"), lwd=2, col=c("gray", "darkcyan"), notch=F, outline=F)
par(op)
wilcox.test(Z3$vol.time0,Z3$vol.late, paired=T)
dev.off()

pdf("ABlength_changed.pdf", 5, 5)
op <- par(mar=(c(3, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
boxplot(Z3$AB.len.time0,Z3$AB.len.late, ylim=c(25,50), ylab=expression(apical~basal~length~microns), names=c("0 sec", "800 sec"), lwd=2, col=c("gray", "purple"), notch=F, outline=F)
par(op)
dev.off()

collect.data <- function(name, file, time0, basal.cells ) {
  traj <- read.delim(file)

  times <- c()
  AB.length <- c()
  basal.z <- c()
  nuc.apical.dist <- c()
  vol.above.nuc <- c()
  sa.above.nuc <- c()
  
  for(id in basal.cells) {
    cur.traj <- traj[traj$trajectory.id == id,]
    cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]

    times <- c(times, cur.traj$time - time0)
    AB.length <- c(AB.length, cur.traj$AB.length)

    sub.z <- cur.traj$apical.z[which(cur.traj$time == time0)]
    if(length(sub.z) == 0) { basal.z <- c(basal.z, rep(NA,length(cur.traj$basal.z)) )  }
    else { basal.z <- c(basal.z, cur.traj$basal.z - sub.z)  }

    nuc.apical.dist <- c(nuc.apical.dist, cur.traj$nuc.apical.dist)
    vol.above.nuc <- c(vol.above.nuc, cur.traj$vol.above.nuc)
    sa.above.nuc <- c(sa.above.nuc, cur.traj$sa.above.nuc)
  }
  data.frame(times, AB.length, basal.z, nuc.apical.dist, vol.above.nuc, sa.above.nuc)
}



Z0 <- collect.data(data1.name, data1.file, time0.data1, c(AF.basal.cells.data1, PF.basal.cells.data1))
Z1 <- collect.data(data1.name, data1.file, time0.data1, c(AF.1.neighbors.data1, PF.1.neighbors.data1))
Z2 <- collect.data(data1.name, data1.file, time0.data1, c(AF.2.neighbors.data1, PF.2.neighbors.data1))
Z3 <- collect.data(data1.name, data1.file, time0.data1, c(AF.3.neighbors.data1, PF.3.neighbors.data1))

tsteps <- sort(unique(Z0$times))
nuc.dist.mu0 <- c()
for(t in tsteps ) {  nuc.dist.mu0 <- c(nuc.dist.mu0, median(Z0[Z0$times == t,]$nuc.apical.dist, na.rm=T) ) }

tsteps <- sort(unique(Z1$times))
nuc.dist.mu1 <- c()
for(t in tsteps ) {  nuc.dist.mu1 <- c(nuc.dist.mu1, median(Z1[Z1$times == t,]$nuc.apical.dist, na.rm=T) ) }

tsteps <- sort(unique(Z2$times))
nuc.dist.mu2 <- c()
for(t in tsteps ) {  nuc.dist.mu2 <- c(nuc.dist.mu2, median(Z2[Z2$times == t,]$nuc.apical.dist, na.rm=T) ) }

tsteps <- sort(unique(Z3$times))
nuc.dist.mu3 <- c()
for(t in tsteps ) {  nuc.dist.mu3 <- c(nuc.dist.mu3, median(Z3[Z3$times == t,]$nuc.apical.dist, na.rm=T) ) }

pdf("neighbor_nuc_analysis.pdf")
op <- par(mar=(c(4.3, 5, 1, 1.25) + 0.1), lwd=3, ps=22)
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
plot(c(-700, 200), c(8,26), type="n", xlab="time (seconds)", ylab="distance of nucleus from cell apex, microns")
points(tsteps, nuc.dist.mu0, type="b", col="green2")
points(tsteps, nuc.dist.mu1, type="b", col="magenta2")
points(tsteps, nuc.dist.mu2, type="b", col="blue")
points(tsteps, nuc.dist.mu3, type="b", col="red")
dev.off()


traj <- read.delim(data1.file)

pdf("vol_traces.pdf")
for(id in c(PF.basal.cells.data1, AF.basal.cells.data1) ) {
  cur.traj <- traj[traj$trajectory.id == id,]
  cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
  cur.traj$time <- cur.traj$time - 850
  cur.traj <- cur.traj[ 0 <= cur.traj$time & cur.traj$time <= 800,]
  plot(cur.traj$time, cur.traj$volume, ylim=c(500,2000), type="b")
  title(id)
}  
dev.off()


pdf("ABlength_traces.pdf")
for(id in c(PF.basal.cells.data1, AF.basal.cells.data1) ) {
  cur.traj <- traj[traj$trajectory.id == id,]
  cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
  cur.traj$time <- cur.traj$time - 850
  cur.traj <- cur.traj[0 <= cur.traj$time & cur.traj$time <= 800,]
  plot(cur.traj$time, cur.traj$AB.length, ylim=c(25,50), type="b")
  title(id)
}  
dev.off()

pdf("ABlength_example.pdf", 5, 5)
op <- par(mar=(c(4.8, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
plot(c(0,850), c(25,50), type="n", lwd=2, ylab=expression(apical~basal~length~microns), xlab="time (sec)")
for(id in c(120) ) {
  cur.traj <- traj[traj$trajectory.id == id,]
  cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
  cur.traj$time <- cur.traj$time - 850
  cur.traj <- cur.traj[0 <= cur.traj$time & cur.traj$time <= 800,]
  points(loess.smooth(cur.traj$time, cur.traj$AB.length, family="gaussian"), type="l", lwd=3, col="purple")
  points(cur.traj$time, cur.traj$AB.length, ylim=c(25,50), type="p", pch=21, bg="gray", col="black", cex=1.2)
}
par(op)
dev.off()

pdf("volume_example.pdf", 5, 5)
op <- par(mar=(c(4.8, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
plot(c(0,850), c(500,2000), type="n", lwd=2, xlab="time (sec)", ylab=expression(volume~microns^3))
for(id in c(120) ) {
  cur.traj <- traj[traj$trajectory.id == id,]
  cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
  cur.traj$time <- cur.traj$time - 850
##  cur.traj <- cur.traj[0 <= cur.traj$time & cur.traj$time <= 800,]
  points(loess.smooth(cur.traj$time, runmed(cur.traj$volume, 4), family="gaussian"), type="l", lwd=3, col="darkcyan")
  points(cur.traj$time, cur.traj$volume, type="p", pch=21, bg="gray", col="black", cex=1.2)
}
par(op)
dev.off()


