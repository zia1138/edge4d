library(plyr) 
contact <- read.delim("live_data_output/072112_02_t11_50_contact_analysis.txt", stringsAsFactors=F)
contact$time <- contact$time - 850 ## time0 is 850 
contact <- contact[which(contact$neighbor.trajectory.id >= 0),]
contact$pcent <- contact$neighbor.contact.sa / contact$total.contact.sa


PF.basal.cells.data1 <- c(0, 32, 38, 41, 47, 70, 79, 112, 117, 120, 144, 124, 152, 156, 183, 195, 202, 203)
PF.1.neighbors.data1 <- c(25, 34, 57, 74, 75, 87, 89, 98, 108, 109, 114, 116, 124, 146, 148, 177, 181, 198, 206, 221, 238)
PF.2.neighbors.data1 <- c(2, 10, 20, 30, 60, 64, 67, 86, 94, 100, 125, 136, 155, 162, 175, 182, 194, 215, 222, 315)

AF.basal.cells.data1 <- c(137,110,35,40,36,3,140,141,142,178,26)
AF.1.neighbors.data1 <- c(6, 9, 23, 27, 42, 48, 65, 71, 73, 95, 105, 113, 118, 123, 130, 173, 190, 193, 196, 213)
AF.2.neighbors.data1 <- c(4, 18, 58, 68, 72, 82, 84, 85, 88, 106, 131, 133, 145, 150, 153, 157, 159, 207, 214, 216, 226, 236, 247, 271, 279, 317)

all <- NULL
for(traj.id in c(PF.basal.cells.data1, AF.basal.cells.data1, PF.1.neighbors.data1, AF.1.neighbors.data1, PF.2.neighbors.data1, AF.2.neighbors.data1)) {
  current <- contact[which(contact$trajectory.id == traj.id),]
  proc.traj <- function(traj) {
    if( sum(traj$time == 0) == 0 || sum(traj$time == 800) == 0) {
      c(start.sa=NA,end.sa=NA, dx=NA, dy=NA, cx=NA, cy=NA, nx=NA, ny=NA)
    }
    else {
      start.idx <- which(-100 <= traj$time & traj$time <= 0)
      end.idx <- which(700 <= traj$time & traj$time <= 800)

      cx <- traj$cent.x[which(traj$time==0)]
      cy <- traj$cent.y[which(traj$time==0)]
      nx <- traj$neighbor.cent.x[which(traj$time==0)]
      ny <- traj$neighbor.cent.y[which(traj$time==0)]
      dx <- nx - cx
      dy <- ny - cy
      
      ##start.sa <- median(traj$neighbor.surface.area[start.idx], na.rm=T)
      ##end.sa <- median(traj$neighbor.surface.area[end.idx], na.rm=T)
      start.sa <- median(traj$pcent[start.idx], na.rm=T)
      end.sa <- median(traj$pcent[end.idx], na.rm=T)

      c(start.sa = start.sa, end.sa = end.sa, dx=dx, dy=dy, cx=cx, cy=cy, nx=nx, ny=ny)
    }
  }
  sa.frame <- ddply(current, c("trajectory.id", "neighbor.trajectory.id"), proc.traj)
  sa.frame <- sa.frame[which(!(is.na(sa.frame$start.sa) | is.na(sa.frame$end.sa))),]
  sa.frame$incr <- sa.frame$start.sa < sa.frame$end.sa
  all <- rbind(all, sa.frame)
}


pdf("contact_area_dist.pdf")
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
op <- par(mar=(c(4, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
plot(type="n", c(-8,8), c(-8,8), xlab="anterior-posterior span, microns", ylab = "dorsal-ventral span, microns")
points(all$dx[all$incr], all$dy[all$incr], xlim=c(-8,8), ylim=c(-8,8), col="green3")
points(all$dx[!all$incr], all$dy[!all$incr], xlim=c(-8,8), ylim=c(-8,8), col="magenta")
points(0,0, col="black", cex=4, pch=20)
par(op)
dev.off()          


all.incr <- all[which(all$incr),]
all.decr <- all[which(!all$incr),]

pdf("contacts.pdf", 8, 4.25)
## ‘c(bottom, left, top, right)’ ‘c(5, 4, 4, 2) + 0.1’.
op <- par(mar=(c(4, 4.8, 1, 1.25) + 0.1), lwd=2, ps=16)
plot(c(0, 140),c(0,54), type="n", xlab="anterior-posterior span, microns", ylab="dorsal-ventral span, microns")
contact.0 <- contact[which(contact$time == 0),]
segments(contact.0$cent.x, contact.0$cent.y, contact.0$neighbor.cent.x, contact.0$neighbor.cent.y, col="lightgray")

segments(all.incr$cx, all.incr$cy, all.incr$nx, all.incr$ny, col="green", lwd=4)
segments(all.decr$cx, all.decr$cy, all.decr$nx, all.decr$ny, col="magenta", lwd=4)


points(c(contact.0$cent.x, contact.0$neighbor.cent.x),c(contact.0$cent.y, contact.0$neighbor.cent.y),  pch=20, col="lightgray" )
points(c(all$cx, all$nx),c(all$cy, all$ny),  pch=20, col="black" )


focal.x <- unique(contact.0[contact.0$trajectory.id == 141,]$cent.x)
focal.y <- unique(contact.0[contact.0$trajectory.id == 141,]$cent.y)
points(focal.x, focal.y,  pch=22, col="yellow4" , cex=2, lwd=3, bg="yellow")

focal.x <- unique(contact.0[contact.0$trajectory.id == 23,]$cent.x)
focal.y <- unique(contact.0[contact.0$trajectory.id == 23,]$cent.y)
points(focal.x, focal.y,  pch=22, col="darkred" , cex=2, lwd=3, bg="red")

focal.x <- unique(contact.0[contact.0$trajectory.id == 40,]$cent.x)
focal.y <- unique(contact.0[contact.0$trajectory.id == 40,]$cent.y)
points(focal.x, focal.y,  pch=22, col="blue" , cex=2, lwd=3, bg="cyan")

## focal.x <- unique(contact.0[contact.0$trajectory.id == 26,]$cent.x)
## focal.y <- unique(contact.0[contact.0$trajectory.id == 26,]$cent.y)
## points(focal.x, focal.y,  pch=20, col="orange" , cex=2)

par(op)
dev.off()


# 141
# decreasing 26, 40, 178
# increasing 23, 142

pdf("incr_example.pdf", 6,6)
op <- par(mar=(c(4.8, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
plot(c(0,800), c(5,250), type="n", lwd=2, xlab="time (sec)", ylab=expression(cell~cell~contact~surface~area~microns^2))
cur.traj <- contact[which(contact$trajectory.id == 141 & contact$neighbor.trajectory.id == 23),]
cur.traj <- cur.traj[0 <= cur.traj$time &  cur.traj$time <= 800,]
cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
points(loess.smooth(cur.traj$time, runmed(cur.traj$neighbor.contact.sa, 3), family="gaussian"), type="l", lwd=3, col="red")
points(cur.traj$time, cur.traj$neighbor.contact.sa, type="p", pch=21, bg="red", col="darkred", cex=1.2)
par(op)
dev.off()


pdf("decr_example.pdf", 6,6)
op <- par(mar=(c(4.8, 4.8, 1, 1.25) + 0.1), lwd=2, ps=22)
plot(c(0,800), c(5,250), type="n", lwd=2, xlab="time (sec)", ylab=expression(cell~cell~contact~surface~area~microns^2))
cur.traj <- contact[which(contact$trajectory.id == 141 & contact$neighbor.trajectory.id == 40),]
cur.traj <- cur.traj[0 <= cur.traj$time &  cur.traj$time <= 800,]
cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
points(loess.smooth(cur.traj$time, runmed(cur.traj$neighbor.contact.sa, 3), family="gaussian"), type="l", lwd=3, col="cyan")
points(cur.traj$time, cur.traj$neighbor.contact.sa, type="p", pch=21, bg="cyan", col="blue", cex=1.2)
par(op)
dev.off()


## incr.contacts <- all[all$incr, c(1,2)]
## decr.contacts <- all[!all$incr, c(1,2)]

## pdf("incr.pdf")
## for(i in 1:dim(incr.contacts)[1]) {
##   A <- incr.contacts$trajectory.id[i]
##   B <- incr.contacts$neighbor.trajectory.id[i]
##   cur.traj <- contact[which(contact$trajectory.id == A & contact$neighbor.trajectory.id == B),]
##   cur.traj <- cur.traj[0 <= cur.traj$time &  cur.traj$time <= 800,]
##   cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
##   plot(cur.traj$time, cur.traj$neighbor.contact.sa, type="b")
##   title(c(A,B))
## }
## dev.off()

## pdf("decr.pdf")
## for(i in 1:dim(decr.contacts)[1]) {
##   A <- decr.contacts$trajectory.id[i]
##   B <- decr.contacts$neighbor.trajectory.id[i]
##   cur.traj <- contact[which(contact$trajectory.id == A & contact$neighbor.trajectory.id == B),]
##   cur.traj <- cur.traj[0 <= cur.traj$time &  cur.traj$time <= 800,]
##   cur.traj <- cur.traj[sort(cur.traj$time, index.return=T)$ix,]
##   plot(cur.traj$time, cur.traj$neighbor.contact.sa, type="b")
##   title(c(A,B))
## }
## dev.off() 

