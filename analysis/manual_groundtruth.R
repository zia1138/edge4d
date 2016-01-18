 #x <- read.delim("fixed_final_pass1/OR_OliGreen_BP106_mid04_manual_groundtruth.txt", stringsAsFactors=F)

plot.gtdata <- function(data.set, dim) {
  filename <- paste("fixed_final_pass1/OR_OliGreen_BP106_", data.set, "_manual_groundtruth.txt", sep="")
  x <- read.delim(filename, stringsAsFactors=F)

  y <- x[x$seg.size > 0 & x$status == 0,]
  CELLS <- grep(paste("c0_",dim,sep=""), y$data.id, fixed=T)
  B <- grep(paste("c1_",dim,sep=""), y$data.id, fixed=T)

#  print(wilcox.test(y$RMS[CELLS], y$RMS[B], alternative="greater"))

  pdf(paste(data.set,"_",dim,".pdf", sep=""), 3.5,5)
  op <- par(ps=18, mar=c(1, 2.1, 1, 2) +  0.1)  # ‘c(bottom, left, top, right)’
  boxplot(y$RMS[B], y$RMS[CELLS], ylim=c(0,1.25), col=c("magenta","green"), lwd=2)
  par(op)
##  print(paste("N(cell)=", length(A) , " N(nuclei)=", length(B), sep=""))
  dev.off()
}

plot.gtdata("early03", "xy")
plot.gtdata("early03", "xz")
plot.gtdata("early03", "yz")

plot.gtdata("mid04", "xy")
plot.gtdata("mid04", "xz")
plot.gtdata("mid04", "yz")


plot.gtdata("late01", "xy")
plot.gtdata("late01", "xz")
plot.gtdata("late01", "yz")
