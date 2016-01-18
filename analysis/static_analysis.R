early.fixed <- read.delim("fixed_final_pass1/OR_OliGreen_BP106_early03_static_analysis.txt")
mid.fixed <- read.delim("fixed_final_pass1/OR_OliGreen_BP106_mid04_static_analysis.txt")
late.fixed <- read.delim("fixed_final_pass1/OR_OliGreen_BP106_late01_static_analysis.txt")
results <- t(data.frame(##detected = c(early.fixed$detected, mid.fixed$detected, late.fixed$detected ),
                        split = c(early.fixed$split, mid.fixed$split, late.fixed$split),
                        merged = c(early.fixed$merged, mid.fixed$merged, late.fixed$merged),
                        ##                        junk = c(early.fixed$junk, mid.fixed$junk, late.fixed$junk),
                        missed = c(early.fixed$missed,  mid.fixed$missed, late.fixed$missed)))

pdf("nucs_fixed.pdf", 3,5)
barplot(c(early.fixed$nucs, mid.fixed$nucs, late.fixed$nucs), col=c("purple"),
        names.arg = c("early", "mid", "late"), ylab="Reconstructed Nuclei", ylim=c(0,700), legend=F)
dev.off()

pdf("edge2_errors_fixed.pdf", 3,5)
barplot(as.matrix(results), col=c("gray", "red", "yellow"),##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = c("early", "mid", "late"),  ylab="Error Count", ylim=c(0,700), legend=F)
dev.off()



early.fixed <- read.delim("fixed_final_pass1_nm2010/OR_OliGreen_BP106_early03_static_analysis.txt")
mid.fixed <- read.delim("fixed_final_pass1_nm2010/OR_OliGreen_BP106_mid04_static_analysis.txt")
late.fixed <- read.delim("fixed_final_pass1_nm2010/OR_OliGreen_BP106_late01_static_analysis.txt")
results <- t(data.frame(##detected = c(early.fixed$detected, mid.fixed$detected, late.fixed$detected ),
                        split = c(early.fixed$split, mid.fixed$split, late.fixed$split),
                        merged = c(early.fixed$merged, mid.fixed$merged, late.fixed$merged),
                        ##                        junk = c(early.fixed$junk, mid.fixed$junk, late.fixed$junk),
                        missed = c(early.fixed$missed,  mid.fixed$missed, late.fixed$missed)))


pdf("nm2010_errors_fixed.pdf", 3,5)
barplot(as.matrix(results), col=c("gray", "red", "yellow"),##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = c("early", "mid", "late"),  ylab="Error Count", ylim=c(0,700), legend=F)
dev.off()



early.fixed <- read.delim("fixed_final_pass1_acme/OR_OliGreen_BP106_early03_static_analysis.txt")
mid.fixed <- read.delim("fixed_final_pass1_acme/OR_OliGreen_BP106_mid04_static_analysis.txt")
late.fixed <- read.delim("fixed_final_pass1_acme/OR_OliGreen_BP106_late01_static_analysis.txt")
results <- t(data.frame(##detected = c(early.fixed$detected, mid.fixed$detected, late.fixed$detected ),
                        split = c(early.fixed$split, mid.fixed$split, late.fixed$split),
                        merged = c(early.fixed$merged, mid.fixed$merged, late.fixed$merged),
                        ##                        junk = c(early.fixed$junk, mid.fixed$junk, late.fixed$junk),
                        missed = c(early.fixed$missed,  mid.fixed$missed, late.fixed$missed)))


pdf("acme_errors_fixed.pdf", 3,5)
barplot(as.matrix(results), col=c("gray", "red", "yellow"),##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = c("early", "mid", "late"),  ylab="Error Count", ylim=c(0,700), legend=F)
dev.off()



live.data <- read.delim("live_data_output/072112_02_t11_50_static_analysis.txt", stringsAsFactors=F)
live.data <- live.data[live.data$time - 850 <= 950,]
results <- t(data.frame(##detected = live.data$detected, #green
                        split = live.data$split, #gray
                        merged = live.data$merged, # red
##                        junk = live.data$junk, # blue
                        missed = live.data$missed)) # yellow


pdf("EDGE2_cells.pdf",8.5,3.5)
op <- par(ps=14, mar=c(4, 4, 2.5, 2) +  0.1)  # ‘c(bottom, left, top, right)’
barplot(as.matrix(results), col=c("gray", "red", "yellow"), ##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = live.data$time - 850, ylab="Error Count", ylim=c(0,275), legend=F, xlab="Time (seconds)")
par(op)
dev.off()


pdf("nucs_live.pdf",8.5,3.5)
op <- par(ps=14, mar=c(4, 4, 2.5, 2) +  0.1)  # ‘c(bottom, left, top, right)’
barplot(live.data$nucs, col=c("magenta"),
        names.arg = live.data$time - 850, ylab="Reconstructed Nuclei", ylim=c(0,275), legend=F, xlab="Time (seconds)")
par(op)
dev.off()


##live.data <- read.delim("live_data_output/071812_06_t6-37_static_analysis.txt", stringsAsFactors=F)
live.data <- read.delim("live_data_output/072112_06_t7_51_static_analysis.txt", stringsAsFactors=F)
live.data <- live.data[live.data$time - 1100 <= 950,]
results <- t(data.frame(##detected = live.data$detected, #green
                        split = live.data$split, #gray
                        merged = live.data$merged, # red
##                        junk = live.data$junk, # blue
                        missed = live.data$missed)) # yellow


pdf("EDGE2_cells_rep2.pdf",8.5,3.5)
op <- par(ps=14, mar=c(4, 4, 2.5, 2) +  0.1)  # ‘c(bottom, left, top, right)’
barplot(as.matrix(results), col=c("gray", "red", "yellow"), ##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = live.data$time - 1100, ylab="Error Count", ylim=c(0,275), legend=F, xlab="Time (seconds)")
par(op)
dev.off()
pdf("nucs_live_rep2.pdf",8.5,3.5)
op <- par(ps=14, mar=c(4, 4, 2.5, 2) +  0.1)  # ‘c(bottom, left, top, right)’
barplot(live.data$nucs, col=c("magenta"),
        names.arg = live.data$time - 1100, ylab="Reconstructed Nuclei", ylim=c(0,275), legend=F, xlab="Time (seconds)")
par(op)
dev.off()



live.data <- read.delim("live_data_output/072112_02_t11_50_static_analysis_NM2010.txt", stringsAsFactors=F)
live.data <- live.data[live.data$time - 850 <= 950,]
results <- t(data.frame(##detected = live.data$detected, #green
                        split = live.data$split, #gray
                        merged = live.data$merged, # red
                        ##junk = live.data$junk, # blue
                        missed = live.data$missed)) # yellow


pdf("NM2010_cells.pdf",8.5,3.5)
op <- par(ps=14, mar=c(4, 4, 2.5, 2) +  0.1)  # ‘c(bottom, left, top, right)’
barplot(as.matrix(results), col=c("gray", "red", "yellow"), ##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = live.data$time - 850, ylab="Error Count", ylim=c(0,275), legend=F, xlab="Time (seconds)")
par(op)
dev.off()


live.data <- read.delim("live_data_output/072112_02_t11_50_static_analysis_ACME.txt", stringsAsFactors=F)
live.data <- live.data[live.data$time - 850 <= 950,]
results <- t(data.frame(##detected = live.data$detected, #green
                        split = live.data$split, #gray
                        merged = live.data$merged, # red
                        ##junk = live.data$junk, # blue
                        missed = live.data$missed)) # yellow


pdf("ACME_cells.pdf",8.5,3.5)
op <- par(ps=14, mar=c(4, 4, 2.5, 2) +  0.1)  # ‘c(bottom, left, top, right)’
barplot(as.matrix(results), col=c("gray", "red", "yellow"), ##col=c("green", "gray", "red", "blue", "yellow"),
        names.arg = live.data$time - 850, ylab="Error Count", ylim=c(0,275), legend=F, xlab="Time (seconds)")
par(op)
dev.off()



