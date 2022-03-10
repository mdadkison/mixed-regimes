require(ggplot2)

fnames <- c(
"mixedResults/RPreciseOverhSameHistartSims.txt",
"mixedResults/RPreciseOverhSameLostartSims.txt",
"mixedResults/RPreciseOverhDownHistartSims.txt",
"mixedResults/RPreciseNormhSameHistartSims.txt",
"mixedResults/RPreciseNormhSameLostartSims.txt",
"mixedResults/RPreciseNormhDownHistartSims.txt",
"mixedResults/RPreciseUnderhSameHistartSims.txt",
"mixedResults/RPreciseUnderhSameLostartSims.txt",
"mixedResults/RPreciseUnderhDownHistartSims.txt",
"mixedResults/RImpreciseOverhSameHistartSims.txt",
"mixedResults/RImpreciseOverhSameLostartSims.txt",
"mixedResults/RImpreciseOverhDownHistartSims.txt",
"mixedResults/RImpreciseNormhSameHistartSims.txt",
"mixedResults/RImpreciseNormhSameLostartSims.txt",
"mixedResults/RImpreciseNormhDownHistartSims.txt",
"mixedResults/RImpreciseUnderhSameHistartSims.txt",
"mixedResults/RImpreciseUnderhSameLostartSims.txt",
"mixedResults/RImpreciseUnderhDownHistartSims.txt")

######## have to manually increment because dev.off too slow to keep up with loop? ######################
i <- 18 
xfile <- read.table(fnames[i],header=TRUE)
df <- data.frame(Smsy=ifelse(xfile$Smsy>200000,200000,xfile$Smsy))
pp <- ggplot(df,) + geom_histogram(aes(Smsy),breaks=seq(from=200000/30,to=200000,length.out=30),color="dark red", fill="dark red")
pp <- pp + ylim(0,500)+ theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                           panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                                           text = element_text(size=20))
pp
outname <- paste(c("mixedHist/",substr(fnames[i], 14, nchar(fnames[i])-8),"_Smsy.jpg"),collapse="")
dev.copy(jpeg,outname)
dev.off()

hcount <- rep(0,30)
kinc <- 200000/30
breaks=seq(from=10000,to=200000,length.out=30)
for (j in 1:1000) {
  k <- trunc(df$Smsy[j]/kinc)
  hcount[k] <- hcount[k]+1
}
hcount
sum(hcount)
############### end of manual step through section #################################

ssort <- sort(df$Smsy)
