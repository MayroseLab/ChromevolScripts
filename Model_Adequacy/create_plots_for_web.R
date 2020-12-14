args = commandArgs(trailingOnly=TRUE)
setwd(args[1])

stats_ind = c(1,2,3,4)
stats_names = c("Variance","Entropy","Parsimony score","Parsimony score vs. time slope")

# read distributions
dists = read.table("stats_dist_sims",sep=c(","),stringsAsFactors = F)
dists$V1 = as.numeric(gsub( "\\[", "", as.character(dists$V1)) )
dists[,ncol(dists)] = as.numeric(gsub( "\\]", "", as.character(dists[,ncol(dists)])))
dists = as.matrix(dists)

# read original stats
orig_stats = read.table("orig_stats",sep=",",stringsAsFactors = F)
orig_stats = orig_stats[stats_ind]
names(orig_stats) = stats_names

# handle missing statistics
none_ind = which(orig_stats=="None")
if (length(none_ind) == length(stats_ind)){ # 
	pdf("dists.pdf")
	hist(0, main = "No statistics to show", xlab="", ylab="")
	dev.off()
	quit()
}
if (length(none_ind) > 0){ 
	stats_ind = setdiff(stats_ind,none_ind)
	stats_names = stats_names[stats_ind]
}

# read percentiles
percentiles = read.table("percentiles_limits", stringsAsFactors = F, sep=",")

# read true percentiles
true_p = read.table("true_percentiles", sep=",", stringsAsFactors = F)
true_p = true_p[stats_ind]
names(true_p)  = stats_names


pdf("dists.pdf")
par(mfrow = c(2,2))
for (i in stats_ind){
  y = hist(dists[i,], main = paste(stats_names[i]), ylab = "", xlab = "")
  mtext(paste("95% confidence bounds: ","[",percentiles[i,1],", ",percentiles[i,2],"]","\n",
        "True statistic: ",orig_stats[i],"  Percentile: ",true_p[1,i],sep=""), line = -0.5, cex = 0.8)
  segments(x0=as.numeric(orig_stats[i]),y0=0,x1=as.numeric(orig_stats[i]),y1=max(y$counts)*0.9,col="red",lwd=3)
  segments(x0=as.numeric(percentiles[i,1]),y0=0,x1=as.numeric(percentiles[i,1]),y1=max(y$counts)*0.9,col="blue",lty=2)
  segments(x0=as.numeric(percentiles[i,2]),y0=0,x1=as.numeric(percentiles[i,2]),y1=max(y$counts)*0.9,col="blue",lty=2)
}
dev.off()

#### tails of distribution