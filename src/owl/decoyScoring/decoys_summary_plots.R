
DECOYSETS<-c("4state_reduced", "fisa", "fisa_casp3", "hg_structal", "ig_structal", "ig_structal_hires", "lattice_ssfit", "lmds", "vhp_mcmd","means")
TYPES<-c("atomtype","atomcount","atomcomb","restype","rescount","rescomb")
COLORS<-c(2,5,3,6,4,8)

postscript("plots-means.ps")

for (decoySet in DECOYSETS) {
	data<-read.table(paste(decoySet,".summary",sep=""),header=FALSE)
	data<-data[order(data$V3),] # sorting by cutoff ascending
	
	# zscores
	xlim<-c(min(data$V3),max(data$V3))
	#ylim<-c(min(data$V7),max(data$V7))
	ylim<-c(-2,6.5)
	yspacing<-(ylim[2]-ylim[1])/25
	i=1
	for (type in TYPES) {
		if (i==1) {
			plot(data$V3[data$V2==type],data$V7[data$V2==type],xlim=xlim,ylim=ylim,
					xlab="cutoff",ylab="zscore",main=paste("Mean native z-scores for decoy set '",decoySet,"' (size ",data$V5[1],")",sep=""),
					type="b",col=COLORS[i])
		} else {
			points(data$V3[data$V2==type],data$V7[data$V2==type],xlim=xlim,ylim=ylim,type="b",col=COLORS[i])
		}
		legend(x=xlim[2]-3,y=ylim[2]-yspacing*i, type, col=COLORS[i], lty=1 , bty="n", pch=1)
		i=i+1
		
	}
	abline(v=6,col=8,lty=2)
	abline(h=seq(-2,6,2),col=8)
	
	# proportion native ranked 1
	ylim<-c(0,1)
	yspacing<-(ylim[2]-ylim[1])/25
	i=1
	for (type in TYPES) {
		if (i==1) {
			plot(data$V3[data$V2==type],data$V6[data$V2==type],xlim=xlim,ylim=ylim,
					xlab="cutoff",ylab="native ranked 1st",main=paste("Proportion of ranked-1st natives for decoy set '",decoySet,"' (size ",data$V5[1],")",sep=""),
					type="b",col=COLORS[i])
		} else {
			points(data$V3[data$V2==type],data$V6[data$V2==type],xlim=xlim,ylim=ylim,type="b",col=COLORS[i])
		}
		legend(x=xlim[2]-3,y=ylim[2]-yspacing*i, type, col=COLORS[i], lty=1 , bty="n", pch=1)
		i=i+1
		
	}
	abline(v=6,col=8,lty=2)
	
	# rank correlation
	ylim<-c(0,1)
	yspacing<-(ylim[2]-ylim[1])/25
	i=1
	for (type in TYPES) {
		if (i==1) {
			plot(data$V3[data$V2==type],-data$V8[data$V2==type],xlim=xlim,ylim=ylim,
					xlab="cutoff",ylab="spearman",main=paste("Mean rank correlation for decoy set '",decoySet,"' (size ",data$V5[1],")",sep=""),
					type="b",col=COLORS[i])
		} else {
			points(data$V3[data$V2==type],-data$V8[data$V2==type],xlim=xlim,ylim=ylim,type="b",col=COLORS[i])
		}
		legend(x=xlim[2]-3,y=ylim[2]-yspacing*i, type, col=COLORS[i], lty=1 , bty="n", pch=1)
		i=i+1
		
	}
	abline(v=6,col=8,lty=2)
	
}

dev.off()


DECOYSETS<-c("4state_reduced", "fisa", "fisa_casp3", "hg_structal", "ig_structal", "ig_structal_hires", "lattice_ssfit", "lmds", "vhp_mcmd")
TYPES<-c("atomcomb","atomcount","atomtype","rescomb","rescount","restype")
COLORS<-c(2,2,2,4,4,4)
PCHS<-c(1,4,5,1,4,5)

postscript("plots-individual.ps")
all.atType.scores<-c()
all.atCount.scores<-c()
all.resType.scores<-c()
all.resCount.scores<-c()

for (decoySet in DECOYSETS) {
	data<-read.table(paste(decoySet,".indiv.summary",sep=""),header=FALSE)
	# zscores
	xlim<-c(1,dim(data)[1])
	#ylim<-c(min(data$V6),max(data$V6))
	ylim<-c(-2,6.5)
	yspacing<-(ylim[2]-ylim[1])/25
	i=1
	for (type in TYPES) {
		if (i==1) {
			plot(data[,3*i],ylim=ylim,
					xlab="decoys",ylab="zscore",main=paste("z-scores for decoy set '",decoySet,"'",sep=""),
					type="p",pch=PCHS[i],col=COLORS[i],
					axes=FALSE)
			axis(2)
			axis(1, at=c(1:dim(data)[1]), labels=data$V1)
			box()
		} else {
			points(data[,3*i],ylim=ylim,type="p",pch=PCHS[i],col=COLORS[i])
		}
		legend(x=xlim[2]-(xlim[2]-xlim[1])/10,y=ylim[2]-yspacing*(i-1), type, col=COLORS[i], lty=0 , bty="n", pch=PCHS[i])
		i=i+1
		
	}
	
	abline(v=c(1:dim(data)[1]),col=8,lty=3)
	
	all.atCount.scores<-c(all.atCount.scores,data$V3)
	all.atType.scores<-c(all.atType.scores,data$V6)
	all.resCount.scores<-c(all.resCount.scores,data$V9)
	all.resType.scores<-c(all.resType.scores,data$V12)
	
	
}

# correlation of z-scores atom type vs atom count
corr<-cor(all.atCount.scores,all.atType.scores,method="spearman")
plot(all.atCount.scores,all.atType.scores,xlab="z-score atom count",ylab="z-score atom type",
		main="correlation of z-scores atom type vs count")
legend(x=max(all.atCount.scores)-2,y=min(all.atType.scores)+1,corr,col=3,bty="n")
# correlation of z-scores res type vs res count
corr<-cor(all.resCount.scores,all.resType.scores,method="spearman")
plot(all.resCount.scores,all.resType.scores,xlab="z-score res count",ylab="z-score res type",
		main="correlation of z-scores res type vs count")
legend(x=max(all.resCount.scores)-2,y=min(all.resType.scores)+1,corr,col=2,bty="n")
dev.off()
