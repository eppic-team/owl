

DECOYSETS<-c("4state_reduced", "fisa", "fisa_casp3", "hg_structal", "ig_structal", "ig_structal_hires", "lattice_ssfit")


for (set in DECOYSETS) {
	files<-list.files(path = ".", pattern = paste(set,"_.*.scores",sep=""))
	
	cols=0	
	if (length(files)%%2==0) {
		cols=length(files)/2
	} else {
		cols=(length(files)+1)/2
	}
	par(mfrow = c(cols, cols))
	
	for (file in files) {
		data<-read.table(file)
		plot(data$V2,data$V3,xlab="score",ylab="rmsd",main=file)
	}

}