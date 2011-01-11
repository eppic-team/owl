## plotting CASP7 benchmarks


### function to plot CASP7 benchmark data
# Parameters:
# mat: data.frame, 1st column target and 2nd method (bla, pb3, gtg), rest of columns the data series.
# title: main title of the plot
# seriestitles: the title for each of the series corresponding to the data columns of mat that will appear in the legend 
# datacols: to use only some selected data columns from mat (otherwise every data column from 3 to the end will be taken)
# NOTE: recipe for rotated x-axis labels taken from:
# http://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-I-create-rotated-axis-labels_003f

plotCASP7benchmark<- function (mat, title, seriestitles, datacols=c(3:ncol(mat)))
{
	
	if (length(datacols)!=length(seriestitles)) {
	   stop("Number of datacols and seriestitles length differ! exiting")
	}
	
	# increase margins to make room for rotated labels
	par(mar=c(7,4,4,2) +0.1) 
	
	# creating the series vector (to be used for colors and pchs)
	# as we have 2 columns before the data columns start, we want to map 3 to 1, 4 to 2 and so on
	# we don't want to use directly the indices of the data cols: we rather go for lower colors 
	# and lower pchs to have simple looking things like circles, triangles and so on
	series<-c()
	for (col in datacols) {
		i<-col-2 
	  	series<-c(series,i)
	}
	
	# plotting matrix: first 2 columns are target and method, rest are data
	matplot(mat[,datacols], ylim=c(0,100),axes=FALSE, type="b", lty=1,ylab="gdt_ts",pch=series,col=series, main=title)
	
	# no labels, we put them after with text
	axis(1,at=1:nrow(mat), labels=FALSE) 
	#axis(1,at=1:nrow(mat),labels=mat[,1]) # for normal labels
	
	# these are the rotated labels
	text(1:nrow(mat), par("usr")[3] - 2.00, srt = 90, adj = 1,labels=mat[,1], xpd=TRUE)
	
	# and the y axis labels
	axis(2,at=seq(0,100,10), labels=seq(0,100,10))
	
	
	startypos<-90
	ypos<-startypos
	i<-1
	for (ser in series) {
	  legend(x=nrow(mat)-10,y=ypos, seriestitles[i], lty=1, bty="n",  pch=ser,col=ser)
	  ypos<-ypos-5
	  i<-i+1
	}

}


### PLOTS


postscript("plots.ps", horizontal=TRUE)

## nopp vs pp20 vs pp30 vs pp20mm0
# blast
bla<-read.table("nopp_pp20_pp30_pp20mm0_bla")
title<-c("CASP7 targets: nopp vs pp20 vs pp30 vs pp20mm0 (blast)")
series<-c("nopp","pp20","pp30","pp20mm0")
plotCASP7benchmark(bla,title,series)

title<-c("CASP7 targets: pp20 vs pp20mm0 (blast)")
series<-c("pp20", "pp20mm0")
plotCASP7benchmark(bla,title,series,datacols=c(4,6))

title<-c("CASP7 targets: pp20 vs pp30 (blast)")
series<-c("pp20", "pp30")
plotCASP7benchmark(bla,title,series,datacols=c(4,5))

#psi blast
pb3<-read.table("nopp_pp20_pp30_pp20mm0_pb3")
title<-c("CASP7 targets: nopp vs pp20 vs pp30 vs pp20mm0 (psi-blast)")
series<-c("nopp","pp20","pp30","pp20mm0")
plotCASP7benchmark(pb3,title,series)

title<-c("CASP7 targets: pp20 vs pp20mm0 (psi-blast)")
series<-c("pp20", "pp20mm0")
plotCASP7benchmark(pb3,title,series,datacols=c(4,6))

title<-c("CASP7 targets: pp20 vs pp30 (psi-blast)")
series<-c("pp20", "pp30")
plotCASP7benchmark(pb3,title,series,datacols=c(4,5))

#gtg
gtg<-read.table("nopp_pp20_pp30_pp20mm0_gtg")
title<-c("CASP7 targets: nopp vs pp20 vs pp30 vs pp20mm0 (GTG)")
series<-c("nopp","pp20","pp30","pp20mm0")
plotCASP7benchmark(gtg,title,series)

title<-c("CASP7 targets: pp20 vs pp20mm0 (GTG)")
series<-c("pp20", "pp20mm0")
plotCASP7benchmark(gtg,title,series,datacols=c(4,6))

title<-c("CASP7 targets: pp20 vs pp30 (GTG)")
series<-c("pp20", "pp30")
plotCASP7benchmark(gtg,title,series,datacols=c(4,5))

## blast vs psiblast vs gtg
# bla vs pb3 vs gtg for pp20, max10, cct4
all<-read.table("bla_vs_pb3_vs_gtg_pp20_max10_cct4")
title<-c("CASP7 targets: blast vs psi-blast vs gtg \n (pp20, max 10, cct 4)")
series<-c("blast","psi-blast","GTG")
plotCASP7benchmark(all,title,series)

# pb3 vs gtg for pp20, max10, cct4
all<-read.table("pb3_vs_gtg_pp20_max10_cct4")
title<-c("CASP7 targets: psi-blast vs gtg \n (pp20, max 10, cct 4)")
series<-c("psi-blast","GTG")
plotCASP7benchmark(all,title,series)


## nopp vs pp20mm0 vs pp20mm2_om178_5
#blast
all<-read.table("nopp_pp20mm0_pp20mm2om_bla")
title<-c("CASP7 targets: nopp vs pp20mm0 vs pp20mm2_om178_5 \n(blast)")
series<-c("nopp","pp20","pp20mm0","pp20mm2_om178_5")
plotCASP7benchmark(all,title,series)
#psi-blast
all<-read.table("nopp_pp20mm0_pp20mm2om_pb3")
title<-c("CASP7 targets: nopp vs pp20mm0 vs pp20mm2_om178_5 \n(psi-blast)")
series<-c("nopp","pp20","pp20mm0","pp20mm2_om178_5")
plotCASP7benchmark(all,title,series)
#GTG
all<-read.table("nopp_pp20mm0_pp20mm2om_gtg")
title<-c("CASP7 targets: nopp vs pp20mm0 vs pp20mm2_om178_5 \n(GTG)")
series<-c("nopp","pp20","pp20mm0","pp20mm2_om178_5")
plotCASP7benchmark(all,title,series)

dev.off()
