#!/bin/sh
# Visualize a similarity matrix written by aglappe.BlastUtils.writeClusterGraph
# Henning Stehr, h.stehr@molgen.mpg.de
# 10/Apr/2008

if [ -z "$2" ]
then
  echo "Usage: $0 <in.matrix> <out.ps>"
  exit 1
fi

infile=$1
outfile=$2

cmd="x <- read.table(\"$infile\"); library(lattice); levelplot(as.matrix(x));plot(hclust(as.dist(x)))#;heatmap(as.matrix(x))"
echo $cmd | R -q --vanilla
mv -f Rplots.ps $outfile
