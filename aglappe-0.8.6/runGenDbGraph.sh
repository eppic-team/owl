#!/bin/sh

if [ -z "$4" ]
then
	echo "usage: $0 <listfile> <cutoff> <ct> <outdb>"
        echo "where listfile contains two columns, pdb code and chain code"
	exit 1
fi

listfile=$1
cutoff=$2
ct=$3
outdb=$4


exec 0<$listfile

while read pdbCode pdbChainCode
do
  java -Xmx1024m genDbGraph -p $pdbCode -c $pdbChainCode -d $cutoff -t $ct -o $outdb
done

