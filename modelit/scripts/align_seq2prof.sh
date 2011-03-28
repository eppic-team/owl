#!/bin/sh
# Convenience script for aligning a sequence against a profile using T-Coffee
# Henning Stehr, h.stehr@molgen.mpg.de
# 10/Apr/2008

if [ -z "$3" ]
then
  echo "Usage: $0 <seq.fasta> <alignm.fasta> <outfile.fasta> <logfile>"
  exit 1
fi

seq=$1
prof=$2
outfile=$3
logfile=$4

scriptsdir=/project/StruPPi/CASP8/scripts
psiblastFile=`dirname $prof`/`basename $prof .templates_aln.fasta`.pdb-psiblast.classic.out

numTempl=`grep -c "^>" $prof`
if [ "$numTempl" -eq "1" ]
then
	echo "Single template, using alignment from psiblast from file $psiblastFile"
	if [ -e "$psiblastFile" ]
	then 
		pdbId=`cat $prof | grep "^>" | sed "s/>//"`
		pdbCode=`echo $pdbId | cut -c1-4`
		pdbChainCode=`echo $pdbId | cut -c5`
		perl $scriptsdir/parse_blast.pl -i $seq -a $psiblastFile -p $pdbCode -c $pdbChainCode > $outfile
		if [ "$?" -ne "0" ]
		then
			echo "ERROR in parsing the psi-blast output, no target to template alignment!"
			exit 1
		fi
	else
		echo "ERROR: couldn't find $psiblastFile, no target to template alignment!"
		exit 1
	fi
	exit 0
fi 


# t-coffee seems to have a bug in outputting directly to fasta_aln format: it will add 
# one extra character to all sequences except for the target. This is why we first write 
# to clustalw format and convert to fasta_aln

# yet another bug in t-coffee. This one is an annoyance rather than a bug:
# it takes the file name of the profile file as the sequence name, but maximum sequence name length is limited to 100, 
# so it fails if path of the file name is really long. To avoid this we first copy the $prof to /tmp

cp $prof /tmp
tmpprof=/tmp/`basename $prof`

tmpoutfile=/tmp/`basename $outfile`.tmp
tmpoutfile2=/tmp/`basename $outfile`.tmp2

t_coffee $seq -profile $tmpprof -profile_comparison=full50 -output=clustalw -outfile=$tmpoutfile -quiet $logfile

t_coffee -other_pg seq_reformat -in $tmpoutfile -output fasta_aln > $tmpoutfile2

cat $tmpoutfile2 | perl -ne 'if(/^>/){print}else{print uc($_)}' > $outfile

rm -f $tmpoutfile $tmpoutfile2 $tmpprof
