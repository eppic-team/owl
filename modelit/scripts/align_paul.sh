#!/bin/sh
# Convenience script for structurally aligning templates using paul
# Henning Stehr, h.stehr@molgen.mpg.de
# 11/Apr/2008

if [ -z "$3" ]
then
  echo "Usage: $0 <list.templates> <outfile.fasta> <tempdir> <veryfast|fast|slow|veryslow>"
  exit 1
fi

# input
templates=$1
outfile=$2
tempdir=$3
paul_mode=$4

# parameters
ct="Cb"
d="8.0"

# directories
scriptsdir=`dirname $0`

# output
basename=`basename $templates`

paul_listfile="$tempdir/$basename.paul"
paul_seqfile="$tempdir/$basename.paul.seqs"

# if only one template is given, we skip the alignment steps and simply output the sequence of the template
numTemps=`cat $templates | grep -v "^#" | wc -l`
if [ "$numTemps" -eq "0" ]
then
echo "NO TEMPLATES FOUND, CAN'T ALIGN!"
exit 1
fi
if [ "$numTemps" -eq "1" ]
then
	echo "SINGLE TEMPLATE, SKIPPING ALIGNMENT STEP!"
	for template in `cat $templates | grep -v "^#"`	# this loop should stop after 1 iteration unless something is wrong
	do
		pdb=`echo $template | cut -c1-4`
		chain=`echo $template | cut -c5`
		dumpseq -p ${pdb}${chain} -s > $outfile
		echo "Template sequence written to $outfile"
		exit 0
	done
fi

# make sure files are empty
rm -f $paul_listfile
rm -f $paul_seqfile

echo "Getting graphs"
# parse template list
for template in `cat $templates | grep -v "^#"`
do
  echo $template
  pdb=`echo $template | cut -c1-4`
  chain=`echo $template | cut -c5`
  
  # generate graph directly in paul format (with -P)
  # we capture the output (number of successful structures) so that we can know if there were problems getting the info from pdbase
  numPassed=`genGraph -p ${pdb}${chain} -d $d -t $ct -o $tempdir -P | grep "^Number" | sed -e "s/Number of structures done successfully: \(.*\)/\1/"`
  graphfile="${tempdir}/${pdb}${chain}_${ct}_${d}.cm"

  if [ "$numPassed" -ne "1" ]
  then
      echo "No structure info could be found for $pdb, skipping"
  else
      dumpseq -p ${pdb}${chain} -s > $graphfile.seq
      # add to paul list
      echo `basename $graphfile` >> $paul_listfile
      echo `basename $graphfile.seq` >> $paul_seqfile
  fi
done

# run paul
cur=`pwd`
cd $tempdir
paul -i `basename $paul_listfile` -s `basename $paul_seqfile` -m $paul_mode
if [ "$?" -ne "0" ]
then
    # paul fails sometimes: better to catch it soon and exit
    echo "ERROR: paul failed, exit code was $?"
    exit 1
fi
cd $cur

# convert result to fasta format
t_coffee -other_pg seq_reformat -in ${paul_listfile}.prog.aln -output fasta_aln > ${outfile}.temp

# convert headers and upper casing the sequences
cat ${outfile}.temp | sed 's/_Cb_8.0//' | sed 's/_//' | perl -ne 'if(/^>/){print}else{print uc($_)}' > $outfile

# rename/delete temporary files
rm -f ${outfile}.temp
for template in `cat $templates`
do
  pdb=`echo $template | cut -c1-4`
  chain=`echo $template | cut -c5`  
  graphfile="${tempdir}/${pdb}${chain}_${ct}_${d}.cm"
  rm -f $graphfile $graphfile.seq
done
rm -f $paul_listfile $paul_seqfile $tempdir/$basename.paul.*


echo "Alignment of templates written to $outfile"
