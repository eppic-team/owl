#!/bin/sh
# Script to run the whole CASP prediction pipeline
# Cluster ready: to submit it with qsub, just make sure you use absolute paths


## CONSTANTS
bindir=/project/StruPPi/CASP8/scripts
tmpdir=/tmp

DEFAULT_PDBBLASTDB=seqs_pdbase_20080903.fix.reps.fa
DEFAULT_MAX_HITS=20
DEFAULT_PSIBLAST_ITER=5

paulSpeed=slow

cts="Ca,Cg"
dists="8.0,8.0"
numDistgeomModels=40
DEFAULT_CCT="0.4"

DEFAULT_GTG_DIR="/project/StruPPi/CASP8/gtg"

CASP_AUTHOR="3757_3251_6331"

CASP_METHOD_FILE="/project/StruPPi/CASP8/scripts/method.txt"

## OPTIONS 

targetseqfile=""
outdir=""
pdbblastdb="$DEFAULT_PDBBLASTDB"
psiblastopt=""
psiblastIter="$DEFAULT_PSIBLAST_ITER"
templatefile=""
skiptinkeropt=""
CCT="$DEFAULT_CCT"
selectTemplates=""
gtgDir=$DEFAULT_GTG_DIR
maxHits=$DEFAULT_MAX_HITS
noPhiPsi=""
noTransOmega=""

while getopts i:o:Kj:DP:G:t:s:l:x:FO opt
do
  case "$opt" in
    i) targetseqfile="$OPTARG";;
    o) outdir="$OPTARG";;
    K) psiblastopt="-K";;
    j) psiblastIter="$OPTARG";;
    D) skiptinkeropt="-D";;
    P) pdbblastdb="$OPTARG";;
    G) gtgDir="$OPTARG";;
    t) templatefile="$OPTARG";;
    s) CCT="$OPTARG";;
    l) selectTemplates="$OPTARG";;
    x) maxHits="$OPTARG";;
    F) noPhiPsi="-F";;
    O) noTransOmega="-O";;
  esac
done


if [ -z "$targetseqfile" ] || [ -z "$outdir" ]
then
    echo "Usage: $0 "
    echo " -i <file>       : the target sequence file"
    echo " -o <dir>        : the output directory"
    echo " [-K]            : no psiblast, just blast"
    echo " [-j <number>]   : number of psiblast iterations, default: $DEFAULT_PSIBLAST_ITER"
    echo " [-D]            : no Tinker reconstruction, just output average graph"
    echo " [-P <blast db>] : PDB blast db, default: $DEFAULT_PDBBLASTDB"
    echo " [-G <dir>]      : directory containing GTG output, default: $DEFAULT_GTG_DIR" 
    echo " [-t <file>]     : specifying a templates list file, no template selection step is performed (skip STEP 1)"
    echo " [-s <number>]   : CCT (contact conservation threshold) for the graph averaging step, default: $DEFAULT_CCT"
    echo " [-l <B, P or G>]: forces selection of (B)LAST template list, (P)SIBLAST template list or (G)TG template list (even if they are empty!)"
    echo " [-x <number>]   : a hard maximum for the number of hits to be taken in the template selection step, default: $DEFAULT_MAX_HITS"
    echo " [-F]            : don't use phi/psi constraints for reconstruction. Default: phi/psi constraints used"
    echo " [-O]            : don't enforce trans conformation for omega torsion angle (peptide bond). Default: enforcing"
    echo ""
    echo "REMEMBER: always use absolute paths to run this script in the cluster"
    echo ""
    exit 1

fi





## file definitions
basename=`basename $targetseqfile .fa`
templatesAlnFile=$outdir/$basename.templates_aln.fasta
target2templAlnFile=$outdir/$basename.target2templ.fasta
tcoffeelog=$outdir/$basename.tcoffee.log
ssCompareFile=$outdir/$basename.ss.compare
psipredFile=$outdir/$basename.horiz #input file!

### STEP 1: template selection
if [ -z "$templatefile" ]
then
    echo "################################################"
    echo "# STEP 1:  TEMPLATE SELECTION"
    echo "################################################"
    
    templatefile=$outdir/$basename.templates
    
    selectTemplatesOpt=""
    if [ -n "$selectTemplates" ]
    then
    	selectTemplatesOpt="-l $selectTemplates"
	fi
    $bindir/templateSelection -i $targetseqfile -o $outdir -P $pdbblastdb -b $basename -G $gtgDir $psiblastopt $selectTemplatesOpt -x $maxHits -j $psiblastIter

    numTemps=`cat $templatefile | wc -l`
    if [ "$numTemps" -eq "0" ]
    then
	echo "NO TEMPLATES FOUND, CAN'T MODEL!"
	exit 1
    fi
fi

### STEP 2: structural alignment of templates
echo "################################################"
echo "# STEP 2:  STRUCTURAL ALIGNMENT OF TEMPLATES"
echo "################################################"
$bindir/align_paul.sh $templatefile $templatesAlnFile $tmpdir $paulSpeed
if [ "$?" -ne "0" ]
then
  echo "ERROR: paul alignment failed. Exiting"
  exit 1
fi


### STEP 3: target sequence to templates profile alignment
echo "################################################"
echo "# STEP 3:  TARGET TO TEMPLATES ALIGNMENT"
echo "################################################"
$bindir/align_seq2prof.sh $targetseqfile $templatesAlnFile $target2templAlnFile $tcoffeelog
echo "Wrote target to template alignment to $target2templAlnFile"
echo "Wrote t-coffee output to $tcoffeelog"
# generate secondary structure comparison, needs psipred output file, searching in output dir
if [ -e $psipredFile ]
then
  $bindir/displaySSMatching -n -a $target2templAlnFile -s $psipredFile > $ssCompareFile
  echo "Wrote secondary structure comparison to $ssCompareFile"
else
  echo "Could not find psipred output, please run displaySSMatching manually to view secondary structure comparison."
fi

### STEP 4: graph average and reconstruction
echo "################################################"
echo "# STEP 4:  GRAPH AVERAGING AND RECONSTRUCTION"
echo "################################################"


if [ -z "$skiptinkeropt" ]
then
  $bindir/averageGraph -f $targetseqfile -P $templatefile -t $cts -d $dists -b $basename -o $outdir -a $target2templAlnFile -s $CCT \
  -r $numDistgeomModels -c $CASP_AUTHOR -m $CASP_METHOD_FILE $noPhiPsi $noTransOmega
else
  $bindir/averageGraph -f $targetseqfile -P $templatefile -t $cts -d $dists -b $basename -o $outdir -a $target2templAlnFile -s $CCT  
  echo "Skipping reconstruction step. No 3D model will be written."
fi
