#!/bin/sh 
# Script to continue an MD simulation given a gro file output of an earlier MD simulation.

# CONSTANTS

MDPDIR=/project/StruPPi/Software/gromacs-3.3.3/mdp_files
# MD simulation parameters
MDP_MD=$MDPDIR/md.mdp

SOLVENTGROUP=12

DEFAULT_NP=1

# END CONSTANTS

# command line parameters
inGro=""
outDir=""
np=$DEFAULT_NP
endT=""

while getopts i:o:p:e: opt
do
  case "$opt" in
    i) inGro="$OPTARG";;
    o) outDir="$OPTARG";;
    p) np="$OPTARG";;
	e) endT="$OPTARG";;
  esac
done

if [ -z "$inGro" ] || [ -z "$outDir" ] || [ -z "$endT" ]
then
	echo ""
    echo "Usage: $0 "
    echo "  -i  <input gro file>"
    echo "  -o  <out dir> "
    echo "  -e  <end time (ns)>"    
    echo " [-p] <number of processes>"
    echo "Use end time (ns) for the end time of the simulation, the start" 
    echo "time will be taken from the name of the input gro file "
    echo "Use absolute paths to run on cluster!!"
    echo ""
    exit 1
fi

# gromacs dir
arch=`uname -m`
gmxDir=/project/StruPPi/Software/gromacs-3.3.3/i686/bin
if [ "$arch" == "x86_64" ] 
then
    gmxDir=/project/StruPPi/Software/gromacs-3.3.3/x86_64/bin
fi

# mdrun command (use mpi or not)
mdrunCmd="$gmxDir/mdrun"
if [ "$np" -gt "1" ]
then
    mdrunCmd="mpirun -np $np $gmxDir/mdrun_mpi"
fi

# mdp files

mdpMD=$MDP_MD
if [ -e "$outDir/md.mdp" ]
then
	mdpMD="md.mdp"
	echo "Using $mdpMD in $outDir as parameter file for MD step" 
fi


outLog=out.log
cmdLog=cmd.log


# FUNCTIONS

# checkExitStatus, takes 2 parameters:
# exit value
# program name
checkExitStatus () {
	exitVal="$1"
	prog="$2"
	if [ "$exitVal" -ne "0" ]
	then
		echo "$prog failed. Revise log file $outLog. Exiting"
		exit 1
	fi 	
}

checkFilesExist () {
	while [ ! -z "$1" ]
	do
		if [ ! -e "$1" ]
		then
			echo "File $1 does not exist in $outDir"
			exit 1
		fi
    	shift
	done	
}

# runSimulation, takes 5 parameters:
# mdp file
# in gro file
# in top file
# suffix of the input gro file (without the .gro), e.g. for T0123.pr.gro the suffix is "pr" 
# suffix of the output files, e.g. if input is T0123.pr.gro and we specify "md" as suffix output, we'd get output files like T0123.md.gro or T0123.md.xtc 
runSimulation () {
   	if [ -z "$5" ]
   	then
    	echo "Not enough parameters were passed to runSimulation. This is a bug in this script... exiting"
    	exit 1
   	fi

	# globals: outLog, cmdLog, gmxDir, mdrunCmd, np
	
	# parameters	
 	mdp="$1"
 	gro="$2"
 	top="$3"
	suffixIn="$4"
	suffixOut="$5"
 	
 	# locals: beware!! these cannot conflict with a name already in use elsewhere in the script
 	bn1=`basename $gro .gro`
 	bn2=`basename $bn1 .$suffixIn`
	tpr=$bn2.$suffixOut.tpr
 	trr=$bn2.$suffixOut.trr
 	outGro=$bn2.$suffixOut.gro
 	edr=$bn2.$suffixOut.edr
 	xtc=$bn2.$suffixOut.xtc
 	log=$bn2.$suffixOut.log
 	
 	# grompp creates a mdout.mdp file by default, which is (I think) useless
 	# If there is already one in the outdir it will back it up giving the backed up file a version number. 
 	# But gromacs won't backup up beyond version 128: because of that we remove the file here before running grompp
 	rm -f mdout.mdp
 	
	echo -e "##########\n# preparing run with input $bn1.gro and output $outGro (grompp)\n##########" >> $outLog
	echo "$gmxDir/grompp -f $mdp -c $gro -p $top -o $tpr -np $np" >> $cmdLog
	$gmxDir/grompp -f $mdp -c $gro -p $top -o $tpr -np $np 1>> $outLog 2>&1
	
	echo -e "##########\n# mdrun with input $bn1.gro and output $outGro (mdrun)\n##########" >> $outLog
	echo "$mdrunCmd -s $tpr -o $trr -c $outGro -e $edr -x $xtc -g $log -np $np" >> $cmdLog
	$mdrunCmd -s $tpr -o $trr -c $outGro -e $edr -x $xtc -g $log -np $np 1>> $outLog 2>&1

	checkExitStatus $? mdrun	
}



# END FUNCTIONS

# gromacs assumes that its output directory is the same as the current, a few things don't work well if that's not the case
# thus we change to the outDir before starting to run any commands, 
# that means we use directly the file names without paths for all input/output files


cd $outDir
echo "" > $outLog
echo "" > $cmdLog



inGro=`basename $inGro` # we strip the directory off inGro, it simply must be in outDir 
# top file
basename1=`echo "$inGro" | perl -ne 'print $1 if /^(\w+)\./'`
topFile=$basename1.top

# we check if files are in $outDir, if not we exit
checkFilesExist $inGro $topFile


begT=`echo $inGro | perl -ne 'print $1 if /^\w+\.(\d\d\d\d)ns.*\.gro/'`
if [ -z "$begT" ]
then
	echo "Input gro file is not in the right format, it must contain the ns end point of the simulation and have the .gro extension, e.g. Txxxx.0001ns.md.gro"
	exit 1
fi

begT=`expr $begT + 1` # this strips off the leading 0s and sum 1 

# MD simulation
echo "Molecular dynamics simulation"
lastSuffix=""

for et in `seq $begT $endT`
do
	bt=`expr $et - 1`
	btStr=`printf %04.0f $bt`ns
	etStr=`printf %04.0f $et`ns
	echo "Running simulation up to $etStr"
		
	runSimulation $mdpMD $basename1.$btStr.md.gro $topFile $btStr.md $etStr.md
	lastSuffix="$etStr.md"
done
	


# 8. converting the final md gro file to pdb
mdFile=$basename1.$lastSuffix.gro
tprMD=$basename1.$lastSuffix.tpr
# final pdb file
mdPdbFile=$basename1.$lastSuffix.pdb


echo "Converting $mdFile to $mdPdbFile"
echo -e "##########\n# converting gro to pdb (trjconv)\n##########" >> $outLog
echo "echo 1 | $gmxDir/trjconv -f $mdFile -o $mdPdbFile -s $tprMD" >> $cmdLog
echo 1 | $gmxDir/trjconv -f $mdFile -o $mdPdbFile -s $tprMD 1>> $outLog 2>&1
checkExitStatus $? trjconv
