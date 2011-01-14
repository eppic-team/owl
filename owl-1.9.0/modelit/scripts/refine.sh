#!/bin/sh
# Script for automated simulated annealing structure refinement using Gromacs

# CONSTANTS
BOXTYPE="dodecahedron" 
BOXSIDES="0.5"

MDPDIR=/project/StruPPi/Software/gromacs-3.3.3/mdp_files
# energy minimization parameters (used for ions too)
MDP_EM=$MDPDIR/em.mdp
# position restrained equilibration parameters
MDP_PR=$MDPDIR/pr.mdp
# MD simulation parameters
MDP_MD=$MDPDIR/md.mdp

SOLVENTGROUP=12

DEFAULT_NP=1

# END CONSTANTS

# command line parameters
inPdb=""
outDir=""
np=$DEFAULT_NP
begT="-1"
endT="-1"

while getopts i:o:p:b:e: opt
do
  case "$opt" in
    i) inPdb="$OPTARG";;
    o) outDir="$OPTARG";;
    p) np="$OPTARG";;
	b) begT="$OPTARG";;
	e) endT="$OPTARG";;
  esac
done

if [ -z "$inPdb" ] || [ -z "$outDir" ]
then
	echo ""
    echo "Usage: $0 "
    echo "  -i  <in pdb>"
    echo "  -o  <out dir> "
    echo " [-p] <number of processes>"
    echo " [-b] <begin time (ns)>"
    echo " [-e] <end time (ns)>"
    echo "Use begin and end time (ns) so that the md simulation will output to different files every "
    echo "1ns, instead of using 1 big file for all output."
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

mdpEM=$MDP_EM
mdpPR=$MDP_PR
mdpMD=$MDP_MD
if [ -e "$outDir/em.mdp" ]
then
	mdpEM="em.mdp"
	echo "Using $mdpEM in $outDir as parameter file for EM step" 
fi
if [ -e "$outDir/pr.mdp" ]
then
	mdpPR="pr.mdp"
	echo "Using $mdpPR in $outDir as parameter file for PR step" 
fi
if [ -e "$outDir/md.mdp" ]
then
	mdpMD="md.mdp"
	echo "Using $mdpMD in $outDir as parameter file for MD step" 
fi


outLog=out.log
cmdLog=cmd.log

basename=`basename $inPdb .pdb`

# gro files
groFile=$basename.gro
boxFile=$basename.box.gro
watFile=$basename.wat.gro
ionFile=$basename.ion.gro
emFile=$basename.em.gro
prFile=$basename.pr.gro

# final pdb file
mdPdbFile=$basename.md.pdb

# top files
topFile=$basename.top


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
 	
 	# locals
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

# addIons, takes 7 parameters:
# positive charge
# negative charge
# mdp file
# in gro file
# in top file
# suffix of the input gro file (without the .gro), e.g. for T0123.pr.gro the suffix is "pr" 
# suffix of the output files, e.g. if input is T0123.pr.gro and we specify "md" as suffix output, we'd get output files like T0123.md.gro or T0123.md.xtc 
addIons () {
   	if [ -z "$7" ]
   	then
    	echo "Not enough parameters were passed to addIons. This is a bug in this script... exiting"
    	exit 1
   	fi
	# globals: SOLVENTGROUP, mdpEM, outLog, cmdLog, gmxDir
	
	# parameters
	poscharge="$1"
	negcharge="$2"
	mdp="$3"
	gro="$4"
	top="$5"
	suffixIn="$6"
	suffixOut="$7"
	
	#locals
 	bn1=`basename $gro .gro`
 	bn2=`basename $bn1 .$suffixIn`
	tpr=$bn2.$suffixOut.tpr
 	outGro=$bn2.$suffixOut.gro
 	log=$bn2.$suffixOut.log 	
	ionOpt="-nname CL- -nn $negcharge -pname NA+ -np $poscharge"

	echo "Adding $negcharge CL- and $poscharge NA+"		
	echo -e "##########\n# preparing to add ions (grompp)\n##########" >> $outLog
	echo "$gmxDir/grompp -f $mdp -c $gro -p $top -o $tpr" >> $cmdLog
	$gmxDir/grompp -f $mdp -c $gro -p $top -o $tpr 1>> $outLog 2>&1
	echo -e "##########\n# adding ions (genion)\n##########" >> $outLog
	echo "echo $SOLVENTGROUP | $gmxDir/genion -s $tpr -o $outGro -p $top $ionOpt -g $log" >> $cmdLog
	echo $SOLVENTGROUP | $gmxDir/genion -s $tpr -o $outGro -p $top $ionOpt -g $log 1>> $outLog 2>&1
	
	checkExitStatus $? genion
}


# END FUNCTIONS

# gromacs assumes that its output directory is the same as the current, a few things don't work well if that's not the case
# thus we change to the outDir before starting to run any commands, 
# that means we use directly the file names without paths for all input/output files, except for $inPdb that can be in another directory
firstChar=`echo $inPdb | cut -c 1` # first character of inPdb
if [ "$firstChar" != "/" ] # i.e. inPdb is not an absoute path
then
	inPdb=`pwd`/$inPdb
fi


cd $outDir
echo "" > $outLog
echo "" > $cmdLog


# 1. convert pdb to gro and topology files
echo "Converting pdb to gro"
echo -e "##########\n# pdb2gmx\n##########" >> $outLog
echo "$gmxDir/pdb2gmx -f $inPdb -o $groFile -p $topFile -ff G43a1 -ignh" >> $cmdLog
$gmxDir/pdb2gmx -f $inPdb -o $groFile -p $topFile -ff G43a1 -ignh 1>> $outLog 2>&1
checkExitStatus $? pdb2gmx
charge=`cat $outLog | grep "^Total charge " | perl -pe 's/^Total charge (-?\d+)\.\d+ e/$1/'`

# 2. center the molecule within a box (with specified type and space on the sides)
echo "Creating box"
echo -e "##########\n# creating box (editconf)\n##########" >> $outLog
echo "$gmxDir/editconf -bt $BOXTYPE -f $groFile -o $boxFile -c -d $BOXSIDES" >> $cmdLog
$gmxDir/editconf -bt $BOXTYPE -f $groFile -o $boxFile -c -d $BOXSIDES 1>> $outLog 2>&1
checkExitStatus $? editconf

# 3. add solvent water to the box (-cp is solute, -cs solvent)
echo "Adding solvent to box"
echo -e "##########\n# adding solvent (genbox)\n##########" >> $outLog
echo "$gmxDir/genbox -cp $boxFile -cs spc216.gro -o $watFile -p $topFile" >> $cmdLog
$gmxDir/genbox -cp $boxFile -cs spc216.gro -o $watFile -p $topFile 1>> $outLog 2>&1
checkExitStatus $? genbox

# 4. adding ions if needed
poscharge=0
negcharge=0
if [ "$charge" -lt "0" ]
then
	# no easy way to get abs value, we use bc instead of expr
	poscharge=`echo "$charge/1*-1 + 1"|bc`
	negcharge=1
else 
	if [ "$charge" -gt "0" ]
	then
		poscharge=1
		negcharge=`expr $charge + 1`
	else
		if [ "$charge" -eq "0" ]
		then
			poscharge=1
			negcharge=1  
		fi
	fi
fi
if [ "$poscharge" -eq "0" ]
then
	# if poscharge is still 0, then none of the if conditions above executed --> $charge was not set 
    echo "Charge value couldn't be parsed. Exiting"
    exit 1
fi
echo "Total charge $charge"
addIons $poscharge $negcharge $mdpEM $watFile $topFile wat ion

# 5. minimization
echo "Energy minimization"
runSimulation $mdpEM $ionFile $topFile ion em

# 6. equilibration 
echo "Position restrained equilibration"
runSimulation $mdpPR $emFile $topFile em pr

# 7. MD simulation 
echo "Molecular dynamics simulation"
lastSuffix=""
if [ "$begT" -eq "-1" ] || [ "$endT" -eq "-1" ]
then
	runSimulation $mdpMD $prFile $topFile pr md
	lastSuffix="md"
else
	echo "Splitting output per ns"
	tStr=`printf %04.0f $begT`ns
	echo "Running simulation up to $tStr"
	runSimulation $mdpMD $prFile $topFile pr $tStr.md
	endTminus1=`expr $endT - 1`
	for bt in `seq $begT $endTminus1`
	do
		btStr=`printf %04.0f $bt`ns
		et=`expr $bt + 1`
		etStr=`printf %04.0f $et`ns
		echo "Running simulation up to $etStr"
		runSimulation $mdpMD $basename.$btStr.md.gro $topFile $btStr.md $etStr.md
		lastSuffix="$etStr.md"
	done
	
fi


# 8. converting the final md gro file to pdb
mdFile=$basename.$lastSuffix.gro
tprMD=$basename.$lastSuffix.tpr

echo "Converting $mdFile to $mdPdbFile"
echo -e "##########\n# converting gro to pdb (trjconv)\n##########" >> $outLog
echo "echo 1 | $gmxDir/trjconv -f $mdFile -o $mdPdbFile -s $tprMD" >> $cmdLog
echo 1 | $gmxDir/trjconv -f $mdFile -o $mdPdbFile -s $tprMD 1>> $outLog 2>&1
checkExitStatus $? trjconv
