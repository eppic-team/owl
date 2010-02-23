#!/bin/sh
# shell script to run the reconstruct program (a java class in the default package of the OWL Java library)

# If the program is to be used on a cluster through the Sun Grid Engine 
# job scheduler (option -A) this variable needs to be set to wherever the
# root of the local installation of SGE is
sgeroot=/usr/local/sge

exedir=`dirname $0`

jarfiles=`ls $exedir/jars/*.jar`
CLASSPATH="."
for jarfile in $jarfiles
do
	CLASSPATH="$jarfile:$CLASSPATH"
done



if [ -e "$sgeroot/util/arch" ] 
then
	arch=`$sgeroot/util/arch`
	libstr="-Djava.library.path=$sgeroot/lib/$arch"
else
	libstr=""
fi

java $libstr reconstruct $@
