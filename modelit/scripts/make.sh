#!/bin/sh
if [ -z "$3" ]
then
    echo "usage: $0 <tempdir> <aglappetag> <casptag>"
    echo "if instead of a tag, you want the code from trunk, just specify 'trunk' instead of the tag name"
    exit
fi


tempdir=$1
aglappetag=$2
casptag=$3

# we force compilation in java 1.5 to keep backwards compatibility
#JAVAVERSION=1.5.0
# now we need java 6 for setExecutable() call in TinkerRunner
JAVAVERSION=1.6.0

CLASSPATH=.:\
/project/StruPPi/jars/java-getopt-1.0.13.jar:\
/project/StruPPi/jars/mysql-connector-java.jar:\
/project/StruPPi/jars/vecmath.jar:\
/project/StruPPi/jars/Jama-1.0.2.jar:\
/project/StruPPi/jars/jaligner.jar:\
/project/StruPPi/jars/jung/collections-generic-4.01.jar:\
/project/StruPPi/jars/jung/jung-api-2.0-beta1.jar:\
/project/StruPPi/jars/jung/jung-graph-impl-2.0-beta1.jar:\
/project/StruPPi/jars/drmaa.jar

cd $tempdir

if [ -e "$casptag" ] || [ -e "$aglappetag" ]
then
    echo "File exists with name $casptag or name $aglappetag, can't create directory"
    exit 1
fi

# exporting from svn
echo "Exporting source from svn"


if [ "$aglappetag" = "trunk" ]
then
    aglappetag="aglappe-trunk"
    svn export file:///project/StruPPi/svn/aglappe/trunk/ $aglappetag
else
    svn export file:///project/StruPPi/svn/aglappe/tags/$aglappetag
fi

if [ "$casptag" = "trunk" ]
then
    casptag="casp-trunk"
    svn export file:///project/StruPPi/svn/Casp/trunk/src $casptag
else
    svn export file:///project/StruPPi/svn/Casp/tags/$casptag/src
fi


# compiling
echo "Compiling..."
cd $casptag
cp -R ../$aglappetag/proteinstructure .
cp -R ../$aglappetag/tools .
cp -R ../$aglappetag/tinker .
cp -R ../$aglappetag/sadp .
cp -R ../$aglappetag/sequence .
cp -R ../$aglappetag/graphAveraging .

javac casp/pipeline/*.java casp/benchmarking/*.java casp/tools/*.java

# creating jar file
echo "Creating jar file: $casptag.jar ..."
jar -cf ../$casptag.jar .

# removing temp directories
cd ..
rm -rf $aglappetag
rm -rf $casptag
