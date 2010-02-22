#!/bin/sh
if [ -z "$3" ]
then
    echo "Usage: make-owl.sh <svn base url> <tempdir> <tag>"
    echo "If instead of a tag, you want the code from trunk, just specify 'trunk' instead of the tag name"
    echo "e.g. for svn base url: svn://black/aglappe or svn://www.bioinformatics.org/svnroot/owl"
    exit
fi

svnbaseurl=$1
tempdir=$2
tag=$3

echo "Compiling with:"
javac -version
echo ""

CLASSPATH=.:\
/project/StruPPi/jars/mysql-connector-java.jar:\
/project/StruPPi/jars/JRclient-RE817.jar:\
/project/StruPPi/jars/java-getopt-1.0.13.jar:\
/project/StruPPi/jars/junit-3.8.1.jar:\
/project/StruPPi/jars/commons-codec-1.3.jar:\
/project/StruPPi/jars/xmlrpc-client-3.1.jar:\
/project/StruPPi/jars/xmlrpc-common-3.1.jar:\
/project/StruPPi/jars/ws-commons-util-1.0.2.jar:\
/project/StruPPi/jars/vecmath.jar:\
/project/StruPPi/jars/Jama-1.0.2.jar:\
/project/StruPPi/jars/jaligner.jar:\
/project/StruPPi/jars/drmaa.jar:\
/project/StruPPi/jars/jung/collections-generic-4.01.jar:\
/project/StruPPi/jars/jung/jung-2.0.1/jung-api-2.0.1.jar:\
/project/StruPPi/jars/jung/jung-2.0.1/jung-graph-impl-2.0.1.jar:\
/project/StruPPi/jars/jung/jung-2.0.1/jung-algorithms-2.0.1.jar

cd $tempdir

if [ -e "$tag" ]
then
    echo "File exists with name $tag, can't create directory"
    exit 1
fi

# exporting from svn
echo "Exporting source from svn"


if [ "$tag" = "trunk" ]
then
    tag="owl-trunk"
    svn export $svnbaseurl/trunk/ $tag
else
    svn export $svnbaseurl/tags/$tag
fi


# compiling
echo "Compiling..."
cd $tag/src
javac *.java proteinstructure/*.java tools/*.java tinker/*.java sadp/*.java sequence/*.java actionTools/*.java ppi/*.java

# creating jar file
echo "Creating jar file: $tag.jar ..."
jar -cfm ../../$tag.jar ../Manifest.txt .

# removing $tag temp directory
cd ../..
rm -rf $tag
