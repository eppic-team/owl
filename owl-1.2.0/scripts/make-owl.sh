#!/bin/sh
if [ -z "$2" ]
then
    echo "Usage: make-owl.sh <svn url> <tempdir>"
    echo "e.g. for svn url: svn://black/aglappe/trunk or svn://www.bioinformatics.org/svnroot/owl/tags/owl-1.2.0"
    exit
fi

svnurl=$1
tempdir=$2

echo "Compiling with:"
javac -version
echo ""


cd $tempdir

tag=`basename $svnurl` 


# exporting from svn
echo "Exporting source from svn"
svn export $svnurl
if [ "$?" -ne "0" ]
then
	echo "Couldn't export from svn. Exiting"
	exit 1
fi


# compiling
echo "Compiling..."
cd $tag/src


jarfiles=`ls ../jars/*.jar`
CLASSPATH="."
for jarfile in $jarfiles
do
	CLASSPATH="$jarfile:$CLASSPATH"
done

echo $CLASSPATH


javac *.java proteinstructure/*.java tools/*.java tinker/*.java sadp/*.java sequence/*.java actionTools/*.java ppi/*.java

# creating jar file
echo "Creating jar file: $tag.jar ..."
jar -cfm ../../$tag.jar ../Manifest.txt .

# removing $tag temp directory
cd ../..
rm -rf $tag