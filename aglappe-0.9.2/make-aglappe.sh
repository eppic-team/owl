#!/bin/sh
if [ -z "$2" ]
then
    echo "usage: make-aglappe.sh <tempdir> <aglappetag>"
    echo "if instead of a tag, you want the code from trunk, just specify 'trunk' instead of the tag name"
    exit
fi


tempdir=$1
aglappetag=$2

CLASSPATH=.:/project/StruPPi/jars/mysql-connector-java.jar:/project/StruPPi/jars/JRclient-RE817.jar:/project/StruPPi/jars/java-getopt-1.0.13.jar:/project/StruPPi/jars/junit-3.8.1.jar:/project/StruPPi/jars/commons-codec-1.3.jar:/project/StruPPi/jars/xmlrpc-client-3.0.jar:/project/StruPPi/jars/xmlrpc-common-3.0.jar:/project/StruPPi/jars/ws-commons-util-1.0.1.jar:/project/StruPPi/jars/vecmath.jar:/project/StruPPi/jars/Jama-1.0.2.jar:/project/StruPPi/jars/jaligner.jar:/project/StruPPi/jars/jung/collections-generic-4.01.jar:/project/StruPPi/jars/jung/jung-api-2.0-alpha2.jar:/project/StruPPi/jars/jung/jung-graph-impl-2.0-alpha2.jar

cd $tempdir

if [ -e "$aglappetag" ]
then
    echo "File exists with name $cmviewtag or name $aglappetag, can't create directory"
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


# compiling
echo "Compiling..."
cd $aglappetag
javac proteinstructure/*.java tools/*.java tinker/*.java sadp/*.java

# creating jar file
echo "Creating jar file: $aglappetag.jar ..."
jar -cfm ../$aglappetag.jar Manifest.txt .

# removing $aglappetag temp directory
cd ..
rm -rf $aglappetag
