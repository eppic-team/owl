#!/bin/sh
# script to package the OWL library by exporting trunk/tag from repository, 
# making a jar file for it, adding some scripts and packing it all under a single
# zip file so that the reconstruct program can be distributed to the world 

if [ -z "$3" ]
then
	echo "usage: $0 <tmp dir> <svn url> <version>"
	echo "svn url should be like: svn://black/aglappe/trunk or svn://black/aglappe/tags/owl-1.2.0"
	exit 1
fi

roottmpdir=$1
svnurl=$2
ver=$3


scriptdir=`dirname $0`



# files
tmpdir=$roottmpdir/reconstruct-$ver
mkdir $tmpdir


svn export $svnurl/LICENSE $tmpdir/LICENSE
svn export $svnurl/gpl.txt $tmpdir/gpl.txt
svn export $svnurl/scripts/reconstruct $tmpdir/reconstruct
#svn export  $svnurl/scripts/reconstruct.bat $tmpdir/reconstruct.bat
svn export $svnurl/scripts/reconstruct.cfg $tmpdr/reconstruct.cfg

# jars
svn export $svnurl/jars $tmpdir/jars
cp $thejar $tmpdir/jars/owl.jar

# the owl jar
$scriptdir/make-owl.sh $svnurl $tmpdir
mv $tmpdir/*.jar $tmpdir/jars/owl.jar


# zipping up
cd $roottmpdir
zip -r reconstruct-$ver.zip reconstruct-$ver/