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

if [ ! -e "$roottmpdir" ]
then
	echo "Directory $roottmpdir doesn't exist"
	exit 1
fi

# files
tmpdir=$roottmpdir/reconstruct-$ver
mkdir $tmpdir


svn export $svnurl/LICENSE $tmpdir/LICENSE
svn export $svnurl/gpl.txt $tmpdir/gpl.txt
svn export $svnurl/reconstruct/reconstruct.sh $tmpdir/reconstruct
#svn export  $svnurl/reconstruct/reconstruct.bat $tmpdir/reconstruct.bat
svn export $svnurl/reconstruct/reconstruct.cfg $tmpdir/reconstruct.cfg
svn export $svnurl/reconstruct/README-reconstruct.txt $tmpdir/README.txt
svn export $svnurl/reconstruct/sample.cm $tmpdir/sample.cm

# jars
svn export $svnurl/jars $tmpdir/jars

# the owl jar
$scriptdir/../scripts/make-owl.sh $svnurl $tmpdir
mv $tmpdir/*.jar $tmpdir/jars/owl.jar


# zipping up
cd $roottmpdir
zip -r reconstruct-$ver.zip reconstruct-$ver/