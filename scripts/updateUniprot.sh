#!/bin/sh
# Script to update a local uniprot copy with sequence files and generate blast dbs for them (formatdb)

if [ -z "$1" ]
then
	echo "Usage: $0 <base local dir>"
	exit 1
fi
		

LOCALDIR=$1
CURRENT="$LOCALDIR/current"
DOWNLOAD="$LOCALDIR/download"

FORMATDB=/usr/bin/formatdb

#SITE="ftp://ftp.uniprot.org/pub" # US main ftp
SITE="ftp://ftp.ebi.ac.uk/pub" # UK mirror
# the swiss mirror doesn't seem to update properly, not using it anymore
#SITE="ftp://ftp.expasy.org" # swiss mirror
 

COMPLETEKBDIR="databases/uniprot/current_release/knowledgebase/complete"
UNIREFDIR="databases/uniprot/uniref/uniref100"

SIFTSPDB2UNIPROTFTP="ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/text/pdb_chain_uniprot.lst"


SPROT="uniprot_sprot.fasta"
SPROTGZ="${SPROT}.gz"
TREMBL="uniprot_trembl.fasta"
TREMBLGZ="${TREMBL}.gz"
UNIREF100="uniref100.fasta"
UNIREF100GZ="uniref100.fasta.gz"
ALL="uniprot_all.fasta"
RELDATEFILE="reldate.txt"
SIFTSPDB2UNIPROT="pdb_chain_uniprot.lst"

sproturl="$SITE/$COMPLETEKBDIR/$SPROTGZ"
tremblurl="$SITE/$COMPLETEKBDIR/$TREMBLGZ"
uref100url="$SITE/$UNIREFDIR/$UNIREF100GZ"
reldateurl="$SITE/$COMPLETEKBDIR/$RELDATEFILE"

# remove existing download directory if there was one
rm -rf $DOWNLOAD
# create the download dir
mkdir $DOWNLOAD

# getting the release date file if newer available
release=""
curl -z $CURRENT/$RELDATEFILE $reldateurl > $DOWNLOAD/$RELDATEFILE
if [ -s "$DOWNLOAD/$RELDATEFILE" ]
then
	release=`head -1 $DOWNLOAD/$RELDATEFILE | sed "s/UniProt Knowledgebase Release \(...._..\).*/\1/"`
	echo "New uniprot release $release available. Downloading files."
else
	echo "No new uniprot release available. Exiting"
	rm -rf $DOWNLOAD
	exit 0
fi


# download if newer available
curl -z $CURRENT/$TREMBL $tremblurl > $DOWNLOAD/${TREMBL}.gz
if [ -s "$DOWNLOAD/${TREMBL}.gz" ]
then
	echo "New trembl version downloaded"
else
	echo "Remote trembl file not newer than local one. Something wrong. Exiting."
	exit 1
fi

curl -z $CURRENT/$SPROT $sproturl > $DOWNLOAD/${SPROT}.gz
if [ -s "$DOWNLOAD/${SPROT}.gz" ]
then
	echo "New sprot version downloaded"
else
	echo "Remote sprot file not newer than local one. Something wrong. Exiting."
	exit 1	
fi

curl -z $CURRENT/$UNIREF100 $uref100url > $DOWNLOAD/${UNIREF100}.gz
if [ -s "$DOWNLOAD/${UNIREF100}.gz" ]
then
    echo "New Uniref100 version downloaded"
else
    echo "Remote Uniref100 file not newer than local one. Something wrong. Exiting"
    exit 1
fi

# getting the SIFTS PDB to UNIPROT mapping file
curl $SIFTSPDB2UNIPROTFTP > $DOWNLOAD/$SIFTSPDB2UNIPROT


# uncompressing
gzip -df $DOWNLOAD/${SPROT}.gz
gzip -df $DOWNLOAD/${TREMBL}.gz
gzip -df $DOWNLOAD/${UNIREF100}.gz
# creating the "all" file
cat $DOWNLOAD/$TREMBL $DOWNLOAD/$SPROT > $DOWNLOAD/$ALL
	
# run formatdb
# formatdb appends the path used to run it to the .pal index file, 
# thus if the path used is an absolute path it's effectively hard coding 
# them making the directory not movable. That's why we have to cd to the
# DOWNLOAD dir first, so that there's no hard-coded paths in the .pal file


echo "Running formatdb..."

#formatdb log file
logfile="$DOWNLOAD/formatdb.log"

cd $DOWNLOAD
#$FORMATDB -p T -o T -l $logfile -i $SPROT
#$FORMATDB -p T -o T -l $logfile -i $TREMBL
$FORMATDB -p T -o T -l $logfile -i $ALL
$FORMATDB -p T -o T -l $logfile -i $UNIREF100

#renaming DOWNLOAD dir to uniprot version and updating current symlink
echo "Creating new symlink..."
mv $DOWNLOAD $LOCALDIR/uniprot_$release
rm -f $CURRENT
cd $LOCALDIR
ln -s uniprot_$release current

echo "Done"
