#!/bin/sh
function usage {
	echo
	echo "Script to assign id to each possible single model type. "
	echo
	echo "Usage: $1 -d <db> "
	echo "<db>	the database where single model types will be stored"
	echo
} 

#
#Set default value for variables
#
db=""
table="single_model"
file="single_model.txt"

while getopts d: opt
do
  case "$opt" in
  	d) db="$OPTARG";;	
  esac
done

if [ -z "$db" ]
then
	echo "Missing arguments"
	usage $0
	exit 1
fi

arch=`uname -m`
case "$arch" in
    i686)
	mysqldir=/project/tla/dist/mysql-i686
    ;;
    x86_64)
	mysqldir=/project/tla/dist/mysql
    ;;
    *)
	mysqldir=/project/tla/dist/mysql-i686
    ;;
esac

mysqlbin=$mysqldir/bin/mysql
h=white

$mysqlbin -pnieve -h $h -B -N $db <<ENDSQL
SET sql_mode = "NO_UNSIGNED_SUBTRACTION,TRADITIONAL";
SELECT 'Creating the single model table if not exists ...';
CREATE TABLE IF NOT EXISTS $table (
	single_model_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
  	dist DECIMAL(4,2) NOT NULL,
  	expBB TINYINT NOT NULL,
	CW VARCHAR(16) NOT NULL,
	CT VARCHAR(16) NOT NULL,
	CR VARCHAR(80) NOT NULL,
	d TINYINT(1) NOT NULL,
	UNIQUE INDEX DESCRIPTION_IDX (dist, expBB, CW, CT, CR, d)
) ENGINE=MyISAM;
ENDSQL

contactType="Ca Cb Ca/Cb Ca+Cb+Ca/Cb BB SC BB/SC BB+SC+BB/SC Ca+SC"

#(true) 									: no filtering
#((i_cid!=j_cid)OR(abs(i_num-j_num)>=2))	: remove backbone connectivity 
#((i_cid!=j_cid)OR(abs(i_num-j_num)>=4))	: remove short range edges
#((i_cid!=j_cid)OR(abs(i_num-j_num)>=10))	: remove short and medium range edges
contactRange="(true) \
			((i_cid!=j_cid)OR(abs(i_num-j_num)>=2)) \
			((i_cid!=j_cid)OR(abs(i_num-j_num)>=4)) \
			((i_cid!=j_cid)OR(abs(i_num-j_num)>=10)) \
			((i_sstype!=j_sstype)OR(i_ssid!=j_ssid))"

echo "Outputing the models..."
echo -n "" > $file
first=`seq 2 0.5 4`
fours=`seq 4.1 0.1 4.9`
last=`seq 5 0.5 20`
contactDistance=`echo "$first $fours $last"`

for dist in $contactDistance
do
	for ct in $contactType
	do
		for cr in $contactRange
		do				
			for expBB in -1 0
			do
				#unweighted, undirected
				echo -e "${dist}\t${expBB}\t1\t${ct}\t${cr}\t0" >> $file
				#unweighted, directed
				echo -e "${dist}\t${expBB}\t1\t${ct}\t${cr}\t1" >> $file
				#weighted, undirected
				echo -e "${dist}\t${expBB}\t${ct}\t${ct}\t${cr}\t0" >> $file
				#weighted, directed
				echo -e "${dist}\t${expBB}\t${ct}\t${ct}\t${cr}\t1" >> $file
			done
		done
	done
done

$mysqlbin -pnieve -h $h -B -N $db <<ENDSQL
SELECT 'Loading the data...';
SET sql_mode = "NO_UNSIGNED_SUBTRACTION,TRADITIONAL";
LOAD DATA LOCAL INFILE "$file" IGNORE INTO TABLE $table (dist, expBB, CW, CT, CR, d);
ENDSQL

rm -f $file
