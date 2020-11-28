#!/usr/bin/bash

specfile=/home/rxz/dev/ribxz/ciftools/scripts/species/562.csv
IFS=','

[ ! -f $specfile ] && { echo "$specfile is not found"; exit -1; }


while read pdb
do
	 p3 parallel-test.py ${pdb//\"/}
done < $specfile

