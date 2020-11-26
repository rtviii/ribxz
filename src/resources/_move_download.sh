#!/usr/bin/bash

staticfiles=/home/rtviii/dev/ribxz/static/

file=$1
# for each json file in target direcotry

if [ -f "$file" ]; then
	#get the pdbid, cap
	pdbid=$(basename $file)
	pdbid=${pdbid::4}
	pdbid=${pdbid^^}

	#if its folder in static doesn't yet exist, create
	if [ ! -d $staticfiles$pdbid ]; then
		echo "$staticfiles$pdbid does not exist. Creating.."
		mkdir -p $staticfiles$pdbid 
	fi

#download the cif file for it 
getcif $pdbid $staticfiles$pdbid 
mv $file  "$staticfiles$pdbid"
fi

