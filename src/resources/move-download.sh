#!/usr/bin/bash

staticfiles=/home/rtviii/dev/ribxz/static/

targetDir=$1

for file in "$targetDir"/*.json

do
	if [ -f "$file" ]; then

	pdbid=$(basename $file)
	pdbid=${pdbid::4}
	pdbid=${pdbid^^}

	if [ ! -d $staticfiles$pdbid ]; then
		echo "$staticfiles$pdbid does not exist. Creating.."
		mkdir -p "$staticfiles$pdbid"
		fi

	getcif $pdbid $staticfiles$pdbid 
	cp $file  $staticfiles$pdbid

	fi
done

