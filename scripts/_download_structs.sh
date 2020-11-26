#!/usr/bin/bash


for folder  in ./static/*;do
	pdbid=${folder: -4}
	pdbid=${pdbid^^}
	getcif $pdbid ./static/$pdbid/
done


