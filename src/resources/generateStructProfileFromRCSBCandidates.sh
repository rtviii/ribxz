#!/bin/bash

CANDIDATESFILE=$1

while read structurePDBID; do
	ribxyz  --generateStructureProfile $structurePDBID
done < $1

