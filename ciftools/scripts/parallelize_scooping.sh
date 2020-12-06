#!/bin/bash

SPECIES=$1

specfile=/home/rxz/dev/ribxz/ciftools/scripts/species/$SPECIES.csv

IFS=','

[ ! -f $specfile ] && { echo "$specfile is not found"; exit -1; }

parallel 'python3 pymol_scoop.py {1} {2} {3}' \
::: $( while read pdb;do echo ${pdb//\"/};done < $specfile ) \
::: $SPECIES \
::: 60 70 80


