#!/bin/bash

while read p; do
	ribxz --generateStructureProfile "$p"
done < candidates.txt

