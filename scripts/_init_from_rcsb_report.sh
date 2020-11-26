#!/usr/bin/bash

#cat-pipe your report here


itexists()
{
echo "Here is the $1"
}

awk 'NR>1 {
gsub(/"/,"",$1);
gsub(",","",$1);
system("ribxz -struct " $1)
}'  

