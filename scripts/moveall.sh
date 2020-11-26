#!/usr/bin/bash

for file in ./static_tars/*;do
	name=$(basename $file)
	id=${name::4}
	mv $file ./static/$id/
done
	
