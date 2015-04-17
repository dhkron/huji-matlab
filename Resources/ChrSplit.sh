#!/bin/bash

#This will write 'domains' in the proper directory
for i in `seq 1 22`;
do
	cat CSVs/domains.IMR90.csv.12 | grep chr$i, | awk -F, '{print $2","$3}' > ../TADs/domains.IMR90.chr$i
done

