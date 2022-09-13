#! /bin/bash
ALL="1 2 3 4 5 6 7 8"
for x in $ALL; do
	cat input/output_${x}.out output/output_${x}.out
done
