mkdir -p output/fdm
mkdir -p output/fdm/backup

#if [ $# -ne 0 ]; then
#	index=1
#	ALL="${!index}"
#	for ((index=2; index <= $#; index++))
#	do
#		OTHERS="${OTHERS} ${!index}"
#	done
#else
#	ALL="1 2 3 4 5 6 7 8"
#fi

if [ $# -ne 0 ]; then

	for ((index=0; index <= $1; index++))
	do
		ALL="${ALL} ${index}"
	done
else
	ALL="1 2 3 4 5 6 7 8"
fi

for NUM in $ALL; do
	echo "######################## $NUM ########################"
	if [ -f output/fdm/${NUM}.out ]; then
		mv output/fdm/${NUM}.out output/fdm/backup
	fi
	./build/fieldsolver2d --method fdm --in input/input_${NUM}.data --out output/fdm/output_${NUM}.out ${OTHERS}
	echo "--------------  reference value  -------------"
	cat output/fdm/output_${NUM}.out
done
