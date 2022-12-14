mkdir -p output/bem
mkdir -p output/bem/backup

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

	for ((index=$1; index <= $2; index++))
	do
		ALL="${ALL} ${index}"
	done
else
	ALL="1 2 3 4 5 6 7 8"
fi
#ALL="74"
for NUM
in $ALL; do
	echo "######################## $NUM ########################"
	if [ -f output/bem/${NUM}.out ]; then
		mv output/bem/${NUM}.out output/bem/backup
	fi
	./build/fieldsolver2d --method bem --in input/input_${NUM}.data --out output/bem/output_${NUM}.out ${OTHERS}
	echo "--------------  reference value  -------------"
	cat output/bem/output_${NUM}.out
done
