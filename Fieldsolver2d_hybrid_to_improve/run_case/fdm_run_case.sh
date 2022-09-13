mkdir -p output
mkdir -p output/backup

if [ $# -ne 0 ]; then
	index=1
	ALL="${!index}"
	for ((index=2; index <= $#; index++))
	do
		OTHERS="${OTHERS} ${!index}"
	done
else
	ALL="1 2 3 4 5 6 7 8"
fi

for NUM in $ALL; do
	echo "######################## $NUM ########################"
	if [ -f output/${NUM}.out ]; then
		mv output/${NUM}.out output/backup
	fi
	../build/fieldsolver2d --method fdm --in ../input/input_${NUM}.data --out output/output_${NUM}.out ${OTHERS}
	echo "--------------  reference value  -------------"
	cat ../input/output_${NUM}.out
done
