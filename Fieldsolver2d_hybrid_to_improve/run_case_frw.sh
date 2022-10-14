mkdir -p output/frw
mkdir -p output/frw/backup

if [ $# -ne 0 ]; then

	for ((index=$1; index <= $2; index++))
	do
		ALL="${ALL} ${index}"
	done
else
	ALL="1 2 3 4 5 6 7 8"
fi
#ALL="74"

for NUM in $ALL; do
	echo "######################## $NUM ########################"
	if [ -f output/frw/${NUM}.out ]; then
		mv output/frw/${NUM}.out output/frw/backup
	fi
	./build/fieldsolver2d --method frw --in input/input_${NUM}.data --out output/frw/output_${NUM}.out ${OTHERS}
#	echo "--------------  reference value  -------------"
#	cat output/frw/output_${NUM}.out
done
