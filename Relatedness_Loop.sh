for rep in 1
do
	for start in $(seq 1 10 100)
	do
		sed "s/_REPEAT_/${rep}/g" Relatedness_parts_generic.sh > Relatedness_parts_Fill${rep}_part${start}.sh
		sed -i "s/_START_/${start}/g" Relatedness_parts_Fill${rep}_part${start}.sh
		qsub Relatedness_parts_Fill${rep}_part${start}.sh
	done
done

