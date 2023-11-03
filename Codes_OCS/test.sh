for i in "a" "b"
do
	for rep in $(seq 0 2)
	do
		echo $i$rep
	done
done
