for rep in $(seq 0 20)
do
	num=$(less FillInBurnIn${rep}/IndForGeno_10000.txt | wc -l)
	echo $rep
	echo $num
done
