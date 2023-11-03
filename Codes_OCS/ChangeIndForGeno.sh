for rep in $(seq 0 20)
do
	echo $rep
	mv FillInBurnIn$rep/IndForGeno.txt FillInBurnIn${rep}/IndForGeno_10K.txt
done
