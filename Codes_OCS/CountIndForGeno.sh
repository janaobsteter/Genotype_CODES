for rep in $(seq 0 30)
do
	echo $rep
	less FillInBurnIn$rep/IndForGeno_10K.txt | wc -l
done
