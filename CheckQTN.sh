for rep in $(seq 10 19)
do
	echo Gen$rep
	echo $(head Gen$rep/SimulatedData/QTNGenotype_Last20Gen.txt | wc -l)
done

