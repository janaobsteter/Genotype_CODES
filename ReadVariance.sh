for rep in $(seq 0 10) 
do
	echo $rep
	tail FillInBurnIn${rep}/SimulatedData/TraitVarianceComponents.txt | grep VarA
        tail FillInBurnIn${rep}/SimulatedData/TraitVarianceComponents.txt | grep VarE

done
