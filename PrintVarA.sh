for rep in 0 1 2 3 4 5
do
	echo $(grep VarA FillInBurnIn${rep}/SimulatedData/TraitVarianceComponents.txt)
done
