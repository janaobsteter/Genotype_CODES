rep=$1
for pg in "1_1" "1_2" "2_1"
do
	for control in 10 9 8 5 2 1
	do
		qsub SimulationEddie_pheno_${rep}_True${control}_${pg}.sh
	done
done
