#the first argument is the replication number

for rep in $(seq $1 $1)
do
	sed "s/_REP_/$rep/g" SimulationEddie_Import_FillIn.sh > SimulationEddie_${rep}Import_FillIn.sh
	qsub SimulationEddie_${rep}Import_FillIn.sh
done

