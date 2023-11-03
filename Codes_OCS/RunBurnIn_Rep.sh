#the first argument is the replication number
for rep in $(seq 0 19)
do
	sed "s/_Replicate_/$rep/g" SimulationEddie_BurnIn.sh > SimulationEddie_${rep}BurnIn.sh
	qsub SimulationEddie_${rep}BurnIn.sh
done

