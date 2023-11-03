#the first argument is the replication number
for rep in $(seq 21 30)
do
	sed "s/_Replicate_/$rep/g" SimulationEddie_BurnIn.sh > SimulationEddie_${rep}BurnIn.sh
	qsub SimulationEddie_${rep}BurnIn.sh
done

