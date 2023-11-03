#the first argument is the replication number
for rep in $(seq 0 0)
do
	sed "s/_Replicate_/$rep/g" SimulationEddie_BurnInHerd.sh > SimulationEddie_${rep}BurnInHerd.sh
	qsub SimulationEddie_${rep}BurnInHerd.sh
done

