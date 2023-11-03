#the first argument is the replication number
for rep in $(seq 1 3)
do
	sed "s/_REP_/$rep/g" SimulationEddie_permEnv_BurnInNoFillIn.sh > SimulationEddie_${rep}BurnIn_permEnv.sh
	qsub SimulationEddie_${rep}BurnIn_permEnv.sh
done

