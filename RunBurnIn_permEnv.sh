#the first argument is the replication number

for rep in $(seq 9 9)
do
	sed "s/_REP_/$rep/g" SimulationEddie_permEnv_BurnIn.sh > SimulationEddie_${rep}permEnv_BurnIn.sh
	qsub SimulationEddie_${rep}permEnv_BurnIn.sh

done

