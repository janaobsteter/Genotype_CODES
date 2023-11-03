repeats=$1

for strategy in "SU55"
do
	for scenario in "Class" 
	do
		for rep in $(seq 6 6)
		do
			sed "s/_REPEATS_/$repeats/g" SimulationEddie_permEnv_generic.sh > PermEnv_${strategy}_${scenario}${rep}_${repeats}.sh
			sed -i "s/_REP_/$rep/g" PermEnv_${strategy}_${scenario}${rep}_${repeats}.sh
			sed -i "s/_STRATEGY_/$strategy/g" PermEnv_${strategy}_${scenario}${rep}_${repeats}.sh
			sed -i "s/_SCENARIO_/$scenario/g" PermEnv_${strategy}_${scenario}${rep}_${repeats}.sh
			qsub PermEnv_${strategy}_${scenario}${rep}_${repeats}.sh
		done
	done
done

