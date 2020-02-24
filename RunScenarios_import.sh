perimport=$1

for strategy in "SU55"
do
	for scenario in "Class" 
	do
		for rep in $(seq 1 1)
		do
			sed "s/_PERIMPORT_/$perimport/g" Import_generic.sh > Import_${strategy}_${scenario}${rep}_${perimport}.sh
			sed -i "s/_REP_/$rep/g" Import_${strategy}_${scenario}${rep}_${perimport}.sh
			sed -i "s/_STRATEGY_/$strategy/g" Import_${strategy}_${scenario}${rep}_${perimport}.sh
			sed -i "s/_SCENARIO_/$scenario/g" Import_${strategy}_${scenario}${rep}_${perimport}.sh
			qsub Import_${strategy}_${scenario}${rep}_${perimport}.sh
		done
	done
done

