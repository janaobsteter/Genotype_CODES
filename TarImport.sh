for strategy in "GenGen" "ClassGen"
do
	for scenario in "100_100" "50_50" "25_25" "10_10" "0_100" "0_0"
	do
		for rep in $(seq 7 9)
		do
			for trait in 2 3
			do
				sed "s/_STRATEGY_/${strategy}/g" CompressImport_generic.sh > CompressImport_${strategy}${rep}_${scenario}_${trait}.sh
				sed -i "s/_SCENARIO_/${scenario}/g" CompressImport_${strategy}${rep}_${scenario}_${trait}.sh
				sed -i "s/_REP_/${rep}/g" CompressImport_${strategy}${rep}_${scenario}_${trait}.sh
				sed -i "s/_TRAIT_/${trait}/g" CompressImport_${strategy}${rep}_${scenario}_${trait}.sh
				qsub CompressImport_${strategy}${rep}_${scenario}_${trait}.sh
			done
		done
	done
done
