for strategy in "SU15"
do
        for scenario in "Class"
	do
		for rep in $(seq 0 4)
		do
			sed "s/_strategy_/$strategy/g" Heterozygosity_Replace_REP.R > tmp.R
			sed -i "s/_REPEAT_/$rep/g" tmp.R 
			sed "s/_scenario_/$scenario/g" tmp.R > Heterozygosity_${strategy}_${scenario}.R
			sed "s/_heterozygosityScript_/Heterozygosity_${strategy}_${scenario}.R/g" Heterozygosity_qsub_Rep.sh  > tmp.sh
			sed "s/_SCENARIO_/${strategy}_${scenario}/g" tmp.sh > Heterozygosity_${strategy}_${scenario}.sh
			qsub Heterozygosity_${strategy}_${scenario}.sh
		done

	done
done

