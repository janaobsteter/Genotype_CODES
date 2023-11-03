for strategy in "SU55" "SU51" "SU15"
do
        for scenario in  "Class" "GenSLO" "BmGen" "Gen"
	do
		sed "s/_strategy_/$strategy/g" Heterozygosity_ReplaceQTNs.R > tmp.R
		sed "s/_scenario_/$scenario/g" tmp.R > Heterozygosity_${strategy}_${scenario}_QTN.R
		sed "s/_heterozygosityScript_/Heterozygosity_${strategy}_${scenario}_QTN.R/g" Heterozygosity_qsub_Rep.sh  > tmp.sh
		sed "s/_SCENARIO_/${strategy}_${scenario}/g" tmp.sh > Heterozygosity_${strategy}_${scenario}_QTN.sh
		qsub Heterozygosity_${strategy}_${scenario}_QTN.sh
	done
done

