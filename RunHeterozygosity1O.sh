for strategy in "SU55" "SU51" "SU15"
do
        for scenario in "OtherCowsGen"
	do
		sed "s/_strategy_/$strategy/g" Heterozygosity_Replace1.R > tmp.R
		sed "s/_scenario_/$scenario/g" tmp.R > Heterozygosity_${strategy}_${scenario}1.R
		sed "s/_heterozygosityScript_/Heterozygosity_${strategy}_${scenario}1.R/g" Heterozygosity_qsub_Rep.sh  > tmp.sh
		sed "s/_SCENARIO_/${strategy}_${scenario}1/g" tmp.sh > Heterozygosity_${strategy}_${scenario}.sh
		qsub Heterozygosity_${strategy}_${scenario}.sh
	done
done
