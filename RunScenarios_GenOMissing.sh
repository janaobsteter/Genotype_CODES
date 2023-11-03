for strategy in "SU55" "SU15" "SU51"
do
	for rep in 10 12
	do
		sed "s/_Replicate_/$rep/g" SimulationEddie_Gen.sh > SimulationEddie_$repGen.sh
		sed "s/_Strategy_/$strategy/g" SimulationEddie_$repGen.sh > SimulationEddie_${rep}Gen_${strategy}.sh
		qsub SimulationEddie_${rep}Gen_${strategy}.sh
	done
done

