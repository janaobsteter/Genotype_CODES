for strategy in "SU55"
do
	for rep in $(seq 0 2) ##$(seq 0 19) 
	do
		#sed "s/_Replicate_/$rep/g" SimulationEddie_Class_FemaleGS.sh > SimulationEddie_$repClass.sh
		#sed "s/_Strategy_/$strategy/g" SimulationEddie_$repClass.sh > SimulationEddie_${rep}Class_${strategy}.sh
		#qsub SimulationEddie_${rep}Class_${strategy}.sh	

		sed "s/_Replicate_/$rep/g" SimulationEddie_GenSLO_FemaleGS.sh > SimulationEddie_$repGenSLO.sh
		sed "s/_Strategy_/$strategy/g" SimulationEddie_$repGenSLO.sh > SimulationEddie_${rep}GenSLO_${strategy}.sh
		qsub SimulationEddie_${rep}GenSLO_${strategy}.sh
		
		sed "s/_Replicate_/$rep/g" SimulationEddie_OtherCowsGen_FemaleGS.sh > SimulationEddie_$repOtherCowsGen.sh
		sed "s/_Strategy_/$strategy/g" SimulationEddie_$repOtherCowsGen.sh > SimulationEddie_${rep}OtherCowsGen_${strategy}.sh
		qsub SimulationEddie_${rep}OtherCowsGen_${strategy}.sh

		sed "s/_Replicate_/$rep/g" SimulationEddie_BmGen_FemaleGS.sh > SimulationEddie_$repBmGen.sh
		sed "s/_Strategy_/$strategy/g" SimulationEddie_$repBmGen.sh > SimulationEddie_${rep}BmGen_${strategy}.sh		
		qsub SimulationEddie_${rep}BmGen_${strategy}.sh

		sed "s/_Replicate_/$rep/g" SimulationEddie_Gen_FemaleGS.sh > SimulationEddie_$repGen.sh
		sed "s/_Strategy_/$strategy/g" SimulationEddie_$repGen.sh > SimulationEddie_${rep}Gen_${strategy}.sh
		qsub SimulationEddie_${rep}Gen_${strategy}.sh
	done
done

