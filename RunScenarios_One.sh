for rep in  12 13 18 19
do
	#sed "s/_Replicate_/$rep/g" SimulationEddie_Class.sh > SimulationEddie_$repClass.sh
	#qsub SimulationEddie_$repClass.sh
	sed "s/_Replicate_/$rep/g" SimulationEddie_GenSLO.sh > SimulationEddie_$repGenSLO.sh
	qsub SimulationEddie_$repGenSLO.sh
	#sed "s/_Replicate_/$rep/g" SimulationEddie_OtherCowsGen.sh > SimulationEddie_$repOtherCowsGen.sh
	#qsub SimulationEddie_$repOtherCowsGen.sh
	#sed "s/_Replicate_/$rep/g" SimulationEddie_BmGen.sh > SimulationEddie_$repBmGen.sh
	#qsub SimulationEddie_$repBmGen.sh
	#sed "s/_Replicate_/$rep/g" SimulationEddie_Gen.sh > SimulationEddie_$repGen.sh
	#qsub SimulationEddie_$repGen.sh
done

