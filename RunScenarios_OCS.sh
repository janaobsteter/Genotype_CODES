for degree in  15 30 45 60 75
do
	sed "s/DEGREE/$degree/g" SimulationEddie_0Gen_OCS.sh > SimulationEddie_${degree}GenOCS.sh
	qsub SimulationEddie_${degree}GenOCS.sh
done

