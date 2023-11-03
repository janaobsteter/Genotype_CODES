for degree in "maxGain" "minCoan"
do
	sed "s/DEGREE/$degree/g" SimulationEddie_0Gen_OCS_maxMin.sh > SimulationEddie_${degree}GenOCS_maxMin.sh
	qsub SimulationEddie_${degree}GenOCS_maxMin.sh
done

