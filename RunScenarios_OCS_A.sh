for rep in 0 ##$(seq 0 1)
	do
	for degree in  15 75 ##15 30 45 60 75
	do
		sed "s/DEGREE/$degree/g" SimulationEddie_0GenOCS_SU55_A.sh > SimulationEddie_${degree}GenOCS_A.sh
		sed "s/_REPEAT_/$rep/g" SimulationEddie_${degree}GenOCS_A.sh > SimulationEddie_${degree}_${rep}GenOCS_A.sh
		qsub -pe sharedmem 8 -l h_vmem=4G SimulationEddie_${degree}_${rep}GenOCS_A.sh
	done
done

