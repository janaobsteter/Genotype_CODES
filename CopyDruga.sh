for scenario in 50 55
do
	for rep in $(seq 1 10)
	do
		cp SimulationEddie_${scenario}_${rep}GenOCS_drugaPol.sh  Gen${rep}_${scenario}OCS/
	done
done
