for rep in $(seq 0 9)
do
  	for degree in 15 30 45 60 75
        do
          	cp SimulationEddie_${degree}_${rep}GenOCS.sh Gen${rep}_${degree}OCS/
        done
done
