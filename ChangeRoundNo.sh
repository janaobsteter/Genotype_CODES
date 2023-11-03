for rep in $(seq 0 9)
do
  	for degree in 15 30 45 60 75
        do
		cd Gen${rep}_${degree}OCS/
		gen=$(echo $(tail GenPed_EBV.txt | cut -f1 -d"," | sort | uniq))
		gen1=$(($gen - 20 + 1))
                sed -i "s/range(21,/range(${gen1},/g" Eddie_SelectionRep_Selection_OCS.py 
		qsub SimulationEddie_${degree}_${rep}GenOCS.sh
		cd /home/v1jobste/JanaO/

        done
done

