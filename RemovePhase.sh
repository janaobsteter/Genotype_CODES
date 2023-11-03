for strategy in "GenGen" "ClassGen"
do
	for scenario in "100_100" "50_50" "25_25" "10_10" "0_100" "0_0"
	do
		for rep in $(seq 0 2)
		do
			for trait in 2 3
			do
				rm -r /home/v1jobste/jobsteter/10K/SU55_import/${strategy}${rep}_${scenario}1${trait}/Chromosomes
			done
		done
	done
done
