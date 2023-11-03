for strategy in "ClassGen" #"ClassGen" "GenGen"
do
	for scenario in "100_100" "50_50" "25_25" "10_10" "0_100" "0_0"
	do
		for rep in $(seq 7 7)
		do
			for trait in 2 3
			do
				echo ${strategy}_${scenario}_${rep}_${trait}
				head 10K/SU55_import/${strategy}${rep}_${scenario}1${trait}/GenicVariance_import.csv
				echo """"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
			done
		done
	done
done
