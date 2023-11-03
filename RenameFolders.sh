criteria=$1

for rep in $(seq 0 10)
do
	for scenario in Class GenSLO OtherCowsGen BmGen Gen
	do
		mv $scenario$rep $scenario$rep_$criteria
	done
done
