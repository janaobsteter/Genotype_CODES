
for str in "SU55" "SU51" "SU15"
do
	for sc in "Class" "GenSLO" "BmGen" "Gen" #"OtherCowsGen"
	do
		for rep in $(seq 0 19)
		do
			cp $str/$sc$rep/MeanHet_Neutral.csv ./Heterozygosity/MeanHet_Neutral_$str_$sc$rep.csv
		done
	done
done
