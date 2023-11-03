
#for rep in $(seq 0 19)
#do#
#	cp ./Essentials/selection10.py FillInBurnIn$rep
#done

for str in "SU55" "SU51" "SU15"
do
	for sc in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
	do
		for rep in $(seq 0 19)
		do
			cp ./Essentials/selection10.py $str/$sc$rep/
		done
	done
done
