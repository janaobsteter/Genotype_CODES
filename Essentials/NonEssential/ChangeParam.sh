for scenario in "Class" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
	sed -i 's/pbn,4/pbn,1/g' SelectionParam_${scenario}.csv
	sed -i 's/pozitivnoTestDoz,400/pozitivnoTestDoz,2000/g' SelectionParam_${scenario}.csv
	sed -i 's/genpbn,5/genpbn,1/g' SelectionParam_${scenario}.csv
done
