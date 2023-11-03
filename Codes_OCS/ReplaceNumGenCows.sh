for scenario in "Clas" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
  	sed -i "s/'k', 25.0/'k', 2.5/g" SelectionParam_${scenario}.csv 
done

