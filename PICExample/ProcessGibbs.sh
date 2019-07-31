trait=$1

mv all_solutions gibbsSol${trait}
bash Match_AFTERRenum_Gibbs${trait}.sh 
mv renumbered_Solutions Gibbs${trait}.csv
paste -d" " Gibbs${trait}.csv SampleSeq.csv > GIBBS${trait}.csv
sort -k4,4 -n GIBBS${trait}.csv > GIBBS${trait}_sort.csv

