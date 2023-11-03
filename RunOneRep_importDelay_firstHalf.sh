rep=$1
trait2=$2
bash RunScenarios_import_firstHalf.sh $rep 0 0 "Class" "Gen" 1 $trait2
bash RunScenarios_import_firstHalf.sh $rep 10 10 "Class" "Gen" 1 $trait2
bash RunScenarios_import_firstHalf.sh $rep 25 25 "Class" "Gen" 1 $trait2
bash RunScenarios_import_firstHalf.sh $rep 50 50 "Class" "Gen" 1 $trait2
bash RunScenarios_import_firstHalf.sh $rep 100 100 "Class" "Gen" 1 $trait2
bash RunScenarios_import_firstHalf.sh $rep 0 100 "Class" "Gen" 1 $trait2
