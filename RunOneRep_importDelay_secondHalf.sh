rep=$1
trait=$2
bash RunScenarios_import_secondHalf.sh $rep 0 0 "Gen" "Gen" 1 $trait
bash RunScenarios_import_secondHalf.sh $rep 10 10 "Gen" "Gen" 1 $trait
bash RunScenarios_import_secondHalf.sh $rep 25 25 "Gen" "Gen" 1 $trait
bash RunScenarios_import_secondHalf.sh $rep 50 50 "Gen" "Gen" 1 $trait
bash RunScenarios_import_secondHalf.sh $rep 100 100 "Gen" "Gen" 1 $trait
bash RunScenarios_import_secondHalf.sh $rep 0 100 "Gen" "Gen" 1 $trait
