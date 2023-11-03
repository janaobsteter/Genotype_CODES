rep=$1
trait2=$2
bash RunScenarios_import.sh ${rep} 0 0 "Gen" "Gen" 1 ${trait2}
bash RunScenarios_import.sh ${rep} 10 10 "Gen" "Gen" 1 ${trait2}
bash RunScenarios_import.sh ${rep} 25 25 "Gen" "Gen" 1 ${trait2}
bash RunScenarios_import.sh ${rep} 50 50 "Gen" "Gen" 1 ${trait2}
bash RunScenarios_import.sh ${rep} 100 100 "Gen" "Gen" 1 ${trait2}
bash RunScenarios_import.sh ${rep} 0 100 "Gen" "Gen" 1 ${trait2}
