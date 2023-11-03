#the first argument is the replication number
sed "s/_Replicate_/$1/g" SimulationEddie_BurnIn.sh > SimulationEddie_$1BurnIn.sh
qsub SimulationEddie_$1BurnIn.sh

