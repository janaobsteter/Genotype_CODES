#the first argument is the replication number
#sed "s/_Replicate_/$1/g" SimulationEddie_Class.sh > SimulationEddie_$1Class.sh
#qsub SimulationEddie_$1Class.sh
sed "s/_Replicate_/$1/g" SimulationEddie_GenSLO.sh > SimulationEddie_$1GenSLO.sh
qsub SimulationEddie_$1GenSLO.sh
sed "s/_Replicate_/$1/g" SimulationEddie_OtherCowsGen.sh > SimulationEddie_$1OtherCowsGen.sh
qsub SimulationEddie_$1OtherCowsGen.sh
sed "s/_Replicate_/$1/g" SimulationEddie_BmGen.sh > SimulationEddie_$1BmGen.sh
qsub SimulationEddie_$1BmGen.sh
sed "s/_Replicate_/$1/g" SimulationEddie_Gen.sh > SimulationEddie_$1Gen.sh
qsub SimulationEddie_$1Gen.sh

