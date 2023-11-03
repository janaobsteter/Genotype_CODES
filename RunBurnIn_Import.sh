#the first argument is the replication number
#second argument is traitHome
#third argument is traitImport

for rep in $(seq $1 $1)
do
	sed "s/_REP_/$rep/g" SimulationEddie_Import_BurnIn.sh > SimulationEddie_${rep}Import_BurnIn_$2$3.sh
	sed -i "s/traitHome/$2/g" SimulationEddie_${rep}Import_BurnIn_$2$3.sh
	sed -i "s/traitImport/$3/g" SimulationEddie_${rep}Import_BurnIn_$2$3.sh
	qsub SimulationEddie_${rep}Import_BurnIn_$2$3.sh
done

