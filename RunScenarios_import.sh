repNo=$1
perimportK=$2
perimportBM=$3
scenarioHome=$4
scenarioImport=$5
traitHome=$6
traitImport=$7

for strategy in "SU55"
do
		for rep in $(seq $repNo $repNo)
		do
			sed "s/_PERIMPORTK_/$perimportK/g" Import_generic.sh > Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_PERIMPORTBM_/$perimportBM/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_REP_/$rep/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_STRATEGY_/$strategy/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_SCENARIOHOME_/$scenarioHome/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_SCENARIOIMPORT_/$scenarioImport/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_TRAITHOME_/$traitHome/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			sed -i "s/_TRAITIMPORT_/$traitImport/g" Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
			qsub Import_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimportK}${perimportBM}.sh
	done
done

