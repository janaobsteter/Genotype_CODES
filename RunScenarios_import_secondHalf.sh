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
			sed "s/_PERIMPORTK_/$perimportK/g" Import_generic_secondHalf.sh > Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_PERIMPORTBM_/$perimportBM/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_REP_/$rep/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_STRATEGY_/$strategy/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_SCENARIOHOME_/$scenarioHome/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_SCENARIOIMPORT_/$scenarioImport/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_TRAITHOME_/$traitHome/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			sed -i "s/_TRAITIMPORT_/$traitImport/g" Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
			qsub Import_secondhalf_${strategy}_${scenarioHome}${scenarioImport}${rep}_${traitHome}${traitImport}_${perimport}.sh
	done
done

