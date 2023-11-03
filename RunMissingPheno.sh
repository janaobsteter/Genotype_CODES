rep=$1

qsub SimulationEddie_pheno_${rep}_True1_1_1.sh
qsub SimulationEddie_pheno_${rep}_True2_1_1.sh
qsub SimulationEddie_pheno_${rep}_True5_1_1.sh
qsub SimulationEddie_pheno_${rep}_True8_1_1.sh
qsub SimulationEddie_pheno_${rep}_True1_1_2.sh
qsub SimulationEddie_pheno_${rep}_True2_1_2.sh
qsub SimulationEddie_pheno_${rep}_True5_1_2.sh
qsub SimulationEddie_pheno_${rep}_True8_1_2.sh
qsub SimulationEddie_pheno_${rep}_True9_1_2.sh
qsub SimulationEddie_pheno_${rep}_True1_2_1.sh
qsub SimulationEddie_pheno_${rep}_True2_2_1.sh
qsub SimulationEddie_pheno_${rep}_True5_2_1.sh

qsub SimulationEddie_pheno_${rep}_False1_1_1.sh
qsub SimulationEddie_pheno_${rep}_False2_1_1.sh
qsub SimulationEddie_pheno_${rep}_False5_1_1.sh
qsub SimulationEddie_pheno_${rep}_False1_1_2.sh
qsub SimulationEddie_pheno_${rep}_False2_1_2.sh
qsub SimulationEddie_pheno_${rep}_False5_1_2.sh
qsub SimulationEddie_pheno_${rep}_False8_1_2.sh
qsub SimulationEddie_pheno_${rep}_False9_1_2.sh
qsub SimulationEddie_pheno_${rep}_False1_2_1.sh
qsub SimulationEddie_pheno_${rep}_False2_2_1.sh
qsub SimulationEddie_pheno_${rep}_False5_2_1.sh

