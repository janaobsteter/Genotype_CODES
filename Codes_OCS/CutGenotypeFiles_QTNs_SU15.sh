#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N CutGenotypeFiles_QTNs_SU15
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/CUT_Genotype
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -Cut_GENOTYPES
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
bash /home/v1jobste/JanaO/CutGenotypes_QTNs_SU15.sh

