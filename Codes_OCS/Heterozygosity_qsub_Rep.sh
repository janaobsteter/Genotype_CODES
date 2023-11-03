#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Het_SCENARIO_
#$ -cwd       
#$ -l h_vmem=64G
#$ -pe sharedmem 1
#$ -o /home/v1jobste/JanaO/Heterozygosity__SCENARIO_.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -_SCENARIO_
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 

/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript _heterozygosityScript_ > /home/v1jobste/JanaO/Heterozygosity__SCENARIO_.txt

