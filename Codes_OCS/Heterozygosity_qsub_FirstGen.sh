#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Heterozygosity
#$ -cwd       
#$ -l h_vmem=64G
#$ -pe sharedmem 1
#$ -o /home/v1jobste/JanaO/Heterozygosity.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -Heterozygosity_Firstgen
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 

/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript Heterozygosity_FirstGen.R > /home/v1jobste/JanaO/Heterozygosity.txt

