#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Hmatrix_prof
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/Hmatrix_prof.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -janaSimulacije
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
#load R module
module load R/3.3.2

/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript CreateHmatrix.R > /home/v1jobste/JanaO/Prof_RHmatrix.txt

