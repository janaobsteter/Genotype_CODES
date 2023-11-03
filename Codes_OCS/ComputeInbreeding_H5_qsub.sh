#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Inbreeding_h5
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/Inbreeding_h5_Error.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -janaSimulacije
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 

/exports/cmvm/eddie/eb/groups/tier2_hickey_external/R-3.4.2/bin/Rscript CreateInbreeding_h5.R

