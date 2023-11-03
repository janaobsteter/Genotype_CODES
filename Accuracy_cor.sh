#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Accuracy_Orig
#$ -cwd        
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/Accuracy__strategy_.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -janaSimulacije
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
# Load Python
module load python/2.7.10

 
# Run the program
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python ComputeAccuracies_EBVandPA_Eddie.py "_strategy_"
#/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python PandasVersion.py
