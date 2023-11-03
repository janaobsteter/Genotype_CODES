#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Jana_SU55_GenOCS_1_15
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/ErrorSU55Gen115.txt
#$ -l h_rt=80:00:00
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
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python Eddie_SelectionRep_Selection_OCS.py 1 "Gen" 15 "SU55" "10K"
#/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python PandasVersion.py
