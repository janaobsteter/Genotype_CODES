#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N PermEnv_BurnIn__REP_
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/BurnInPermEnv__REP_.txt
#$ -l h_rt=100:00:00
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
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python Eddie_SelectionRep_BurnINNoFillIn_herd.py _REP_ 11 0.5,0.33333,0.33333,0.33333,1.5
