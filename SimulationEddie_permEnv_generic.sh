#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N _STRATEGY__permEnv_SCENARIO__REP___REPEATS_repeats
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/_STRATEGY__permEnv__SCENARIO__REP___REPEATS_repeats.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -janaSimulacije
#  runtime limit of 10 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
# Load Python
module load python/2.7.10

 
# Run the program
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python Eddie_SelectionRep_Selection_repeatedPheno.py _REP_ "_SCENARIO_" "_STRATEGY_" "10K" _REPEATS_ 1.0,0.333333,0.333333,0.333333,1.0

