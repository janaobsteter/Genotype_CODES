#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N PermEnv_REP__SCENARIO__PHENOTYPE_
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/PermEnv__REP__SCENARIO__PHENOTYPE_.txt
#$ -l h_rt=120:00:00
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
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python Eddie_SelectionRep_Selection_repeatedPheno.py _REP_ _SCENARIO_ "SU55" "10K" _PHENOTYPE_ 0.5,0.333333,0.333333,0.333333,1.5
