#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Pheno_Gen_REP__NAME_
#$ -cwd
#$ -l h_vmem=64G
#$ -o /home/v1jobste/JanaO/Pheno_Gen_REP___NAME_.txt
#$ -l h_rt=50:00:00
#$ -j yes
#$ -P roslin_hickey_group
#  These options are:
#  job name: -N
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem



# Initialise the environment modules
. /etc/profile.d/modules.sh

# Load Python
module load python/2.7.10

# Run the program
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/anaconda2/bin/python Eddie_SelectionRep_Selection_repeatedPheno_Male.py _REP_ Gen SU55 _REFERENCE_ _REPEATS_ 1.0,0.333333,0.333333,0.333333,1.0 _NAME_


