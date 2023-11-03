#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Import_SCENARIOHOME__SCENARIOIMPORT___PERIMPORTK__PERIMPORTBM__TRAITHOME__TRAITIMPORT___REP_
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/Import_SCENARIOHOME__SCENARIOIMPORT___PERIMPORTK__PERIMPORTBM__TRAITHOME__TRAITIMPORT___REP_.txt
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
/exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/anaconda2/bin/python Eddie_SelectionRep_Selection_import.py _REP_ "_SCENARIOHOME_" "_SCENARIOIMPORT_" "_STRATEGY_" "10K" _TRAITHOME_ _TRAITIMPORT_  _PERIMPORTK_ _PERIMPORTBM_
