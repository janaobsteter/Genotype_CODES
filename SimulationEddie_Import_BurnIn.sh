#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N ImportBurnIn__REP_
#$ -cwd      
#$ -l h_vmem=50G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab/jobsteter/ImportBurnIn__REP_.txt
#$ -l h_rt=48:00:00
#$ -j yes
#$ -P roslin_HighlanderLab

# Initialise the environment modules
. /etc/profile.d/modules.sh
 
# Load Python
module load anaconda
source activate jana_python2.7


# Run the program
python2.7 -u Eddie_TwoPopulations_BurnINNoFillIn.py _REP_ "Class" "SU55" "10K" traitHome traitImport
