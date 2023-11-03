#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N ImportFillIn__REP_
#$ -cwd       
#$ -l h_vmem=16G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/HighlanderLab//jobsteter/ImportFillIn__REP_.txt
#$ -l h_rt=24:00:00
#$ -j yes
#$ -P roslin_HighlanderLab


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
# Load Python
module load anaconda
source activate jana_python2.7


# Run the program
python2.7 Eddie_TwoPopulations_FillIN.py _REP_ 
