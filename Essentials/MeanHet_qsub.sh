#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N MeanHet:_AllGen
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o MeanHet.txt
#$ -l h_rt=01:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -janaSimulacije
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
module load igmm/apps/R/3.4.1
 
# Run the program
Rscript MeanHetMarker_Neutral_QTN_import_AllGen.R 0 "Class" "Class" 1 3 "SU55"
