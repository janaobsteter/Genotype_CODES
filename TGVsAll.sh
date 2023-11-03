#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Jana_OCSAlphaMate
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o /exports/cmvm/eddie/eb/groups/tier2_hickey_external/JanaO/OCS_AlphaMate.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_hickey_group
#  These options are:
#  job name: -OCS_AlphaMate
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
#load R module
module load R/3.3.2

R --file=/home/v1jobste/JanaO/TGVsAll.R > /home/v1jobste/JanaO/TGVsAll_AlphaMate.txt

