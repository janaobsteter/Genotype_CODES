#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N Compress_STRATEGY__SCENARIO__TRAIT__REP_
#$ -cwd       
#$ -l h_vmem=32G
#$ -pe sharedmem 1
#$ -o Compress_STRATEGY__SCENARIO__TRAIT__REP_.txt
#$ -l h_rt=48:00:00
#$ -j yes
# -P roslin_HighlanderLab
#  These options are:
#  job name: -Cut_GENOTYPES
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh
 
tar -zcvf _STRATEGY__REP__SCENARIO__TRAIT_.tar.gz /home/v1jobste/jobsteter/10K/SU55_import/_STRATEGY__REP___SCENARIO_1_TRAIT_
rm -rf /home/v1jobste/jobsteter/10K/SU55_import/_STRATEGY__REP___SCENARIO_1_TRAIT_

