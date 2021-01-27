#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N NameOfTheJob
#$ -cwd
#$ -l h_vmem=4G
#$ -pe sharedmem 1
#$ -o NameOfTheOutputFile
#$ -l h_rt=04:00:00
#$ -j yes
#  These options are:
#  job name: -These lines are commented out
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem


# Initialise the environment modules
. /etc/profile.d/modules.sh

module load LoadTheModulesYouWant (you can check the modules with module avail on Eddie)

# Run the program
Here you put the line to run your program (Rscript or something else)
