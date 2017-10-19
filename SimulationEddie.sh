#!/bin/sh

#$ -N Jana_SimulationRep
#$ -cwd
#$ -l h_vmem=32G
#$ -e /mnt/ris-fas1a/linux_groups2/jhickey_grp/ggorjanc/Jana_Simulacije/Simulacije_error.txt


# Initialise the module environment
source /etc/profile.d/modules.sh

# Load Python
module load python/2.7.10

# Run the program
python Eddie_SelectionRep.py

