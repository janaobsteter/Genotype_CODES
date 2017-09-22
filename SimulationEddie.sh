
#!/bin/sh
# Grid Engine options (lines prefixed with #$)
#$ -N hello              
#$ -/mnt/ris-fas1a/linux_groups2/jhickey_grp/ggorjanc/Jana_Simulacije                  
#$ -l h_rt=00:05:00 
#$ -l h_vmem=16G
#$ -pe sharedmem 8
#$ -e /mnt/ris-fas1a/linux_groups2/jhickey_grp/ggorjanc/Jana_Simulacije /SimError.txt

#  These options are:
#  job name: -janaSimulacije
#  use the current working directory: -cwd
#  runtime limit of 5 minutes: -l h_rt
#  memory limit of 1 Gbyte: -l h_vmem

# Initialise the environment modules
. /etc/profile.d/modules.sh
 
# Load Python
module load python/2.7.1

 
# Run the program
./Eddie_SelectionRep.py

