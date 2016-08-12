#!/bin/bash
NUMIT=1000 # nunmber of iterations

for i in 1 2 3 5 9 15 16 18 19 20 21 23
do
read NUM_SI < <(head -n1 Chr${i}_${NUMIT}it.Chr${i}ImpBase.bgl.phased |  sed 's/ /\n/g' | grep ^SI[0-9] | wc -l)
NUMSII=$(($NUM_SI + 2))
cut -d " " -f-$NUMSII Chr${i}_${NUMIT}it.Chr${i}ImpBase.bgl.phased > SICHR${i}.bgl
done


