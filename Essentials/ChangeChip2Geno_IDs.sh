cd ./SimulatedData/AllIndividualsSnpChips/
sed 's/^ *//' Chip2Genotype.txt | cut -f1 -d" " > Individuals.txt && mv Individuals.txt Chip2Genotype.txt
sed 's/^ *//' Chip2Phase.txt | cut -f1 -d" " > Individuals.txt && mv Individuals.txt Chip2Phase.txt
