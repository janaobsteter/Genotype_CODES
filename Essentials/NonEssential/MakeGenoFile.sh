sed 's/^ *//' ./SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt > ChipFile.txt
cut -f1 -d ' ' ChipFile.txt > Individuals.txt
awk '{$1=""; print $0}' ChipFile.txt | sed 's/ //g' > Snps.txt
paste Individuals.txt Snps.txt | awk '{printf "%- 10s %+ 15s\n",$1,$2}' > GenoFile.txt

