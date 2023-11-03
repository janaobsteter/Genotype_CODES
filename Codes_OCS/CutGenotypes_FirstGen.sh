for rep in $(seq 0 19)
do
	cd /home/v1jobste/JanaO/FillInBurnIn${rep}/SimulatedData/AllIndividualsSnpChips/
	head -n8640 Chip1Genotype.txt > FirstGen.txt
done
