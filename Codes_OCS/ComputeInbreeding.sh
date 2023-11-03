WorkingDir=/home/v1jobste/JanaO/10K_Ref_1Pb/
for scenario in "GenSLO" "OtherCowsGen" "BmGen" "Gen" "Class"
do
	for rep in $(seq 11 19)
	do
		awk -F" " '{print $2, $3, $4}' $WorkingDir/${scenario}${rep}/SimulatedData/PedigreeAndGeneticValues_cat.txt > Pedigree${scenario}${rep}.txt
		awk -F" " '{print $2, $1}' $WorkingDir${scenario}${rep}/SimulatedData/PedigreeAndGeneticValues_cat.txt > Generation${scenario}${rep}.txt
		cp AlphaRelateSpec_generic.txt AlphaRelateSpec.txt
		sed -i "s%PEDIGREEFILE%Pedigree${scenario}${rep}.txt%g" AlphaRelateSpec.txt
		sed -i "s%GENERATIONFILE%Generation${scenario}${rep}.txt%g" AlphaRelateSpec.txt
		./AlphaRelate_Linux
		mv PedigreeInbreeding.txt Inbreeding_${scenario}${rep}.txt
	done
done
	
	
