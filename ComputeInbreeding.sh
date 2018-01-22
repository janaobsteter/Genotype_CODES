WorkingDir=/home/jana/bin/AlphaSim1.05Linux/REAL20GenSel_
for scenario in "Class1" "GenSLO" "OtherCowsGen" "BmGen" "Gen"
do
	awk -F" " '{print $2, $3, $4}' $WorkingDir${scenario}/SimulatedData/PedigreeAndGeneticValues_cat.txt > Pedigree${scenario}.txt
	awk -F" " '{print $2, $1}' $WorkingDir${scenario}/SimulatedData/PedigreeAndGeneticValues_cat.txt > Generation${scenario}.txt
	cp AlphaRelateSpec_generic.txt AlphaRelateSpec.txt
	sed -i "s%PEDIGREEFILE%Pedigree${scenario}.txt%g" AlphaRelateSpec.txt
	sed -i "s%GENERATIONFILE%Generation${scenario}.txt%g" AlphaRelateSpec.txt
	./AlphaRelate_Linux
	mv PedigreeInbreeding.txt Inbreeding_${scenario}.txt
	
done
	
	
