################################################################################
#script to convert chromosome PED files to BGL files 
###############################################################################3



for each in 1 2 3 5 9 15 16 18 19 20 21 23 #required chromosomes with MS
do 
	/home/janao/Downloads/fcgene-1.0.7/fcgene --file Imputed_chr${each} --oformat beagle --out Chr${each} #convert each chromosome PED to BGL
		
done

