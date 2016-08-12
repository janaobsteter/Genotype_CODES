################################################################################
#script to convert combined imputed chromosome PED files to BGL files and simultaneously impute MS
#run from Imputation_GP4ref/OUTPUT
###############################################################################3

REFPATH=/home/janao/Genotipi/MS_impute_phased_Ref+Marker_files/ #path to reference files (McClure)
for each in 1 2 3 5 9 15 16 18 19 20 21 23 #required chromosomes with MS
	do 
		/home/janao/Downloads/fcgene-1.0.7/fcgene --file GP4Imputed_chr${each} --oformat beagle --out Chr${each} #convert each chromosome PED to BGL
		java -Xmx1000m -jar ~/Downloads/beagle3.jar unphased=Chr${each}.bgl phased=$REFPATH/p_All_1kb+MS_BT_ref_chr${each}.txt markers=$REFPATH/MinSNP+MS_chr${each}.txt missing=? out=ImputedMS_chr${each} #run each BGL with reference chromosome file
done


