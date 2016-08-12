for i in Ped*.ped; #Ped*.ped are combined PED files for all individuals genotyped on the same chip
do
FNAME=${i%.ped} #create variable for chip (e.g. Ped50K, PedGP4)
REFPATH=/home/janao/Genotipi/MS_impute_phased_Ref+Marker_files/ #path to reference files (McClure)
cd ./$FNAME #move to chip directory
	for each in 2 #run only for chromosome 2
	do 
		/home/janao/Downloads/fcgene-1.0.7/fcgene --file ${FNAME}chr${each} --oformat beagle --out ${FNAME}chr${each} #convert each chromosome PED to BGL
		java -Xmx1000m -jar ~/Downloads/beagle.jar unphased=${FNAME}chr2_sorted.bgl phased=$REFPATH/p_All_1kb+MS_BT_ref_chr${each}.txt markers=$REFPATH/MinSNP+MS_chr${each}.txt missing=? out=${FNAME}chr${each} #run each BGL with reference chromosome file
	done
cd ..
done
