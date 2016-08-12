#!/bin/bash
#RUN FROM /home/janao/Genotipi/Genotipi1_12042016/


################################################################
#first set the set of chips used for concordance check (common SNPs)

CHIPS=GSChips #Options AllChips, GSChips, GS30KChips - merged PLINK files
cd /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS

for REF in GGP GGPv03 GP3 GP4 HD #50K  
do 
	mkdir Imp_$REF



	###############################################################
	#script to mask SNPs for 5x cross validation imputation
	#make 5 files for imputation with 1/5 missing SNPs in each - then impute separately and check concordance
	################################################################
	#first mask 1/5 SNPs each time and produce 5 bed/bim/fam files
	#All chips are merged in the GP4/OUTPUT directory
	for i in 1 2 3 4 5
	do
		MERGEDPATH=/home/janao/Genotipi/Genotipi1_12042016/MS_Imputation/Imputation_GP4ref/OUTPUT/ #set a path to the merged CHIPS file
		~/bin/plink --file PLINK_merged_exc --cow --zero-cluster /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/ClusterSNPs_CV5_$CHIPS.txt --within /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/IndClus$i.txt --make-bed --out /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/Merged_Masked$i
		~/bin/plink --bfile /home/janao/Genotipi/Genotipi1_12042016/Imp_$REF/Merged_Masked$i --recode --cow --out /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/Merged_Masked$i
	done;

	####################################################3
	#script for imputation of each of the 5 ped and map files with 1/5 missing SNPs
	#######################################################
	

	for i in 1 2 3 4 5
	do
		cd /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/
		cp /home/janao/Genotipi/Genotipi1_12042016/MS_Imputation/PARAMFILE.txt /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/
		MASKEDPATH=/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/
		REFLIST=/home/janao/Genotipi/Genotipi1_12042016/Combined_AB
		sed -i "s%PathToPed%$MASKEDPATH/Merged_Masked${i}.ped%g" PARAMFILE.txt #change ped file input 
		sed -i "s%PathToMap%$MASKEDPATH/Merged_Masked${i}.map%g" PARAMFILE.txt #change map file input
		sed -i "s%OutputName%Merged_Masked${i}%g" PARAMFILE.txt #change output name
		sed -i "s%FImputeOptions%ref=$REFLIST/${REF}_ind.txt;%g" PARAMFILE.txt #change reference for the imputation - list of individuals
		python ~/Genotipi/Zanardi/Zanardi.py --fimpute --save
		rm PARAMFILE.txt
		cd /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/
	done;

	

	##########################################################################
	#a script to extract SNPs fro 5times cross validation - from a single genotype file (merged) and 5 imputed files
	#extract a different subset of SNPs each time
	#also check concordance simultaneously
	#############################################################################
	
	for i in 1 2 3 4 5
	do
		
		cd /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/
		GENPATH=/home/janao/Genotipi/Genotipi1_12042016/MS_Imputation/Imputation_GP4ref/OUTPUT #merged and excluded (0,0 SNPs) PLINK file
		IMPPATH=/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/OUTPUT #Masked and imputed PLINK file (5 times, 5 files)
		SNPPATH=/home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/Imp_$REF/ #lists of SNPs masked in each imputation
		~/bin/plink --file $GENPATH/PLINK_merged_exc --extract $SNPPATH/CVSNPset${i}.txt --cow --recode --out Genotyped$i #extract masked SNPs from merged genotype file
		~/bin/plink --file $IMPPATH/FIMPUTE_Merged_Masked$i --extract $SNPPATH/CVSNPset${i}.txt --cow --recode --out Imputed$i #extract masked and imputed SNPs from a file
		~/bin/plink --file Genotyped$i --merge Imputed$i.ped Imputed$i.map --merge-mode 7 --cow --recode --out Concordance$i 
		cd /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/
	done;
	
	cd /home/janao/Genotipi/Genotipi1_12042016/SNP_Imputation/$CHIPS/
done;

