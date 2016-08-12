#!/bin/bash

####################################################
#SCRIPT TO RUN IMPUTATION FROM SMALLER TO LARGER CHIP
####################################################
CHIP1=GGP
CHIP2=GGPv03
EXCLUDE=/home/janao/Genotipi/Genotipi1_12042016/StepImputation/Exclude_from_GGP.txt #GGP exclusive SNPs
SMALLPLINK=~/Genotipi/Genotipi1_12042016/PLINK_genotypeFiles/$CHIP1/OUTPUT/PLINK_MERGED
LARGEPLINK=~/Genotipi/Genotipi1_12042016/PLINK_genotypeFiles/$CHIP2/OUTPUT/PLINK_MERGED
LARGEPLINK2=/home/janao/Genotipi/Genotipi1_12042016/PLINK_genotypeFiles/GP3v02/OUTPUT/PLINK_MERGED #this is for 26K chip, because it has two chip names
SNPSETS=/home/janao/Genotipi/Genotipi1_12042016/StepImputation/19K26K_CVSNPset
CLUSTERFILE=/home/janao/Genotipi/Genotipi1_12042016/StepImputation/19K26K_Cluster.txt
INDCLUSTER=/home/janao/Genotipi/Genotipi1_12042016/StepImputation/

#Rscript 
echo "Did you run the R script to obtain exclude and cluster files?"
echo "Did you change the input file variables?"
echo "Did you comment out metge two large plink files?"

#first exclude the exclusive GGP SNPs from GGP file
~/bin/plink --file $SMALLPLINK --cow --exclude $EXCLUDE  --chr 1-29 --recode --out 19K_excF

#merge GP3 and GGP3 maps together
~/bin/plink --file  $LARGEPLINK --cow --merge $LARGEPLINK2.ped $LARGEPLINK2.map --chr 1-29 --recode --out 26KF

#create a file with variants with the same physical position
cut -d ";" ./OUTPUT/Error_SNP_position.txt -f1 | tail -n +2 > SamePosSNP.txt

#Exclude these SNPs from PLINK files
~/bin/plink --file 19K_excF --cow --exclude SamePosSNP.txt --recode --out 19K_exc
~/bin/plink --file 26KF --cow --exclude SamePosSNP.txt --recode --out 26K

#Extract individuals FID and ID from the genotype files - in order to create a cluster individuals file
cut -d " " -f1,2 19K_exc.ped > 19K_ind.txt
cut -d " " -f1,2 26K.ped > 26K_ind.txt


#mask every tenth SNP 
#first mask 1/10 SNPs each time and produce 10 bed/bim/fam files
for i in 1 2 3 4 5 6 7 8 9 10
do

#create individual cluster file with all the inviduals on the lower chip masked (assigned to a cluster 1-10)
#awk script to add Cluster CX number to the third column
	awk -v column=3 -v value="C$i" '
	    BEGIN {
	        FS = OFS = " ";
	    }
	    {
		for ( i = NF + 1; i > column; i-- ) {
		    $i = $(i-1);
		}
		$i = value;
		print $0;
	    }
	' 19K_ind.txt > 19K_ClustInd$i.txt


	~/bin/plink --file 19K_exc --cow --zero-cluster $CLUSTERFILE --within 19K_ClustInd$i.txt --make-bed --out 19K_Masked$i
	~/bin/plink --bfile 19K_Masked$i --recode --cow --out 19K_Masked$i


#############################################################
#when ran for the first time you get an error - a list of duplicated SNPs from FImpute OUTPUT
#remove them using PLINK -when u run for the second time
############################################################################
	~/bin/plink --bfile 19K_Masked$i --recode --cow --out 19K_Masked$i


####################################################3
#Impute all 10 MergMasked files
#######################################################

	cp /home/janao/Genotipi/Genotipi1_12042016/PARAMFILE.txt .
	sed -i "s%PathToPed%$PWD/19K_Masked${i}.ped,$PWD/26K.ped%g" PARAMFILE.txt #change ped file input 
	sed -i "s%PathToMap%$PWD/19K_Masked${i}.map,$PWD/26K.map%g" PARAMFILE.txt #change map file input
	sed -i "s%OutputName%ImpMasked${i}%g" PARAMFILE.txt #change output name
	python ~/Genotipi/Zanardi/Zanardi.py --fimpute --save
	rm PARAMFILE.txt
#########################################################################
#a script to extract SNPs from 10 imputed MergeMasked files and one original genotype file
#extract a different subset of SNPs each time
#also check concordance simultaneously
#############################################################################

	~/bin/plink --file 19K_exc --extract 19K26K_${SNPSETS}${i}.txt --cow --recode --out Genotyped$i #extract masked SNPs from merged genotype file
	~/bin/plink --file OUTPUT/FIMPUTE_ImpMasked$i --extract ${SNPSETS}${i}.txt --cow --recode --out Imputed$i #extract masked and imputed SNPs from a file
	~/bin/plink --file Genotyped$i --merge Imputed$i.ped Imputed$i.map --merge-mode 7 --cow --recode --out Concordance$i 
	grep "for a concordance rate" Concordance$i.log > Conc_num$i.txt
done;

cat Conc_num* > Concordance_num.txt
grep -o ' 0.*' Concordance_num.txt > Conc_num.txt
sed -i 's/\.$//' Conc_num.txt
awk '{ sum += $1 } END { if (NR > 0) print sum / NR }' Conc_num.txt > Concordance_avg.txt

##############################################
#impute non-masked files to obtain imputated genotypes for all the SNPs
###############################################
cp /home/janao/Genotipi/Genotipi1_12042016/PARAMFILE.txt .
sed -i "s%PathToPed%$PWD/${SCHIP}_exc.ped,$PWD/${LCHIP}.ped%g" PARAMFILE.txt #change ped file input 
sed -i "s%PathToMap%$PWD/${SCHIP}_exc.map,$PWD/${LCHIP}.map%g" PARAMFILE.txt #change map file input
sed -i "s%OutputName%Imputed%g" PARAMFILE.txt #change output name
#python ~/Genotipi/Zanardi/Zanardi.py --fimpute --save


