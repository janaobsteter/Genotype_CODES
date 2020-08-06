#!/bin/bash

####################################################
#SCRIPT TO RUN IMPUTATION FROM SMALLER TO LARGER CHIP
#RUN FROM STEPX DIRECTORY
####################################################
#all chips with be imputed but the concordance check will be for the SCHIP
#LCHIP is the chip you want to impute to
CONCCHIP=Versa50K #<-- CHANGE!
REFCHIP=IDBv03 #<-- CHANGE!
IMPDATE=13072020
ImpDir=/home/jana/Genotipi/Genotipi_WORK/SNPImp_$IMPDATE
mkdir $ImpDir
ImpCipDir=/home/jana/Genotipi/Genotipi_WORK/SNPImp_$IMPDATE/Ref${REFCHIP}_Conc${CONCCHIP}
mkdir $ImpCipDir
cd $ImpCipDir


ADDCHIP1=GGPv02
ADDCHIP2=GGPv03
ADDCHIP3=GGPv04
ADDCHIP4=HD
ADDCHIP5=HDv02
ADDCHIP6=50Kv01
ADDCHIP7=50Kv02

#REFERENCE CHIP = chip you are imputing to
#CONCPLINK is the chip you are checking concordance for (each one has a different number of missing SNPs)
REFPLINK=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$REFCHIP/PLINK_MERGED
CONCPLINK=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$CONCCHIP/PLINK_MERGED


#
ADDPLINK1=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP1/PLINK_MERGED 
ADDPLINK2=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP2/PLINK_MERGED
ADDPLINK3=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP3/PLINK_MERGED
ADDPLINK4=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP4/PLINK_MERGED
ADDPLINK5=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP5/PLINK_MERGED
ADDPLINK6=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP5/PLINK_MERGED
ADDPLINK7=/home/jana/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava/Top/$ADDCHIP5/PLINK_MERGED

#remove potential Concordance files left behind
rm Conc*

#remove duplicated SNPs from REFPLINK file
~/bin/plink --file $REFPLINK --cow --list-duplicate-vars ids-only --out REFDUPL
~/bin/plink --file $REFPLINK --cow --exclude REFDUPL.dupvar --recode --out $REFCHIP

#keep only the REFPLINK SNPs in CONCFILE and other PLINK files
~/bin/plink --file ${CONCPLINK} --cow --extract $REFCHIP.map --chr 1-29 --recode --out CONCPLINK_REF
~/bin/plink --file ${ADDPLINK1} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP1}_REF
~/bin/plink --file ${ADDPLINK2} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP2}_REF
~/bin/plink --file ${ADDPLINK3} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP3}_REF
~/bin/plink --file ${ADDPLINK4} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP4}_REF
~/bin/plink --file ${ADDPLINK5} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP5}_REF
~/bin/plink --file ${ADDPLINK5} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP6}_REF
~/bin/plink --file ${ADDPLINK5} --cow --extract $REFCHIP.map --chr 1-29 --recode --out ${ADDCHIP7}_REF


#prepare SNP_Cluster file - a list of common SNP on the REF and CONCHIP
#Split the list into 10x datasets - CVSNPset#
cp /home/jana/Genotipi/Genotipi_CODES/Cluster_SNPs.R .
sed -i "s/CONCMAP_REF/CONCPLINK_REF.map/g" Cluster_SNPs.R
sed -i "s%CurrentImputationDir%$PWD%g" Cluster_SNPs.R

#!/usr/bin/enc Rscript
Rscript $PWD/Cluster_SNPs.R $ImpCipDir


#SNPSETS=$PWD/CVSNPset
#CLUSTERFILE=/home/jana/Genotipi/Genotipi1_12042016/StepImputation/Step$STEP/SNP_Cluster.txt
#NDCLUSTER=/home/jana/Genotipi/Genotipi1_12042016/StepImputation/Step$STEP

#create a merge-list txt file for PLINK to merge all the files onto REFCHIP
echo "CONCPLINK_REF.ped" "CONCPLINK_REF.map"$'\n'"${ADDCHIP1}_REF.ped" "${ADDCHIP1}_REF.map"$'\n'"${ADDCHIP2}_REF.ped" "${ADDCHIP2}_REF.map"$'\n'"${ADDCHIP3}_REF.ped" "${ADDCHIP3}_REF.map"$'\n'"${ADDCHIP4}_REF.ped" "${ADDCHIP4}_REF.map"$'\n'"${ADDCHIP5}_REF.ped" "${ADDCHIP5}_REF.map" > MERGELIST.txt


#merge all PLINK files
~/bin/plink --file  $REFPLINK --cow --merge-list MERGELIST.txt --chr 1-29 --recode --out MERGEDOnto_$REFCHIP

#Remove SNPs with different positions on different chromosomes
grep Warning: MERGEDOnto_$REFCHIP.log | grep -o "'[^']\+'" | sed "s/'//g" > RemoveSNPs.txt



#Extract individuals FID and ID from the genotype files - in order to create a cluster individuals file
cut -d " " -f1,2 CONCPLINK_REF.ped > CONC_ind.txt
cut -d " " -f1,2 ${REFCHIP}.ped > REF_ind.txt


#mask every tenth SNP 
#first mask 1/10 SNPs each time and produce 10 bed/bim/fam files
for i in 1 2 3 4 5 6 7 8 9 10
do

#create individual cluster file with all the inviduals on the CONC chip masked (assigned to a cluster 1-10) - all individuals on the chip masked 1/10 small-common SNPs - 10times
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
	' CONC_ind.txt > CONC_ClustInd$i.txt


	~/bin/plink --file CONCPLINK_REF --cow --exclude RemoveSNPs.txt --zero-cluster SNP_Cluster.txt --within CONC_ClustInd$i.txt --make-bed --out CONC_Masked$i
	~/bin/plink --bfile CONC_Masked$i --recode --cow --out CONC_Masked$i



####################################################3
#Impute all 10 Masked files
#######################################################
#REFCHIP is the REFPLINK file with 0,0 SNPs and sex chromosomes excluded
	echo "This is the working directory"
	echo $PWD
	#ls $PWD
	cp -r ~/Genotipi/Genotipi_CODES/Zanardi/* .
        rm -r OUTPUT Example_data
	cp /home/jana/Genotipi/Genotipi_CODES/PARAMFILE.txt .
	sed -i "s%PathToPed%$PWD/CONC_Masked$i.ped,$PWD/$REFCHIP.ped%g" PARAMFILE.txt #change ped file input 
	sed -i "s%PathToMap%$PWD/CONC_Masked$i.map,$PWD/$REFCHIP.map%g" PARAMFILE.txt #change map file input
	sed -i "s%OutputName%ImpMasked${i}%g" PARAMFILE.txt #change output name
	cp PARAMFILE.txt PARAMFILE_Current.txt
	python Zanardi.py --fimpute --save
#	python ~/Genotipi/Genotipi_CODES/Zanardi/Zanardi.py --save
	
	rm PARAMFILE.txt
#########################################################################
#a script to extract SNPs from 10 imputed MergeMasked files and one original genotype file
#extract a different subset of SNPs each time
#also check concordance simultaneously
#############################################################################

	~/bin/plink --file CONCPLINK_REF --extract CVSNPset${i}.txt --cow --recode --out Genotyped$i #extract masked SNPs from merged genotype file
	~/bin/plink --file FIMPUTE_ImpMasked$i --extract CVSNPset${i}.txt --cow --recode --out Imputed$i #extract masked and imputed SNPs from a file
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
#cp /home/jana/Genotipi/Genotipi_CODES/PARAMFILE.txt .
#sed -i "s%PathToPed%$PWD/${SCHIP}_exc.ped,$PWD/${LCHIP}.ped%g" PARAMFILE.txt #change ped file input 
#sed -i "s%PathToMap%$PWD/${SCHIP}_exc.map,$PWD/${LCHIP}.map%g" PARAMFILE.txt #change map file input
#sed -i "s%OutputName%Imputed%g" PARAMFILE.txt #change output name
#python ~/Genotipi/Zanardi/Zanardi.py --fimpute --save


