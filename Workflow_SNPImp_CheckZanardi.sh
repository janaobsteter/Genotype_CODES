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
#r#m Conc*



#mask every tenth SNP 
#first mask 1/10 SNPs each time and produce 10 bed/bim/fam files

#create individual cluster file with all the inviduals on the CONC chip masked (assigned to a cluster 1-10) - all individuals on the chip masked 1/10 small-common SNPs - 10times
#awk script to add Cluster CX number to the third column



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
	
#	rm PARAMFILE.txt
#########################################################################
#a script to extract SNPs from 10 imputed MergeMasked files and one original genotype file
#extract a different subset of SNPs each time
#also check concordance simultaneously
#############################################################################





