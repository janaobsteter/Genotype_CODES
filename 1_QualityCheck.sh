#!/bin/bash

export G=~/Dropbox/Sluzba_2017/F4F/Genotipi/Genotipi_latest/Rjava/IDBv03
export W=~/Dropbox/Sluzba_2017/F4F/Work
export S=~/Dropbox/Sluzba_2017/F4F/Scripts
export P=~/Dropbox/Sluzba_2017/F4F/Programs

cd ${W}

## Check for genotype call, heterozygosity
${P}/plink_mac/plink --file ${G}/F4F_11072017_Top --missing --het --cow --recode --out F4F_11072017_Top

## Extrac the individuals that have genotype call less than 0.9
tail -n +2 F4F_11072017_Top.imiss | awk '$6 >= 0.10 { print $1,$2 }' > IDsWithCallRateLessThan0.90F4F_11072017_Top.txt

## Draw plot for percentage of genotype call rate and heterozygosity rate for individuals
## You need geneplotter package to run next script. Install geneplotter with the following command
# source("https://bioconductor.org/biocLite.R")
# biocLite("geneplotter")
${S}/DrawHetMissPlotInd.R F4F_11072017_Top

## Remove individuals with low genotype call rate and heterozygosity rate out of range +-6SD of heterozygosity
${P}/plink_mac/plink --file F4F_11072017_Top --remove IndividualsToExlcudeMissHet.txt --missing --recode --cow --out F4F_11072017_Top_CleanInds
mv IndividualsToExlcudeMissHet.txt IndividualsToExlcudeMissHet_F4F_11072017_Top.txt

## Extract SNP with missing information on more than 10% of genotypes
tail -n +2 F4F_11072017_Top_CleanInds.lmiss | awk '$5 > 0.10 { print $1,$2 }' > MarkersWithCallRateLessThan0.90Plink_F4F_11072017_Top.txt

## Draw plot for percentage of genotype call rate and heterozygosity rate for markers
${S}/DrawMissPlotMarkers.R F4F_11072017_Top_CleanInds

## Remove SNP with missing information on more than 10% of genotypes
${P}/plink_mac/plink --file F4F_11072017_Top_CleanInds --recode --geno 0.10 --cow --out F4F_11072017_Top_CleanIndsMarkers

