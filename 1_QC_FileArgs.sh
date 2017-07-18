#!/bin/bash

#first argument is the PLINK ped and map file name, second is the chip
#Edit path to script! 
#Run within the directory of the files
export S=/home/jana/Genotipi/Genotipi_CODES


## Check for genotype call, heterozygosity
plink --file $1 --missing --het --cow --recode --out $1_$2

## Extrac the individuals that have genotype call less than 0.9
tail -n +2 $1_$2.imiss | awk '$6 >= 0.10 { print $1,$2 }' > IDsWithCallRateLessThan0.90$1_$2.txt

## Draw plot for percentage of genotype call rate and heterozygosity rate for individuals
## You need geneplotter package to run next script. Install geneplotter with the following command
# source("https://bioconductor.org/biocLite.R")
# biocLite("geneplotter")
${S}/DrawHetMissPlotInd.R $1_$2
echo "Created plot:$1_$2_imiss-vs-het.pdf"

## Remove individuals with low genotype call rate and heterozygosity rate out of range +-6SD of heterozygosity
plink --file $1_$2 --remove IndividualsToExlcudeMissHet.txt --missing --recode --cow --out $1_$2_CleanInds
mv IndividualsToExlcudeMissHet.txt IndividualsToExlcudeMissHet_$1_$2.txt

## Extract SNP with missing information on more than 10% of genotypes
tail -n +2 $1_$2_CleanInds.lmiss | awk '$5 > 0.10 { print $1,$2 }' > MarkersWithCallRateLessThan0.90Plink_$1_$2.txt

## Draw plot for percentage of genotype call rate and heterozygosity rate for markers
${S}/DrawMissPlotMarkers.R $1_$2_CleanInds
echo "Created plot:$1_$2_CleanInds_lmiss.pdf"

## Remove SNP with missing information on more than 10% of genotypes
plink --file $1_$2_CleanInds --recode --geno 0.10 --cow --out $1_$2_CleanIndsMarkers