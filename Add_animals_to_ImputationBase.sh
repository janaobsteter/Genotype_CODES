#!/bin/bash

DIRS=/home/janao/Genotipi/Genotipi1_12042016/Imputation_HDref/OUTPUT/Beagle_imputedMS
DIRB=/home/janao/Genotipi/beagle
for i in 1 2 3 5 9 15 16 18 19 20 21 23
do

#sort the sample beagle and base beagle file
#sort -k2,2 -n $DIRS/Chr$i.bgl > Chr${i}_sorted.bgl
#sort -k2,2 -n $DIRB/MS_imputation_base_chr$i.bgl > Chr${i}Base.bgl
#capitalise both files
sed 's/.*/\U&/' $DIRS/Chr$i.bgl > Chr${i}_sorted.bgl
sed 's/.*/\U&/' $DIRB/MS_imputation_base_chr$i.bgl > Chr${i}Base.bgl
#extract Min CHR SNPs from sample file
awk 'FNR==NR{a[$2];next}($2 in a){print}' Chr${i}Base.bgl Chr${i}_sorted.bgl > Chr${i}Imp.bgl
#remove M and markername columns from Base File
cut -f3- -d " " Chr${i}Base.bgl > Chr${i}Base_.bgl
#paste together Sample bgl and Base bgl files
paste Chr${i}Imp.bgl Chr${i}Base_.bgl > Chr${i}ImpBase.bgl
done




