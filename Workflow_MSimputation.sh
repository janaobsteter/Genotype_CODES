#!/bin/bash

#1) first transform all the Raw GeneSeek files to PLINK files (use the option of top/bot, A/B, F/R allele)
#SNPchimpRepo PEDDA_ROW python script
#peddar.param parameter file
bash /home/janao/Genotipi/Raw_to_Ped.sh

#2) create chip directories and move PLINK ped and map files to the corresponding directories - manually
#name chip directories 50K GGP GGPv03 GP3v02 GP4 HD

#3) Run python script Zanardi to combine all ped and map 
#parameter file PARAMFILE.txt
bash /home/janao/Genotipi/Combine_Zanardi.sh

#4) Exclude SNPs with chr0 pos0 
#list of all - Illumina + GeneSeek - 
#/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/Exclude_SNPs.txt
sudo ~/bin/plink --file PLINK_MERGED --recode --cow --exclude /home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/Exclude_SNPs.txt --out PLINK_merged_exc

#5) use Zanardi python script to impute using GP4 reference
#INPUT_PED=/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT/PLINK_merged_exc.ped
#INPUT_MAP=/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips/OUTPUT/PLINK_merged_exc.map
#PARAMFILE in /home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/GeneSeek_chips
#use PLINK_merged_exc for imputation
#outname=GP4_GSImputed
sudo python ~/Genotipi/Zanardi/Zanardi.py --fimpute --save

#6) split imputed file by chr
#input file is outfile from previous step (GP4_GSImputed)
sudo sh /home/janao/Genotipi/Split_by_chr.sh

#7) use fcgene and beagle3 to impute MS from IMPUTED data
sudo ~/Genotipi/Impute_MS_Beagle.sh

#8) Create directory called
mkdir Beagle_imputedMS
mv chr* Beagle_imputedMS
mv Chr* Beagle_imputedMS

#9) problem with chr2 - markers not in the right order
#therefore go and merge Chr2.bgl file and MinSetCHR2 file, order by position and resave in R
#!/bin/R
#change read in directory for the Chr2.bgl file
R 
~/Genotipi/Order_markers_chr2.R
#write the file then go and add first two lines from original Chr2.bgl file to the new Chr2_sorted.bgl

#9) Re-impute Chr2, run from OUTPUT directory
sudo sh ~/Genotipi/Impute_chr2.sh
mv chr* Beagle_imputedMS
mv Chr* Beagle_imputedMS
gunzip *
mv Chr2.Chr2_sorted.bgl.phase chr2.Chr2.bgl.phased

#10) extract microsatelites
#run from Beagle_imputedMS directory
sudo sh ~/Genotipi/extract_MSa.sh

#11) read all the files into R and calculte concordance rate
~/Genotipi/MS_concordance.R



