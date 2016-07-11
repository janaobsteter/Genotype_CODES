################################################################################
#script to convert combined imputed chromosome PED files to BGL files and simultaneously impute MS
#run from Imputation_GP4ref/OUTPUT
###############################################################################3

REFPATH=/home/janao/Genotipi/MS_impute_phased_Ref+Marker_files/ #path to reference files (McClure)
java -Xmx1000m -jar ~/Downloads/beagle3.jar unphased=./Beagle_imputedMS/Chr2_sorted.bgl phased=$REFPATH/p_All_1kb+MS_BT_ref_chr2.txt markers=$REFPATH/MinSNP+MS_chr2.txt missing=? out=Chr2 
