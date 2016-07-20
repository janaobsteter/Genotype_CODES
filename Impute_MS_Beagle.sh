################################################################################
#script to convert combined imputed chromosome PED files to BGL files and simultaneously impute MS
#run from Imputation_GP4ref/OUTPUT
###############################################################################3

REFBGLPATH=/home/janao/Genotipi/MS_Imputation/MS_03062016/AddingToRef #path to reference files (McClure)
REFMARKERPATH=/home/janao/Genotipi/MS_Imputation/beagle
SPATH=/home/janao/Genotipi/MS_Imputation/MS_03062016
for each in 1 2 3 5 9 15 16 18 19 20 21 23 #required chromosomes with MS
	do 
		java -Xmx1000m -jar ~/bin/beagle3.jar unphased=$SPATH/Chr${each}.bgl phased=$REFBGLPATH/Chr${each}RefAdded.bgl markers=$REFMARKERPATH/MinSNP+MS_06112013_chr9${each}.txt missing=? niterations=10 out=AddedRef_10it #run each BGL with reference chromosome file
done


