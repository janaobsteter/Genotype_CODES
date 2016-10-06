#####################################################################################################
#script for adding newly downlaoded genotypes on one or multiple chips or one and multiple genotype "packages"
#suitable for GENESEEK downloads (Rjave) 
#scripts requires PLINK 1.9, SNPchimp and Zanardi software
#script creates temporary directory Genotipi_DATEDOWNLOAD with newly downloaded genotypes
#script adds new genotypes to the corresponding chip directories in Genotipi_latest and PLINK_MERGED files
#script also created a table of newly genotyped ind for Govedo GENOTIPIZIRANE_ZIVALI
#it also extracts the 800 SNP for parental verification and prepares a format for PARENTAL_SNP800
#it also utilises R scripts NewIndForGovedo.R and Add_alleles_Govedo_PARENTAL_SNP800.R
####################################################################################

#set paths
#1) path to the SNPchimpRepo
PATH_PR=/home/janao/Genotipi/Genotipi_CODES/SNPchimpRepo/source_codes/PEDDA_ROW 
#2)path to the newly downloaded genotype files
DOWNLOAD_PATH=~/Downloads
#3) path to the genotype directory - where you want to create you temp directory 
GEN_PATH=~/Genotipi/Genotipi_DATA/
#4) path to the directory where you store your latest and final genotypes
GENLAT_PATH=~/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava
#5) path to the directory where you store info and files for parental verification (800 SNPs) - i.e. names of the required SNPs and previous SNP800 files
PARSNP_PATH=~/Genotipi/ParentalVerification_SNPSNP
#6) path to Zanardi
ZAN_PATH=/home/janao/Genotipi/Genotipi_CODES/Zanardi
#7) make sure your R scripts are in the GEN_PATH directory
#8) make sure your PARAMFILE.txt is in the GEN_PATH directory
#9) path to sifrant for 800 SNPs
SIF_FILE='"~/Genotipi/ParentalVerification_SNPSNP/Sifrant_SNP.csv"'
#10) file that contains IDs and sekvence of animals
IDSEQ_FILE='"~/Genotipi/Genotipi_CODES/Rjave_seq_ID.csv"'
#11)Directory to CODE directory
CODE_DIR=~/Genotipi/Genotipi_CODES


#REQUIRED INPUT
#|||||||||||||||||||||||||||||||||||||||||||||||||
DATEDOWNLOAD=05092016 #date of downloading the file
FILE=Matija_Rigler_BOVUHDV03_20160831.zip # downloaded file
#|||||||||||||||||||||||||||||||||||||||||||||||||


mkdir $GEN_PATH/Genotipi_${DATEDOWNLOAD} #create temp directory for new genotypes
TEMPDIR=$GEN_PATH/Genotipi_${DATEDOWNLOAD}



cp $DOWNLOAD_PATH/$FILE $TEMPDIR #copy downloaded genotyped to the temp directory
cd $TEMPDIR
unzip * #unzip - may contain multiple folders - packages

#create associative array to hold the number of SNPs on each chip
declare -A chip
chip[19720]=GGPv02
chip[26151]=GGPv03
chip[30105]=GGPv04
chip[30106]=GGPv04
chip[76883]=HD
chip[138892]=HDv02
chip[139376]=HDv02
chip[54001]=50Kv01
chip[54609]=50Kv02


##################################################################################################3
#ask the user whether he/she wants to merge this newly obtained genotypes with other files in genotipi Latest

####################################################################################################
#loop through subdirectories - downloaded genotype packages
####################################################################################################

for i in */ #for all subdirectories - i.e. all the genotype "packages" downloaded
do
	cd $i
	##unzip SNP_Map and FinalReport which are needed for raw to ped tranformation
	unzip SNP_Map.zip
	unzip *FinalReport.zip
	##################################################################################
	#PEDDAROW Raw to Ped
	#script for editing parameter file for PEDDAROW
	#PEDDA_ROW = python program for raw Illumina --> ped in map
	##############################################################################
	PASMA="Rjava"
	cp $PATH_PR/pedda_row.py . #copy python script to each subdirectory
	cp $PATH_PR/peddar.param . #coppy parameter file in each subdirectory
	_FINREP=$(find *FinalReport.txt -printf '%f\n') #find FinalReport
	SERNUM=${_FINREP%_FinalReport.txt} #remove sufix to only get the number of the genotyping package
	#change parameter file
	sed -i "s%SI %SI%g" ${_FINREP} #remove blank space from SI #####
	sed -i "s%SI %SI%g" ${_FINREP} #remove blank space from SI #####
	sed -i "s%SI %SI%g" ${_FINREP} #remove blank space from SI #####
	sed -i "s%SI %SI%g" ${_FINREP} #remove blank space from SI #####
	sed -i "s/test_FinalReport.txt/"${_FINREP}"/g" peddar.param #change finrep parameter in place
	sed -i "s/test_outputfile/"${SERNUM}"/g" peddar.param #change output prefix parameter in place
	sed -i "s/TEST/$PASMA/g" peddar.param #change brdcode parameter in place
	python pedda_row.py
	
	################################################################################
	#check whether there are any animal names in the ped files --> if there are, notify and remove them from ped file
	###################################################################################
	read LINES < <(less *.ped | cut -f3 -d " " | grep [A-Z] | wc -l)
	if [ $LINES -gt 0 ]
		then
		echo " "
		echo "WARNING!!!"
		echo "Ped file in ${i} contains animal names in the third column"
		NAMES=$(less *.ped | cut -f3 -d " " | grep [A-Z])
		for n in $NAMES
		do 
			sed -i "s/${n} //g" *.ped
		done
	fi
	
	echo "No animals names found in the ped file, proceding with the script"
	#rm $DOWNLOAD_PATH/$FILE #everything is OK, you can remove the ziped downloaded file

	############################################################################
	#read the number of SNPs on the chip and copy to the right chip directory within Genotipi/Genotipi_latest
	############################################################################
	read NUMSNP < <( less *map | wc -l)
	CH_IP=$(echo ${chip[$NUMSNP]})
	echo "The chip is $CH_IP"
	##the number of genotyped animals in the package
	read ANINUM < <( less *ped | wc -l)
	echo "Number of genotypes animals is $ANINUM"
	## copy the ped and map files to the chip directory in the genotipi latest
	cp *ped $GENLAT_PATH/$CH_IP
	cp *map $GENLAT_PATH/$CH_IP
	echo "Copying ped and map files from $TEMPDIR/$i to ~/Genotipi_latest/Rjava/$CH_IP"
	
	#################################################################################
	#get the list of newly obtained ped and map files for Zanardi PARAMFILE edit
	#obtain them in this directory since in the Genotipi_latest there are all the ped and map file
	#################################################################################
	ls -d $PWD/*.ped > ${CH_IP}_ped.txt #list all new ped and map files - for Zanardi update of the PLINK_MERGED file - 		just to add the new files
	ls -d $PWD/*.map > ${CH_IP}_map.txt
	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${CH_IP}_ped.txt #replace new lines with commas
	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${CH_IP}_map.txt
	pedfiles=$(cat ${CH_IP}_ped.txt)
	mapfiles=$(cat ${CH_IP}_map.txt)

	######################################################################################
	#extract the names of the genotyped individuals and transform them into a format for the Govedo database
	##############################################################################
	ped=$(ls *ped)
	cut -f2 -d " " $ped > ${ped}_ind.txt
			DATE=$(echo $ped | grep -o -E '_[0-9]+.ped')
			DATE=${DATE#_}
			DATE=${DATE%.ped}	
			awk -v column=2 -v value="${CH_IP}" '
			    BEGIN {
				FS = OFS = " ";
			    }
			    {
				for ( i = NF + 1; i > column; i-- ) {
				    $i = $(i-1);
				}
				$i = value;
				print $0;
			    }' ${ped}_ind.txt > ${ped}_IND.txt
			awk -v column=3 -v value="${DATE}" '
				    BEGIN {
					FS = OFS = " ";
				    }
				    {
					for ( i = NF + 1; i > column; i-- ) {
					    $i = $(i-1);
					}
					$i = value;
					print $0;
				    }' ${ped}_IND.txt > ${ped}_INDI.txt

	cat *INDI.txt > ${SERNUM}_INDIVIDUALS.txt
	mv ${SERNUM}_INDIVIDUALS.txt $TEMPDIR
	rm *ind*
	rm *INDI.txt
	rm *IND.*
	

	###########################################################################################
	#extract the 800 SNPs needed for parental verification
	#move them to temporary directory	
	#######################################################################################
	ped1=$(echo $ped |  sed "s/.ped$//")
	~/bin/plink --file $ped1 --cow --extract $PARSNP_PATH/Names_800SNPs.txt  -recode --out ${SERNUM}_chip800
	

	    	
	#########################################################################################
	#ZANARDI to update the PLINK_MERGED files for the chip with addition of new genotypes
	# once you have the names of the new ped and map files cd to the chip directory in the Genotipi_latest 
	##########################################################################################
	cd $GENLAT_PATH/$CH_IP 
	cp $CODE_DIR/PARAMFILE.txt . #copy PARAMFILE for Zanardi to the directory and edit it for adding new files to the MERGED FILE
	sed -i "s%PathToPed%$pedfiles,$PWD/OUTPUT/PLINK_MERGED.ped%g" PARAMFILE.txt
	sed -i "s%PathToMap%$mapfiles,$PWD/OUTPUT/PLINK_MERGED.map%g" PARAMFILE.txt
	sed -i "s%OutputName%${CH_IP}_Merged%g" PARAMFILE.txt
	cp PARAMFILE.txt $ZAN_PATH	
	python $ZAN_PATH/Zanardi.py --mds
	

	cd $TEMPDIR	

done

#you are now in the temporary directory

##################################################################################3
#combine the files with newly genotyped individuals from all the chips once all the loops are done
#################################################################################
cat *_INDIVIDUALS.txt > GenotypedInd
rm *INDIVIDUALS.txt
sed -i 's/.ped//g' GenotypedInd
sed -i 's/_//' GenotypedInd
ALLANI=$(less GenotypedInd | wc -l)
echo "The number of newly genotyped animals on all chips is $ALLANI"

#run R script to tranform into Govedo format
#copy the script and edit the input files
#the output is writen to $TEMPIDR/Genotyped_individuals.csv
cp $CODE_DIR/NewIndForGovedo.R $TEMPDIR
sed -i "s%IDSeqFile%$IDSEQ_FILE%g" $TEMPDIR/NewIndForGovedo.R

#!/usr/bin/enc Rscript
Rscript $TEMPDIR/NewIndForGovedo.R

RANI=$(less Genotyped_individuals.csv | wc -l) #check whether all genotyped animals are found a sequence
RANI=$(($RANI - 1))
echo "The number of animals with found corresponding sequence is $RANI"


########################################################################################
#tranform in a format suitable for Govedo PARENTAL_SNP800
##############################################################################

#the output PLINK_MERGED ped and map are in the $TEMPDIR/OUTPUT
#cd to the OUTPUT directory


cp $CODE_DIR/Add_alleles_Govedo_PARENTAL_SNP800.R .

#11) ped and map SNP800 files of the sample
SNP800GENO_FILE="'$TEMPDIR/$(ls *chip800.ped)'"
SNP800MAP_FILE="'$TEMPDIR/$(ls *chip800.map)'"

#run R script to tranform into Govedo format
#requires files: sifrant SNPov, file with animal IDs and sequences, 
sed -i "s%SifrantFile%$SIF_FILE%g" Add_alleles_Govedo_PARENTAL_SNP800.R
sed -i "s%IDSeqFile%$IDSEQ_FILE%g" Add_alleles_Govedo_PARENTAL_SNP800.R
sed -i "s%SNP800genoFile%$SNP800GENO_FILE%g" Add_alleles_Govedo_PARENTAL_SNP800.R
sed -i "s%SNP800mapFile%$SNP800MAP_FILE%g" Add_alleles_Govedo_PARENTAL_SNP800.R
sed -i "s%Govedo800SNPFile%'$TEMPDIR/GovedoSNP800.csv'%g" Add_alleles_Govedo_PARENTAL_SNP800.R


#!/usr/bin/enc Rscript
Rscript $TEMPDIR/Add_alleles_Govedo_PARENTAL_SNP800.R



	





