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
GEN_PATH=~/Genotipi/Genotipi_DATA
#4) path to the directory where you store your latest and final genotypes
GENLAT_PATH=~/Genotipi/Genotipi_DATA/Genotipi_latest 
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
DATEDOWNLOAD=11082016 #date of downloading the file
FILE=Matija_Rigler.zip # downloaded file
#|||||||||||||||||||||||||||||||||||||||||||||||||



TEMPDIR=$GEN_PATH/Genotipi_${DATEDOWNLOAD}




cd $TEMPDIR

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

	#rm $DOWNLOAD_PATH/$FILE #everything is OK, you can remove the ziped downloaded file

	############################################################################
	#read the number of SNPs on the chip and copy to the right chip directory within Genotipi/Genotipi_latest
	############################################################################
	read NUMSNP < <( less *map | wc -l)
	CH_IP=$(echo ${chip[$NUMSNP]})
	echo "The chip is $CH_IP"

	##########################################################################
	#get the list of newly obtained ped and map files for Zanardi PARAMFILE edit
	#obtain them in this directory since in the Genotipi_latest there are all the ped and map file
	##########################################################################
	ls -d $PWD/*.ped > ${CH_IP}_ped.txt #list all new ped and map files - for Zanardi update of the PLINK_MERGED file - 		just to add the new files
	ls -d $PWD/*.map > ${CH_IP}_map.txt
	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${CH_IP}_ped.txt replace new lines with commas
	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${CH_IP}_map.txt
	pedfiles=$(cat ${CH_IP}_ped.txt)
	mapfiles=$(cat ${CH_IP}_map.txt)

	
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


