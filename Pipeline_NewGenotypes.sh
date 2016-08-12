DATE=07072016 #datum genotipizacije (v imenu filea)
FILE=Matija_Rigler_BOV50KV02_20160306.zip # downloaded file


mkdir ~/Genotipi/Genotipi_${DATE} #create temp directory for new genotypes
TEMPDIR=~/Genotipi/Genotipi_${DATE}


cp ~/Downloads/$FILE $TEMPDIR #copy downloaded genotyped to the temp direcotry
cd $TEMPDIR
unzip * #unzip - may containg multiple folders - packages

#create associative array to hold the number of SNPs on each chip
declare -A chip
chip[19720]=GGPv02
chip[26151]=GGPv03
chip[30105]=GGPv04
chip[76883]=HD
chip[138892]=HDv02
chip[54001]=50Kv01
chip[54609]=50Kv02


#create and empty array to hold the chips of the new files --> merged plinks that need to be modified 
MODCHIP=()

for i in */ #for all subdirectories - i.e. all the genotype "packages" downloaded
do
	cd $i
	#unzip SNP_Map and FinalReport which are needed for raw to ped tranformation
	unzip SNP_Map.zip
	unzip *FinalReport.zip
	#Raw into PEd and Map with PEDDAROW
	#script for editing parameter file for PEDDAROW
	#PEDDA_ROW = python program for raw Illumina --> ped in map
	PATH_PR=/home/janao/Genotipi/SNPchimpRepo/source_codes/PEDDA_ROW
	PASMA="Rjava"
	cp $PATH_PR/pedda_row.py . #copy python script to each subdirectory
	cp $PATH_PR/peddar.param . #coppy parameter file in each subdirectory
	#cd $i
	_FINREP=$(find *FinalReport.txt -printf '%f\n') #find FinalReport
	SERNUM=${_FINREP%_FinalReport.txt} #remove sufix to only get the number of the genotyping package
	sed -i 's/SI /SI/g' $_FINREP #remove blank space from SI #####
	sed -i "s/test_FinalReport.txt/"${_FINREP}"/g" peddar.param #change finrep parameter in place
	sed -i "s/test_outputfile/"${SERNUM}"/g" peddar.param #change output prefix parameter in place
	sed -i "s/TEST/$PASMA/g" peddar.param #change brdcode parameter in place
	python pedda_row.py

	#check whether there are any animal names in the ped files --> if there are, terminate the program and notify
	read LINES < <(less *.ped | cut -f3 -d " " | grep [A-Z] | wc -l)
	if [ $LINES -gt 0 ]
		then
		echo " "
		echo "WARNING!!!"
		echo "Ped file in ${i} contains animal names in the third column"
		NAMES=$(less *.ped | cut -f3 -d " " | grep [A-Z])
		for i in $NAMES
		do 
			sed -i "s/${i} //g" *.ped
		done
	fi
	
	echo "No animals names found in the ped file, proceding with the script"
	#rm ~/Downloads/$FILE #everything is OK, you can remove the ziped downloaded file

	#read the number of SNPs on the chip and move to the right chip directory within Genotipi/Genotipi_latest
	read NUMSNP < <( less *map | wc -l)
	CH_IP=$(echo ${chip[$NUMSNP]})
	echo "The chip is $CH_IP"
	#the number of genotyped animals in the package
	read ANINUM < <( less *ped | wc -l)
	echo "Number of genotypes animals is $ANINUM"

	cp *ped ~/Genotipi/Genotipi_latest/$CH_IP
	cp *map ~/Genotipi/Genotipi_latest/$CH_IP
	echo "Copying ped and map files from $i to Genotipi_latest/$CH_IP"


	MODCHIP+=,$CH_IP # in each round add chip name to the array

	cd ..
done

cd ~/Genotipi/Genotipi_latest
echo $MODCHIP > modchip.txt #write to a file in order to get unique values
sed -i 's/,/\n/g' modchip.txt #replace commas with new line
less modchip.txt | sort | uniq > modchip_uniq.txt #sort and get unique values
file=modchip_uniq.txt 
rm modchip.txt
read -d $'\x04' chips < "$file" #chips variable now holds unique chips than need to be updated

echo "Modifying PLINK_MERGED files for $CH_IP"
for each in $chips
do
	cd $each
	cp /home/janao/Documents/PARAMFILE.txt .
	ls -d $PWD/*.ped > ${each}_ped.txt
	ls -d $PWD/*.map > ${each}_map.txt
	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${each}_ped.txt
	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${each}_map.txt
	pedfiles=$(cat ${each}_ped.txt)
	mapfiles=$(cat ${each}_map.txt)
	sed -i "s%PathToPed%$pedfiles%g" PARAMFILE.txt
	sed -i "s%PathToMap%$mapfiles%g" PARAMFILE.txt
	sed -i "s%OutputName%${each}_Merged%g" PARAMFILE.txt
	cp PARAMFILE.txt /home/janao/Genotipi/Zanardi	
	python /home/janao/Genotipi/Zanardi/Zanardi.py --mds
	cd ..
done
	






