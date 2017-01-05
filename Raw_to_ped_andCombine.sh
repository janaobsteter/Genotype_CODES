#script za ureditev peddar.param file s parametri za tek PEDDA_ROW python programa
#PEDDA_ROW = program za raw Illumina --> ped in map
#za vse subdirektorije v pwd

GENOTIPDIR=/home/janao/Genotipi/Genotipi_DATA/Rjava_TEMP/
PATH_PR=/home/janao/Genotipi/Genotipi_CODES/SNPchimpRepo/source_codes/PEDDA_ROW
PASMA="Rjava"

mkdir $GENOTIPDIR/GCpeds


for i in $(find $GENOTIPDIR -print | egrep *FinalReport.txt$  |cat)
do
	SERNUM=${$i%_FinalReport.txt}_gc #remove sufix to only get the number of the genotyping package
	sed -i 's/SI /SI/g' $i #remove blank space from SI #####
	sed -i "s/test_FinalReport.txt/"${i}"/g" peddar.param #change finrep parameter in place
	sed -i "s/test_outputfile/"${SERNUM}"/g" peddar.param #change output prefix parameter in place
	sed -i "s/TEST/$PASMA/g" peddar.param #change brdcode parameter in place
	python pedda_row.py
	cp $SERNUM.ped $GENOTIPDIR/GCpeds
	cp $SERNUM.map $GENOTIPDIR/GCpeds
done


#	mkdir GCpeds
#	cd GCpeds
#	cp /home/janao/Genotipi/Genotipi_CODES/PARAMFILE.txt .
#	ls -d ../*_a.ped > ${i}_ped.txt
#	ls -d ../*_a.map > ${i}_map.txt
#	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${i}_ped.txt
#	sed -n -i '1h;2,$H;${g;s/\n/,/g;p}' ${i}_map.txt
#	pedfiles=$(cat ${i}_ped.txt)
#	mapfiles=$(cat ${i}_map.txt)
#	sed -i "s%PathToPed%$pedfiles%g" PARAMFILE.txt
#	sed -i "s%PathToMap%$mapfiles%g" PARAMFILE.txt
#	sed -i "s%OutputName%${i}_A_Merged%g" PARAMFILE.txt
#	cp PARAMFILE.txt /home/janao/Genotipi/Zanardi	
#	python /home/janao/Genotipi/Genotipi_CODES/Zanardi/Zanardi.py --mds
#	cd ..
