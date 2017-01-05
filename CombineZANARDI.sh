for i in .
do
	
	cp /home/janao/Genotipi/Genotipi_CODES/PARAMFILE.txt .
	ls -d $PWD/*.ped > Ped.txt
	ls -d $PWD/*.map > Map.txt

#	pedfiles=$(cat Ped.txt)
#	mapfiles=$(cat Map.txt)
#	sed -i "s%PathToPed%$pedfiles%g" PARAMFILE.txt
#	sed -i "s%PathToMap%$mapfiles%g" PARAMFILE.txt
#	sed -i "s%OutputName%$STRAND_Merged%g" PARAMFILE.txt
#	cp PARAMFILE.txt /home/janao/Genotipi/Zanardi	
#	python /home/janao/Genotipi/Genotipi_CODES/Zanardi/Zanardi.py --mds
#	cd ..
done
