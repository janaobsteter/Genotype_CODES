#!/bin/bash

#############################################################
#combine ped files by chip using Zanardi
################################################################


cd /home/janao/Genotipi/Genotipi_DATA/Genotipi_latest/Rjava

for each in 50Kv01 50Kv02 GGPv02 GGPv03 GGPv04 HD HDv02
do
	cd $each
	cp /home/janao/Genotipi/Genotipi_CODES/PARAMFILE.txt .
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
	python /home/janao/Genotipi/Genotipi_CODES/Zanardi/Zanardi.py --mds
	cd ..
done
	

