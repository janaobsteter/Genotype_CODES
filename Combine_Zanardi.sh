#!/bin/bash

#############################################################
#combine ped files by chip using Zanardi
################################################################


cd /home/janao/Genotipi/Genotipi1_12042016

for each in 50K GGP GGPv03 GP3v02 GP4 HD
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
	

