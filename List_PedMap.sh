 #!/bin/bash

find Matija_Rigler_BOV50*/*.ped > ped50Kfiles.txt
find Matija_Rigler_BOVGGP*/*.ped > pedGGPv03files.txt
find Matija_Rigler_BOVGP4*/*.ped > pedGP4v01files.txt
find Matija_Rigler_BOVGP3*/*.ped > pedGP3v02files.txt
find Matija_Rigler_Bovine*/*.ped > pedHDfiles.txt

find Matija_Rigler_BOV50*/*.map > map50Kfiles.txt
find Matija_Rigler_BOVGGP*/*.map > mapGGPv03files.txt
find Matija_Rigler_BOVGP4*/*.map > mapGP4v01files.txt
find Matija_Rigler_BOVGP3*/*.map > mapGP3v02files.txt
find Matija_Rigler_Bovine*/*.map > mapHDfiles.txt

cat ped50Kfiles.txt map50Kfiles.txt | sort -n -r > MapPed50K.txt
cat pedGGPv03files.txt mapGGPv03files.txt | sort -n -r > MapPedGGPv03.txt
cat pedGP4v01files.txt mapGP4v01files.txt | sort -n -r > MapPedGP4v01.txt
cat pedGP3v02files.txt mapGP3v02files.txt | sort -n -r > MapPedGP3v02.txt
cat pedHDfiles.txt mapHDfiles.txt | sort -n -r > MapPedHD.txt
