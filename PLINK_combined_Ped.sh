#!/bin/bash

#########################################################################################################################
#combine ped and map files by chip
#some ped files problem due to animal names in the ped files --> manually removes
#########################################################################################################################

##run from Genotipi/Genotipi1XX directory

#combine 50K bed and map files onto 20160306 file
/home/janao/bin/plink-1.07-x86_64/plink --file /home/janao/Genotipi/Genotipi1_12042016/Matija_Rigler_BOV50KV02_20160306/Matija_Rigler_BOV50KV02_20160306 --merge-list MapPed50K.txt --recode --out Ped50K

#combine GP4v01 bed and map files onto 20160301 file
/home/janao/bin/plink-1.07-x86_64/plink --file /home/janao/Genotipi/Genotipi1_12042016/Matija_Rigler_BOVGP4V01_20160301/Matija_Rigler_BOVGP4V01_20160301 --merge-list MapPedGP4v01.txt --recode --out PedGP4v01

#combine GP3v02 bed and map files onto 20151013 file
/home/janao/bin/plink-1.07-x86_64/plink --file /home/janao/Genotipi/Genotipi1_12042016/Matija_Rigler_BOVGP3V02_20151013/Matija_Rigler_BOVGP3V02_20151013 --merge-list MapPedGP3v02.txt --recode --out PedGP3v02

#combine GGPv03 bed and map files onto 20150719 file
/home/janao/bin/plink-1.07-x86_64/plink --file /home/janao/Genotipi/Genotipi1_12042016/Matija_Rigler_BOVGGPV03_20150719/Matija_Rigler_BOVGGPV03_20150719 --merge-list MapPedGGPv03.txt --recode --out PedGGPv03

#combine HD bed and map files onto GGP_HD_29apr2014 file / 22 april
/home/janao/bin/plink-1.07-x86_64/plink --file /home/janao/Genotipi/Genotipi1_12042016/Matija_Rigler_Bovine_GGP_HD_29apr2014/Matija_Rigler_GGP_HD_22apr2014 --merge-list MapPedHD.txt --recode --out PedHD


