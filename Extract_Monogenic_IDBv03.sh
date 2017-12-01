#first extract the Monogenic Diseases / without milk proteins
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_MonogenicMap.txt --recode --out $1_Monogenic


#then extract the Beta-Casein SNPs
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_BetaA2.txt --recode --out $1_BetaA2
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_BetaB.txt --recode --out $1_BetaB

#extract the kappa-casein SNPs
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_KapaAB_1.txt --recode --out $1_KapaAB_1
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_KapaAB_2.txt --recode --out $1_KapaAB_2
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_KapaC.txt --recode --out $1_KapaC
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_KapaE.txt --recode --out $1_KapaE
plink --file $1 --cow --extract ~/Genotipi/Genotipi_CODES/IDBv03_KapaI.txt --recode --out $1_KapaI

