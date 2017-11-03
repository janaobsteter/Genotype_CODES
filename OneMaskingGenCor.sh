
for i in $(seq 0 9)
do 
	#plink --file MERGEDforImputation --cow --extract ClusterList$i.txt --recodeA --out ORIGA$i
	#plink --file ./OUTPUT/FIMPUTE_Imputed$i --cow --extract ClusterList$i.txt --recodeA --out IMPUTEDA$i
	cp /home/jana/Genotipi/Genotipi_CODES/GenotipicCorrelation.R .
	sed -i "s/IMPUTEDRAWFILE/IMPUTEDA${i}/g" GenotipicCorrelation.R
	sed -i "s/ORIGINALRAWFILE/ORIGA${i}/g" GenotipicCorrelation.R
	sed -i "s/NUMBEROFMASKING/${i}/g" GenotipicCorrelation.R
	Rscript GenotipicCorrelation.R
done
	
