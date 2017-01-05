for i in GGPv03 GGPv04 HD HDv02 50Kv01 50Kv02:
	do 
		cp ./$i/OUTPUT/PLINK_MERGED.ped .
		cp ./$i/OUTPUT/PLINK_MERGED.map .
		~/bin/plink --file OUTPUTs --cow --merge PLINK_MERGED --recode --out OUTPUTs
	done
