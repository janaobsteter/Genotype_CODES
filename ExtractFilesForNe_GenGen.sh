#mkdir MissingFiles
scenario="GenGen"
for rep in $(seq 0 9)
do
    for import in 0_0 10_10 25_25 50_50 100_100 0_100
    do
        for trait in 2 3
        do
	    if [ ! -d "MissingFiles/${scenario}${rep}_${import}${trait}" ] 
	    then
	            mkdir MissingFiles/${scenario}${rep}_${import}${trait}
	    fi

            tar -zxvf ${scenario}${rep}${import}${trait}.tar.gz home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/MeanHet_Marker_import.csv
            mv home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/MeanHet_Marker_import.csv MissingFiles/${scenario}${rep}_${import}${trait}/MeanHet_Marker_import.csv
            tar -zxvf ${scenario}${rep}${import}${trait}.tar.gz home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/MeanHet_Neutral_import.csv
            mv home/v1jobste/jobsteter/10K/SU55_import/${scenario}${rep}_${import}1${trait}/MeanHet_Marker_import.csv MissingFiles/${scenario}${rep}_${import}${trait}/MeanHet_Neutral_import.csv
        done
    done
done
	
