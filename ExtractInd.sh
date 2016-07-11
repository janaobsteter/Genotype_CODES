for each in 50Kv02 GGPv02 GGPv03  GGPv04 HD HDv02
do
cd $each
	for i in Matija*.ped
	do
		cut -f2 -d " " $i > ${i}_ind.txt
		DATE=$(echo $i | grep -o -E '_[0-9]+.ped')
		awk -v column=2 -v value="${each}" ' #this is chip name added to the second column
		    BEGIN {
		        FS = OFS = " ";
		    }
		    {
		        for ( i = NF + 1; i > column; i-- ) {
		            $i = $(i-1);
		        }
		        $i = value;
		        print $0;
		    }' ${i}_ind.txt > ${i}_IND.txt
		awk -v column=3 -v value="${DATE}" '
			    BEGIN {
				FS = OFS = " ";
			    }
			    {
				for ( i = NF + 1; i > column; i-- ) {
				    $i = $(i-1);
				}
				$i = value;
				print $0;
			    }' ${i}_IND.txt > ${i}_INDI.txt
	done
cat *INDI.txt > ${each}_INDIVIDUALS.txt
mv ${each}_INDIVIDUALS.txt ..
rm *ind*
rm *INDI.txt
rm *IND.*
cd ..
done
cat *_INDIVIDUALS.txt > GenotypedInd
rm *INDIVIDUALS.txt
sed -i 's/.ped//g' GenotypedInd
sed -i 's/_//' GenotypedInd
sed -i -e 's/^/SI/' GenotypedInd
sed -i 's/SISI/SI/g' GenotypedInd
sed -i 's/SI4384195/SI04384195/g' GenotypedInd
sed -i 's/SI4574059/SI04574059/g' GenotypedInd
