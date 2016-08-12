#script za ureditev peddar.param file s parametri za tek PEDDA_ROW python programa
#PEDDA_ROW = program za raw Illumina --> ped in map
#za vse subdirektorije v pwd
PATH_PR=/home/janao/Genotipi/SNPchimpRepo/source_codes/PEDDA_ROW
PASMA="Rjava"
for i in $(find $PWD -maxdepth 1 -type d); 
do
cp $PATH_PR/pedda_row.py $i #copy python script to each subdirectory
cp $PATH_PR/peddar.param $i #coppy parameter file in each subdirectory
cd $i
_FINREP=$(find *FinalReport.txt -printf '%f\n') #find FinalReport
SERNUM=${_FINREP%_FinalReport.txt} #remove sufix to only get the number of the genotyping package
sed -i 's/SI /SI/g' $_FINREP #remove blank space from SI #####
sed -i "s/test_FinalReport.txt/"${_FINREP}"/g" peddar.param #change finrep parameter in place
sed -i "s/test_outputfile/"${SERNUM}"/g" peddar.param #change output prefix parameter in place
sed -i "s/TEST/$PASMA/g" peddar.param #change brdcode parameter in place
python pedda_row.py
done
