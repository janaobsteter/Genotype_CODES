homeDir=$PWD

for rep in 4
do
	cp ${homeDir}/CodeDir/Relatedness* FillInBurnIn${rep}/GA/
	cd ${homeDir}/FillInBurnIn${rep}/GA/
	sed -i "s/_REPEATFILLIN_/${rep}/g" Relatedness_Loop.sh
	bash Relatedness_Loop.sh
done
	
