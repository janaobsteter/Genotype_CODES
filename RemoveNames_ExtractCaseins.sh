#####################################################################################################
#script for adding newly downlaoded genotypes on one or multiple chips or one and multiple genotype "packages"

####################################################################################


STRAND=Top

for i in .
do
	sed -i "s%SI %SI%g" $i/*ped #remove blank space from SI #####
	sed -i "s%SI %SI%g" $i/*ped #remove blank space from SI #####
	sed -i "s%SI %SI%g" $i/*ped #remove blank space from SI #####
	sed -i "s%SI %SI%g" $i/*ped #remove blank space from SI #####
	
	################################################################################
	#check whether there are any animal names in the ped files --> if there are, notify and remove them from ped file
	###################################################################################
	read LINES < <(less *.ped | cut -f3 -d " " | grep [A-Z] | wc -l)
	if [ $LINES -gt 0 ]
		then
		echo " "
		echo "WARNING!!!"
		echo "Ped file in ${i} contains animal names in the third column"
		NAMES=$(less *.ped | cut -f3 -d " " | grep [A-Z])
		for n in $NAMES
		do 
			sed -i "s/${n} //g" *.ped
		done
	fi
done


	
	




	





