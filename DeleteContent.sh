for pg in "1_1" "1_2" "2_1"
do
	for ref in "True" "False"
	do
		for ph in 1 2 5 8 9 10
		do
			> Pheno_Gen0_${ref}${ph}_${pg}.txt
		done
	done
done
