for rep in $(seq 0 20)
do
	echo $rep
	less 10K_SireUse/SU55/Gen${rep}/IndForGeno_new.txt | wc -l
done
