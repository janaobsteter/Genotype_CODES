for rep in $(seq 1 19)
do
	echo $(less Gen${rep}/IndForGeno_new.txt | wc -l)
done
