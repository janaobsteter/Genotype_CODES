for rep in $(seq 0 10)
do
	awk -v var=$rep '{print $0, var}' GenerationClass${rep}.txt > GenerationClass${rep}A.txt
done
