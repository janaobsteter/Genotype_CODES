for rep in $(seq 0 19)
do
	cd FillInBurnIn$rep
	rm SelectionParam_*
	cd ..
done
