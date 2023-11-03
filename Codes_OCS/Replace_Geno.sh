
missing=$(perl -E "print '5' x 50")

for i in $(seq 1 950)
do
	missing=$(perl -E "print '5' x 50")
	linenum=${i}s
	sed -Ei "${linenum}/[0,1,2]{50}/$missing/" GenoFile.txt
done

