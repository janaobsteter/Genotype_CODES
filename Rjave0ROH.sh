#!/bin/bash
CHIP=50KImpQC

sudo ~/bin/plink --file FIMPUTE_ImputedQC --cow --homozyg --out Rjave0${CHIP}_ROH
#sudo ~/bin/plink-1.07-x86_64/plink --file FIMPUTE_ImputedQC --cow --homozyg --homozyg-window-kb 5000 --homozyg-gap 250 --homozyg-density 50 --homozyg-window-het 1 --out WOORjave${CHIP}_ROH --noweb

cut -c11-21 Rjave0${CHIP}_ROH.hom.indiv   > Rjave0${CHIP}ROHSeq.txt
cut -c37-44 Rjave0${CHIP}_ROH.hom.indiv > Rjave0${CHIP}ROHKb.txt
cut -c183-188 Rjave0${CHIP}_ROH.hom > Rjave0${CHIP}PHET.txt
paste Rjave0${CHIP}ROHSeq.txt Rjave0${CHIP}ROHKb.txt Rjave0${CHIP}PHET.txt> Rjave0${CHIP}Data.txt -d ","
sed -i 's/ ,/,/g' Rjave0${CHIP}Data.txt
sed -i 's/, /,/g' Rjave0${CHIP}Data.txt



#cut -c11-21 WOORjave${CHIP}_ROH.hom.indiv   > WOORjave${CHIP}ROHSeq.txt
#cut -c37-44 WOORjave${CHIP}_ROH.hom.indiv > WOORjave${CHIP}ROHKb.txt
#paste WOORjave${CHIP}ROHSeq.txt WOORjave${CHIP}ROHKb.txt > WOORjave${CHIP}Data.txt -d ","
#sed -i 's/ ,/,/g' WOORjave${CHIP}Data.txt
#sed -i 's/, /,/g' WOORjave${CHIP}Data.txt
#sed -i 's/, /,/g' WOORjave${CHIP}Data.txt


