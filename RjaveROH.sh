#!/bin/bash
CHIP=50KImpQC

sudo ~/bin/plink --file FIMPUTE_ImputedQC --cow --homozyg --homozyg-density 120 --homozyg-gap 1000 --homozyg-window-missing 2 --homozyg-window-het 1 --homozyg-window-snp 50 --out Rjave${CHIP}_ROH2
sudo ~/bin/plink --file FIMPUTE_ImputedQC --cow --homozyg --homozyg-density 1000 --homozyg-window-missing 2 --homozyg-window-het 1 --homozyg-window-snp 20 --homozyg-kb 10 --out Rjave${CHIP}_ROH1
#sudo ~/bin/plink-1.07-x86_64/plink --file FIMPUTE_ImputedQC --cow --homozyg --homozyg-window-kb 5000 --homozyg-gap 250 --homozyg-density 50 --homozyg-window-het 1 --out WOORjave${CHIP}_ROH --noweb

#cut -c11-21 Rjave${CHIP}_ROH.hom.indiv   > Rjave${CHIP}ROHSeq.txt
#cut -c37-44 Rjave${CHIP}_ROH.hom.indiv > Rjave${CHIP}ROHKb.txt
#cut -c183-188 Rjave${CHIP}_ROH.hom > Rjave${CHIP}PHET.txt
#paste Rjave${CHIP}ROHSeq.txt Rjave${CHIP}ROHKb.txt Rjave${CHIP}PHET.txt> Rjave${CHIP}Data.txt -d ","
#sed -i 's/ ,/,/g' Rjave${CHIP}Data.txt
#sed -i 's/, /,/g' Rjave${CHIP}Data.txt



#cut -c11-21 WOORjave${CHIP}_ROH.hom.indiv   > WOORjave${CHIP}ROHSeq.txt
#cut -c37-44 WOORjave${CHIP}_ROH.hom.indiv > WOORjave${CHIP}ROHKb.txt
#paste WOORjave${CHIP}ROHSeq.txt WOORjave${CHIP}ROHKb.txt > WOORjave${CHIP}Data.txt -d ","
#sed -i 's/ ,/,/g' WOORjave${CHIP}Data.txt
#sed -i 's/, /,/g' WOORjave${CHIP}Data.txt
#sed -i 's/, /,/g' WOORjave${CHIP}Data.txt


