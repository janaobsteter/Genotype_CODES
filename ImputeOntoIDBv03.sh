#!/bin/bash


os.system('plink --file ' + WorkingDir +'/OUTPUT/FIMPUTE_Imputed' + str(masking) + ' --cow --keep ImputedIndsConc.txt --recode --out ImputedMasking' + str(masking))
os.system('plink --file ' + WorkingDir +'/MERGEDforImputation --cow --keep ImputedIndsConc.txt --recode --out OrigMasking' + str(masking))


#Now compute the allelic concordance
chkConc = Concordance('OrigMasking' + str(masking) + '.ped', 'ImputedMasking' + str(masking) + '.ped')
os.system('cut -f1 ClusterSNPs_' + str(masking) + '.txt -d " " > ClusterList' + str(masking) + '.txt' )
chkConc.concSNPs('ClusterList' + str(masking) + '.txt', 'Concordance' + str(masking))
chkConc.extractConc('Concordance' + str(masking) + '.txt')
