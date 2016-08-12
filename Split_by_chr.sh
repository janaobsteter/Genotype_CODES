 for chr in $(seq 1 32) ; do /home/janao/bin/plink --file FIMPUTE_Imputed --cow --chr $chr --recode --out Imputed_chr$chr; done 
