 for chr in $(seq 1 32) ; do /home/janao/bin/plink --file FIMPUTE_GP4_GSImputed --cow --chr $chr --recode --out Imputed_chr$chr; done 
