#!/bin/bash

~/bin/plink --file FIMPUTE_Imputed --cow --mind 0.01 --maf 0.01 --geno 0.1 --hwe 0.0001 --recode --out FIMPUTE_ImputedQC

