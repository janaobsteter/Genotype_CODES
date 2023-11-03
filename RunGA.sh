rep=$1

sed "s/_REPEAT_/${rep}/g" GA_qsub.sh > GA${rep}_qsub.sh
qsub GA${rep}_qsub.sh
