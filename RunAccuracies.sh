for strategy in "SU55" ###"SU51" "SU15"
do
	sed "s/_strategy_/$strategy/g" Accuracy_cor.sh  > Accuracy_cor_${strategy}.py
	qsub Accuracy_cor_${strategy}.py
done

