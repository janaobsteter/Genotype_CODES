
#IF EFFECT IS 1!
awk '{if ($1==1 && $2==1) print $3,$4}' gibbsSol2 | sort +0 -1 > sol.temp
awk '{print $1,$10}' renadd01.ped | sort +0 -1 > ids.temp
join -1 +1 -2 +1 ids.temp sol.temp | sort -n -k2,2 >   renumbered_Solutions
