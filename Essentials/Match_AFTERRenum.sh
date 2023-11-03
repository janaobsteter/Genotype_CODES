
#IF EFFECT IS 1!
awk '{if ($1==1 && $2==2) print $3,$4}' solutions | sort +0 -1 > sol.temp
awk '{print $1,$10}' renadd02.ped | sort +0 -1 > ids.temp
join -1 +1 -2 +1 ids.temp sol.temp | sort -n -k2,2 >   renumbered_Solutions
#the first field is renumf90 ID, the second is origin ID
