trait=$1

./renumf90 < renumParam${trait}
./gibbs1f90 renf90.par --rounds 10000 --burnin 1000 --thin 10 --thinprint 1000
