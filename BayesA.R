#SAME FOR INTERCEPT AND alpha EFFECT
# SAME FOR VARIANCE E AND VARIANCE alpha

# Parameters
# Parameters
nmarkers = 2000;   #number of markers
numiter  = 2000;   #number of iterations
vara     = 1.0/20.0; 

cat("Rohan Fernando's implementation of Bayes A\n")
# input data
data     = matrix(scan("trainData.out0"),ncol=nmarkers+2,byrow=TRUE);
nrecords = dim(data)[1];
startMarker = 1901;
x = cbind(1,data[,startMarker:nmarkers]);  #this is the mean and then the markers
y = data[,nmarkers+1];
a =  data[,nmarkers+2];
# inital values

nmarkers = nmarkers - startMarker + 1;
mean2pq = 0.5;                                         #just an approximation
scalea  = 0.5*vara/(nmarkers*mean2pq);                 # slide 51  0.5 = (v-2)/v

size = dim(x)[2];
b = array(0.0,size);
meanb = b;
b[1] = mean(y);
var  = array(0.0,size);

# adjust y
 ycorr = y - x%*%b;                  

# mcmc sampling
for (iter in 1:numiter){
	
# sample vare
	vare = ( t(ycorr)%*%ycorr )/rchisq(1,nrecords + 3); # slide 30
	
# sample intercept
	ycorr = ycorr + x[,1]*b[1];
	rhs    = sum(ycorr)/vare;
	invLhs = 1.0/(nrecords/vare);
	mean = rhs*invLhs;                            
	b[1] = rnorm(1,mean,sqrt(invLhs));                  # slide 26
	ycorr = ycorr - x[,1]*b[1];
	meanb[1] = meanb[1] + b[1];
	
# sample variance for each locus
	
	for (locus in 2:size){
		var[locus] = (scalea*4+b[locus]*b[locus])/rchisq(1,4.0+1) # slides 35-36
	}
# sample effect for each locus	
	for (locus in 2:size){
		ycorr = ycorr + x[,locus]*b[locus];   #unadjust y for this locus
		rhs = t(x[,locus])%*%ycorr/vare;
		lhs = t(x[,locus])%*%x[,locus]/vare + 1.0/var[locus];
		invLhs = 1.0/lhs;
		mean = invLhs*rhs;
		b[locus]= rnorm(1,mean,sqrt(invLhs));                     # slide 28
		ycorr = ycorr - x[,locus]*b[locus];   #adjust y for the new value of this locus
		meanb[locus] = meanb[locus] + b[locus];		
	}
	if ((iter %% 100)==0) cat ("iteration ",iter,"\n");
	
}

meanb = meanb/numiter;
aHat  = x %*% meanb;
corr = cor(a,aHat);
cat ("corr = ",corr, "\n");
plot(a,aHat)
