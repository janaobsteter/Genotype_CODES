attach(geyser)
plot(duration, waiting, xlim = c(0.5,6), ylim = c(40,100))
f1 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100))
lines(f1)
image(f1, zlim = c(0, 0.05))
f2 <- kde2d(duration, waiting, n = 50, lims = c(0.5, 6, 40, 100),
            h = c(width.SJ(duration), width.SJ(waiting)) )
image(f2, zlim = c(0, 0.05))
persp(f2, phi = 30, theta = 20, d = 5)
plot(duration[-272], duration[-1], xlim = c(0.5, 6),
     ylim = c(1, 6),xlab = "previous duration", ylab = "duration")
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(1.5, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels = c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(0.6, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels = c(0.05, 0.1, 0.2, 0.4) )
f1 <- kde2d(duration[-272], duration[-1],
            h = rep(0.4, 2), n = 50, lims = c(0.5, 6, 0.5, 6))
contour(f1, xlab = "previous duration",
        ylab = "duration", levels = c(0.05, 0.1, 0.2, 0.4) )
detach("geyser")

pedTest <- 
a <- 

plot(TGVs$zMean, TGVs$SDSt)
f1 <- kde2d(TGVs$zMean, TGVs$SDSt)
plot(f1)

library(devtools)
devtools::install_github("aloy/qqplotr")
install.packages("qqplotr")
library(ggplot2)

ped <- read.csv("~/TGVsAll_10KRef_2801_1Pb.csv", sep=" ")
ped <- ped[ped$Generation %in% 40:60,]
ped$group <- paste0(ped$scenario, ped$Rep)
ggplot(data=ped, aes(x=zMean, y=SDGenicSt, group=scenario, colour=scenario, fill=scenario)) + stat_density_2d(aes(fill = ..level..), geom = 'polygon') + 
  scale_fill_continuous() 

avg <- aggregate(ped$zMean)
p <- ggplot(data=ped, aes(y=zMean, x=SDGenicSt, group=group, colour=scenario)) + geom_point()
#p + geom_density2d(aes(y=zMean, x=SDGenicSt, colour=scenario, group=scenario),alpha=.5) + scale_fill_gradient(values=c("yellow", "red", "blue", "green", "brown")) 
p + stat_density_2d(aes(y=zMean, x=SDGenicSt, group=scenario, fill = ..level..), geom = 'polygon', alpha=0.05, bins=2) 

ped$zMean <- as.numeric(ped$zMean)
ped$SDGenicSt <- as.numeric(ped$SDGenicSt)
ggplot(data=ped, aes(y=zMean, x=SDGenicSt, group=scenario, colour=scenario)) + geom_smooth(method = "loess", formula = y ~ x) 
  
  stat_ellipse(inherit.aes = TRUE, level = 0.95)
  
  ##
  stat_summary(aes(y=ped$zMean ,x=ped$SDGenicSt, color=scenario, fill = ..level..), alpha = 0.1, 
                 geom = "polygon",data=ped, size=2, adjust = 1/8, contour=TRUE, args = list(shape1 = 2, shape2 = 2))

  
a <- hdr.2d(cars$speed, cars$dist, kde.package="ks", prob=0.05)
hdrscatterplot(cars$speed, cars$dist, kde.package="ks", levels = 1)  

library(ks)
set.seed(8192)
samp <- 200
mus <- rbind(c(-2,2), c(0,0), c(2,-2))
Sigmas <- rbind(diag(2), matrix(c(0.8, -0.72, -0.72, 0.8), nrow=2), diag(2))
cwt <- 3/11
props <- c((1-cwt)/2, cwt, (1-cwt)/2)
x <- rmvnorm.mixt(n=samp, mus=mus, Sigmas=Sigmas, props=props)
x <- ped[,c("SDGenicSt", "zMean")]
Hpi1 <- Hpi(x=x)
Hpi2 <- Hpi.diag(x=x)
fhat.pi1 <- kde(x=x, H=Hpi1)
fhat.pi2 <- kde(x=x, H=Hpi2)
plot(fhat.pi1)
plot(fhat.pi2)

plot(x[,1], x[,2])

library(car)
fit <- glm(ped$zMean ~ ped$SDGenicSt)
with(p, confidenceEllipse(fit, chisq=TRUE, levels = 0.95))



library(MASS) 

# parameters: 
n<-100 

# generate samples: 
set.seed(138813) 
#seed <- .Random.seed 
x<-rnorm(n); y<-rnorm(n) 

# estimate non-parameteric density surface via kernel smoothing 
den<-kde2d(x, y, n=n) 
# store z values of density estimate in an extra variable 
den.z <-den$z 

# this is the critical block, which I still do not comprehend in detail 
z <- array() 
for (i in 1:n){ 
  z.x <- max(which(den$x < x[i])) 
  z.y <- max(which(den$y < y[i])) 
  z[i] <- den$z[z.x, z.y] 
} 

# store class/level borders of confidence interval in variables 
confidence.border <- quantile(z, probs=1-0.6827, na.rm = TRUE) # +-1sd 

plot(x,y) 
par(new=TRUE) 
contour(den, levels=confidence.border, col = "red", add = TRUE) 


library(mvtnorm) # References rmvnorm()
library(ellipse) # References ellipse()
set.seed(17)

# Set the covariance matrix
sigma2 <- matrix(c(5, 2, 2, 5), ncol=2)

# Set the means
mu <- c(5,5)

# Get the correlation matrix
P <- cov2cor(sigma2)

# Generate the data
p <- rmvnorm(n=50, mean=mu, sigma=sqrt(sigma2))

# Plot the data
plot(p)

# Plot the ellipse
lines( ellipse( P, centre = c(5,5)) , col='red')

evals <- eigen(P)$values
evecs <- eigen(P)$vectors
# Angles of a circle
a <- seq(0, 2*pi, len=100)

# Get critical value
c2 <- qchisq(0.95, 2)
c <- sqrt(c2)

# Get the distances
xT <- c * sqrt(evals[1]) * cos(a)
yT <- c * sqrt(evals[2]) * sin(a)

M <- cbind(xT, yT)

transM <- evecs %*% t(M)
transM <- t(transM)

lines(transM + mu)
