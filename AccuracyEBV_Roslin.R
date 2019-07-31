RoslinBlue    <- rgb(red= 55, green=152, blue=217, maxColorValue=255) ## blue
RoslinBlueL   <- rgb(red=157, green=196, blue=224, maxColorValue=255) ## blue light
RoslinGreen   <- rgb(red= 95, green=168, blue= 59, maxColorValue=255) ## green
RoslinGreenL  <- rgb(red=178, green=227, blue=156, maxColorValue=255) ## green light
RoslinOrange  <- rgb(red=219, green= 79, blue= 15, maxColorValue=255) ## orange
RoslinOrangeL <- rgb(red=223, green=177, blue=156, maxColorValue=255) ## orange light
RoslinViolet  <- rgb(red=193, green= 23, blue=115, maxColorValue=255) ## violet
RoslinVioletL <- rgb(red=222, green=153, blue=191, maxColorValue=255) ## violet light
RoslinGray    <- rgb(red= 83, green= 83, blue= 83, maxColorValue=255) ## gray
RoslinGrayL   <- rgb(red=208, green=205, blue=207, maxColorValue=255) ## gray light
RoslinCol <- c(RoslinBlue, RoslinGreen, RoslinOrange, RoslinViolet, RoslinGray)

n <- 1:20
h2 <- 0.25
p2 <- c(0,0.05,0.10,0.15,0.20)
r2 <- h2+p2
Acc <- matrix(nrow=length(n), ncol=length(p2))
for (i in seq_along(p2)) {
  Acc[,i] <- sqrt(n*h2/(1+(n-1)*r2[i]))
}
matplot(y=Acc, x=n, type="l", col=RoslinCol, lty=1, lwd=2, ylim=c(0.4, 0.98))
legend(x="bottomright", legend=r2, col=RoslinCol, lty=1, lwd=2, bty="n", title="Repeatability")
## Loss in Acc by doing 9 instead of 10 test-days
Acc[10, ] - Acc[9,]

library(Rmisc)



acc <- read.csv("Accuracy_permEnv_selCat.csv")
acc[order(acc$AgeCat),]
