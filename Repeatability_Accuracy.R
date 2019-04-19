#kg mleka
varG <- 8.2
varPE <- 3.9
varTE <- 9

#kg beljakovin
varG <- 0.008
varPE <- 0.004
varTE <- 0.011

varP <- varG + varPE + varTE


t <- (varG + varPE) / varP
t

#t <- 0.49
#t <- 0.35


r <- function (n, h2, t)  {
  sqrt((n * h2)  / (1 + (n-1)*t) )
}

acc <- data.frame(N = NA, h2 = NA, r = NA)

for (her in c(0.03, 0.1, 0.25, 0.35, 0.5)) {
  for (num in c(1:10)) {
    acc <- rbind(acc, c(num, her, r(num, her, t )))
  }
}

acc <- acc[-1,]
acc
acc$h2 <- as.factor(acc$h2)
acc$N <- as.factor(acc$N)
ggplot(data=acc, aes(x=N, y=r, group=h2, colour=h2)) + geom_line() + ggtitle(paste0("Repeatability: ", round(t, 2))) + xlab("Number of records") + 
  ylab("Accuracy") + scale_y_continuous(breaks=c(seq(0, 1.1, by=0.1), 1))
