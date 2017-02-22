##  make my data
vseRjavePDF <- read.csv("~/Documents/F4F/OdbiraZivali/RjaveKrave_937_22022017_PDF.csv",header=T) 
seznamB <- read.csv("~/Documents/F4F/OdbiraZivali/seznamB_22022017.csv", header=T)

## knitr loop
library("knitr")
for (creda in 32717){
  knit2pdf("testingloops.Rnw", output=paste0('report_', creda, '.tex'))
}

## knitr loop

knit2pdf("test.Rnw", output='test.tex')




latex(credaTable,table=F,center='centering',file='',
      booktabs=T,numeric.dollar=F,colheads=c("Col A","Col B","Col C"),
      colnamesTexCmd="bfseries", col.just=c("S","S","S", "S"))



library(xtable)
#Just some random data
x1 <- rnorm(1000); x2 <- rnorm(1000); x3 <- rnorm(1000)
y  <- 2 + 1 *x1 + rnorm(1000)
#Run regressions
reg1 <- summary(lm(y ~ x1))
reg2 <- summary(lm(y ~ x2))
reg3 <- summary(lm(y ~ x3))
#Create data.frame
df <- data.frame(Model = 1:3,
                 Alpha = c(reg1$coef[1,1], reg2$coef[1,1], reg3$coef[1,1]),
                 Beta  = c(reg1$coef[2,1], reg2$coef[2,1], reg3$coef[2,1]),
                 tV    = c(reg1$coef[2,2], reg2$coef[2,2], reg3$coef[2,2]),
                 AdjR  = c(reg1$adj.r.s,  reg2$adj.r.s,   reg3$adj.r.s))
strCaption <- paste0("\\textbf{Table Whatever} This table is just produced with some",
                     "random data and does not mean anything. Just to show you how ",
                     "things work.")
print(xtable(df, digits=2, caption=strCaption, label="Test_table"), 
      size="footnotesize", #Change size; useful for bigger tables
      include.rownames=FALSE, #Don't print rownames
      include.colnames=FALSE, #We create them ourselves
      caption.placement="top", 
      hline.after=NULL, #We don't need hline; we use booktabs
      add.to.row = list(pos = list(-1, 
                                   nrow(df)),
                        command = c(paste("\\toprule \n",
                                          "Model & $\\alpha$ & $\\beta$ & t-value & $R^2$ \\\\\n", 
                                          "\\midrule \n"),
                                    "\\bottomrule \n")
      )
)