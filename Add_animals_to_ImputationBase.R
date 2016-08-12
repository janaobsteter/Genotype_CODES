#add your animals to the imputation base
Chr1Base <- read.csv ("~/Genotipi/beagle/MS_imputation_base_chr1.bgl", sep=" ", header=F)
Chr1Bgl <- read.csv("~/Genotipi/Genotipi1_12042016/Imputation_HDref/OUTPUT/Chr1.bgl", sep=" ", header=F)
Chr1Bgl <- read.csv("/home/janao/Genotipi/Genotipi1_12042016/Imputation_GP4ref/OUTPUT/Chr1.bgl", sep=" ", header=F)
Chr1Bgl <- Chr1Bgl[-2,]

#capitalise in terminal
sed -i 's/.*/\U&/' Chr1.bgl

#69 lines match, only id line does not match - change in terminal
#merge them together
colnames(Chr1Base)[2] <- "Marker" 
colnames(Chr1Bgl)[2] <- "Marker" 
Chr1Imp <- merge(Chr1Base, Chr1Bgl, by="Marker")


setwd("/home/janao/Genotipi/beagle")
Chr1Imp <- read.csv("Chr1Imp.bgl", sep=" ", header=F)

#####################################################################################################################
#correct files for chromosomes where number of SNPs in imputed does not match number of SNPs in Base
####################################################################################################################
setwd("/home/janao/Genotipi/beagle/BasePlusSample/GP4_all")
#Chr19
Imp19 <- read.csv ("Chr19Imp.bgl", sep=" ", header=F)
Imp19 <- Imp19 [, 2:400]
colnames(Imp19)[1] <- "id"
Base19 <- read.csv("Chr19Base.bgl", sep=" ", header=F)
colnames(Base19)[2] <- "id"
Chr19ImpBase <- merge (Base19, Imp19, by="id", all=T, sort=F)
orderRow <- match(Base19$id, Chr19ImpBase$id)
Chr19ImpBase <- Chr19ImpBase[orderRow, ]
#Chr19ImpBase_sorted <- Chr19ImpBase[order(Base19$id),]
#Chr19ImpBase <- Chr19ImpBase_sorted[c(which(Chr19ImpBase_sorted$id=="ID"), 1:((which(Chr19ImpBase_sorted$id=="ID"))-1), (which(Chr19ImpBase_sorted$id=="ID")+1):nrow(Chr19ImpBase_sorted)), c(2,1,3:(ncol(Chr19ImpBase_sorted)))]
Chr19ImpBase <- Chr19ImpBase[,c(2,1,3:ncol(Chr19ImpBase))]
write.table(Chr19ImpBase, "Chr19ImpBase.bgl", row.names=F, col.names=F, sep=" ", quote=F, na="?")

#Chr9
Imp9 <- read.csv ("Chr9Imp.bgl", sep=" ", header=F)
Imp9 <- Imp9 [, 2:400]
colnames(Imp9)[1] <- "id"
Base9 <- read.csv("Chr9Base.bgl", sep=" ", header=F)
colnames(Base9)[2] <- "id"
Chr9ImpBase <- merge (Base9, Imp9, by="id", all=T, sort=F)
orderRow <- match(Base9$id, Chr9ImpBase$id)
Chr9ImpBase <- Chr9ImpBase[orderRow, ]
Chr9ImpBase <- Chr9ImpBase[,c(2,1,3:ncol(Chr9ImpBase))]
#Chr9ImpBase <- Chr9ImpBase[c(nrow(Chr9ImpBase), 1:(nrow(Chr9ImpBase)-1)), c(2,1,3:(ncol(Chr9ImpBase)))]
write.table(Chr9ImpBase, "Chr9ImpBase.bgl", row.names=F, col.names=F, sep=" ", quote=F, na="?")

#Chr20
Imp20 <- read.csv ("Chr20Imp.bgl", sep=" ", header=F)
Imp20 <- Imp20 [, 2:400]
colnames(Imp20)[1] <- "id"
Base20 <- read.csv("Chr20Base.bgl", sep=" ", header=F)
colnames(Base20)[2] <- "id"
Chr20ImpBase <- merge (Base20, Imp20, by="id", all=T, sort=F)
orderRow <- match(Base20$id, Chr20ImpBase$id)
Chr20ImpBase <- Chr20ImpBase[orderRow, ]
Chr20ImpBase <- Chr20ImpBase[,c(2,1,3:ncol(Chr20ImpBase))]
#Chr20ImpBase <- Chr20ImpBase[c(nrow(Chr20ImpBase), 1:(nrow(Chr20ImpBase)-1)), c(2,1,3:(ncol(Chr20ImpBase)))]
write.table(Chr20ImpBase, "Chr20ImpBase.bgl", row.names=F, col.names=F, sep=" ", quote=F, na="?")

#Chr21
Imp21 <- read.csv ("Chr21Imp.bgl", sep=" ", header=F)
Imp21 <- Imp21 [, 2:400]
colnames(Imp21)[1] <- "id"
Base21 <- read.csv("Chr21Base.bgl", sep=" ", header=F)
colnames(Base21)[2] <- "id"
Chr21ImpBase <- merge (Base21, Imp21, by="id", all=T, sort=F)
orderRow <- match(Base21$id, Chr21ImpBase$id)
Chr21ImpBase <- Chr21ImpBase[orderRow, ]
Chr21ImpBase <- Chr21ImpBase[,c(2,1,3:ncol(Chr21ImpBase))]
#Chr21ImpBase <- Chr21ImpBase[c(nrow(Chr21ImpBase), 1:(nrow(Chr21ImpBase)-1)), c(2,1,3:(ncol(Chr21ImpBase)))]
write.table(Chr21ImpBase, "Chr21ImpBase.bgl", row.names=F, col.names=F, sep=" ", quote=F, na="?")
