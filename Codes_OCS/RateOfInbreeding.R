y = log(1 - AvgInbPerGen)
t = Gen - min(Gen) + 1
fit = MASS:::rlm(y ~ t, maxit=2000)
# lahko bi uporabilli tudi lm(), ampak je rlm bolj robusten (pri rastlinskih simulacijah vidim precejsnje spremembe skozi cas, mislim, da pri tebi tukaj ne bi smelo biti vecjih razlik med rlm in lm)
dF = 1 - exp(coef(fit)[2])
Ne = 1 / (2 * dF)