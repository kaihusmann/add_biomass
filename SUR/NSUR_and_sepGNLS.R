#------#
# NSUR # #GITHUB TEST
#------#

rm(list = ls())

library(systemfit)
library(nlme)

#setwd('H:\\Lange')
setwd('/home/khusmann/mnt/H/additive_Biomassefunktionen/')

d <- read.csv2("kompartimentsgewichte.csv")
coef_matrix_dh <- read.csv2('coefficients_dbh_h.csv')
coef_matrix_d  <- read.csv2('coefficients_dbh.csv')



## gnls Regression ##
# Erstellen der separaten Nichtlin. Modelle und Erstellen der Vektoren fuer Gewichtung.
gnls.derbholz <- gnls(derbholz ~ exp(b11) * (BHD ^ b12) * (Hoehe ^ b13), start = list(b11 = -10 ,b12 = 2, b13 = 1), weights=varPower(), data=d)
d$w_derbholz <- attr(gnls.derbholz$modelStruct[['varStruct']], 'weights')
write.csv2(summary(gnls.derbholz)$varBeta, './SUR/var_covar/covar_derbholz_sepGNLS.csv')


gnls.derbholzrinde <- gnls(derbholzrinde ~ exp(b11) * (BHD ^ b12) * (Hoehe ^ b13), start = list(b11 = -10 ,b12 = 2, b13 = 1), weights=varPower(), data=d)
d$w_derbholzrinde <- attr(gnls.derbholzrinde$modelStruct[['varStruct']], 'weights')
write.csv2(summary(gnls.derbholzrinde)$varBeta, './SUR/var_covar/covar_derbholzrinde_sepGNLS.csv')

gnls.ast <- gnls(ast ~ exp(b11) * (BHD ^ b12), start = list(b11 = -10 ,b12 = 1.5), weights=varPower(), data=d)
d$w_ast <- attr(gnls.ast$modelStruct[['varStruct']], 'weights')
write.csv2(summary(gnls.ast)$varBeta, './SUR/var_covar/covar_ast_sepGNLS.csv')

gnls.reisig <- gnls(reisig ~ exp(b11) * (BHD ^ b12), start = list(b11 = -10 ,b12 = 1.5), weights=varPower(), data=d)
d$w_reisig <- attr(gnls.reisig$modelStruct[['varStruct']], 'weights')
write.csv2(summary(gnls.reisig)$varBeta, './SUR/var_covar/covar_reisig_sepGNLS.csv')

gnls.vollbaum <- gnls(vollbaum ~ exp(b11) * (BHD ^ b12), start = list(b11 = -7 ,b12 = 2), weights=varPower(), data=d)
d$w_vollbaum <- attr(gnls.vollbaum$modelStruct[['varStruct']], 'weights')



## NSUR Aufruf Vorbereitung: Funktionen und Startwerte definieren ##
f.derbholz      <- w_derbholz * derbholz           ~ (exp(b11) * (BHD ^ b12) * (Hoehe ^ b13)) * w_derbholz
f.derbholzrinde <- w_derbholzrinde * derbholzrinde ~ (exp(b21) * (BHD ^ b22) * (Hoehe ^ b23)) * w_derbholzrinde
f.ast           <- w_ast * ast                     ~ (exp(b31) * (BHD ^ b32)) * w_ast
f.reisig        <- w_reisig * reisig               ~ (exp(b41) * (BHD ^ b42)) * w_reisig

f.vollbaum      <- w_vollbaum * vollbaum           ~ (exp(b11) * (BHD ^ b12) * (Hoehe ^ b13) + exp(b21) * (BHD ^ b22) * (Hoehe ^ b23) + exp(b31) * (BHD ^ b32) + exp(b41) * (BHD ^ b42)) * w_vollbaum

model <- list(f.derbholz,
              f.derbholzrinde,
              f.ast,
              f.reisig,
              f.vollbaum
)

sv <- c(b11 = -14.2, b12 = 2.2, b13 = 1.0,
        b21 = -14.8, b22 = 2.0, b23 = 0.9,
        b31 = -4.8, b32 = 1.57,
        b41 = -7.3, b42 = 1.7)

## NSUR Funktion Aufrufen ##
OLS <- nlsystemfit(method = 'OLS', eqns = model, startvals = sv, data = d)
summary(OLS)
OLS$b

SUR <- nlsystemfit('SUR', model, sv, data = d)
summary(SUR)

# Comparison: OLS, gnls, SUR
cbind(
  OLS$b,
  SUR$b,
  c(coef(gnls.derbholz), coef(gnls.derbholzrinde), coef(gnls.ast), coef(gnls.reisig))
)

b <-SUR$b
SUR$se # Std. error
SUR$covb # Covar
sqrt(SUR$covb)


## Datenausgabe ##
ret <- data.frame(component = c(
    'Derbholz.b0', 'Derbholz.b1', 'Derbholz.b2', 
    'Derbholzrinde.b0', 'Derbholzrinde.b1', 'Derbholzrinde.b2', 
    'Ast.b0', 'Ast.b1', 
    'Reisig.b0', 'Reisig.b1'
  ),
  parameter.sur = SUR$b,
  parameter.separateGNLS = c(
    coef(gnls.derbholz), 
    coef(gnls.derbholzrinde),
    coef(gnls.ast),
    coef(gnls.reisig)
    )
  )

write.csv2(ret, '/home/khusmann/mnt/H/additive_Biomassefunktionen/SUR/coef_SUR_sepGNLS.csv', row.names = FALSE)
write.csv2(SUR$covb,'./SUR/var_covar/covar_SUR.csv')

save(file = './SUR/nonlin_models.RData', list = c('SUR', 'gnls.derbholz', 'gnls.derbholzrinde', 'gnls.ast', 'gnls.reisig'))
## plots for visual validation ##

# Fitted and observed over bbh
plot(d$derbholz ~ d$BHD)
fitted.derbholz <- exp(b[1]) * (d$BHD  ^ b[2]) * d$Hoehe ^ b[3]
points(fitted.derbholz ~ d$BHD, pch = 3)

plot(d$derbholzrinde ~ d$BHD)
fitted.derbholzrinde <- exp(b[4]) * (d$BHD  ^ b[5]) * d$Hoehe ^ b[6]
points(fitted.derbholzrinde ~ d$BHD, pch = 3)

plot(d$ast ~ d$BHD)
fitted.ast <- exp(b[7]) * (d$BHD ^ b[8])
points(fitted.ast ~ d$BHD, pch = 3)

plot(d$reisig ~ d$BHD)
fitted.reisig <- exp(b[9]) * (d$BHD ^ b[10])
points(fitted.reisig ~ d$BHD, pch = 3)

plot(d$vollbaum ~ d$BHD)
fitted.vollbaum <- apply(cbind(fitted.derbholz, fitted.derbholzrinde, fitted.ast, fitted.reisig), 1, sum)
points(fitted.vollbaum ~ d$BHD, pch = 3)
