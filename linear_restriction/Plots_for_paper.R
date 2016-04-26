#-----------------#
# Plots for paper #
#-----------------#

# Run after running Linear_Restriction_Call

library(ggplot2)
library(reshape2)
library(plot3D)
library(colorspace)
library(stargazer)

#############################
# 2D Plot für Ast und Resig #
#############################
theme_set(theme_bw(20))

Ast <- exp(para2[1]) * d$BHD^para2[5]
Reisig <- exp(para2[2]) * d$BHD^para2[6]
tt <- as.data.frame(cbind(d$reisig, d$ast, d$BHD))
names(tt) <- c('Twigs', 'Branches', 'DBH')
tt <- melt(tt, 'DBH')
fit <- as.data.frame(cbind(Reisig, Ast, d$BHD))
names(fit) <- c('Twigs', 'Branches', 'DBH')
fit <- melt(fit, 'DBH')

P2 <- ggplot() + geom_line(data = fit, aes(x = DBH, y = value, color = variable), size = 1.2) + ylab('Biomass [kg]') + xlab('DBH [mm]') + 
      geom_point(data = tt, aes(x = DBH, y = value, color = variable, shape = variable), size = 2.5) + 
      theme(legend.title = element_blank(), legend.position = 'bottom') + 
      scale_color_manual(values = c('black', 'gray55')) +
      scale_shape_manual(values=c(8, 16)) 
P2

#########################
# QQ-Plot für residuals #
#########################

qqplot <- function(y,
                   distribution=qnorm,
                   title="Normal Q-Q",
                   xlab="Theoretical Quantiles",
                   ylab="Sample Quantiles") {
    x <- distribution(ppoints(y))
  d <- data.frame(x=x, y=sort(y))
  p <- ggplot(d, aes(x=x, y=y)) +
    geom_point() +
    geom_line(aes(x=x, y=x)) +
    #opts(title=title) +
    xlab(xlab) +
    ylab(ylab)
  return(p)
}

Residuals.A <-  sort(d$ast) - sort(Ast)
qqplot(Residuals.A)

Residuals.R <-  sort(d$reisig) - sort(Reisig)
qqplot(Residuals.R)

D.plot.boot <- function(x, y){
  z <-   exp(para2[3]) * x^para2[7] * y^para2[9] 
  return(z)
}
Derbholz <- D.plot.boot(x = d$BHD, y = d$Hoehe)
Residuals.D <- sort(d$derbholz) - sort(Derbholz)
qqplot(Residuals.D)

DR.plot.boot <- function(x, y){
  z <-   exp(para2[4]) * x^para2[8] * y^para2[10] 
  return(z)
}
Derbholzrinde <- DR.plot.boot(x = d$BHD, y = d$Hoehe)
Residuals.DR <-  sort(d$derbholzrinde) - sort(Derbholzrinde)
qqplot(Residuals.DR)

VB.plot.boot <- function(x,y){
  z <- exp(para2[3]) * x^para2[7] * y^para2[9] + exp(para2[4]) * x^para2[8] * y^para2[10] + exp(para2[1]) * x^para2[5] + exp(para2[2]) * x^para2[6]
  return(z)
}
VOllbaum <- VB.plot.boot(x = d$BHD, y = d$Hoehe)
Residuals.V <- sort(d$vollbaum) - sort(VOllbaum) 
qqplot(Residuals.V)

#############
# 3D Plots  #
#############

# Derbholz
farbe <- gray.colors(100)
z.werte.D <- outer(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), D.plot.boot)
x <- seq(min(exp(0)),700, length.out = 30)
y <- seq(min(exp(h)), max(exp(h)), length.out = 30)


pdf('Stem_wood.pdf', h = 7, w = 8)
perspbox(x = x, y = y, z = z.werte.D, bty = "b2", ticktype = "detailed", xlim = c(0, 700),  
         d = 2, plot = F, theta = 20, phi = 15, expand = 0.7, main = 'Stem wood',
         xlab = 'DBH [mm]', ylab = 'Height [cm]', zlab = 'Biomass [kg]', cex.axis = 0.7)
scatter3D(exp(bhd), exp(h), d$derbholz, colkey = F, pch = 16, add = T, xlim = c(0,700),
          surf= list(x = x, y =  y, z =  z.werte.D, border = 'black', facets = NA),
          col = farbe[1:50], ticktype = 'detailed', cex = 2, xlim = c(0, 700)
          )
dev.off()

# Derbholzrinde
farbe <- gray.colors(100)
z.werte.DR <- outer(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), DR.plot.boot)
x <- seq(min(exp(0)),700, length.out = 30)
y <- seq(min(exp(h)), max(exp(h)), length.out = 30)


pdf('Stem_bark.pdf', h = 7, w = 8)
perspbox(x = x, y = y, z = z.werte.DR, bty = "b2", ticktype = "detailed", xlim = c(0, 700),
         d = 2, plot = F, theta = 20, phi = 15, expand = 0.7, main = 'Stem bark',
         xlab = 'DBH [mm]', ylab = 'Height [cm]', zlab = 'Biomass [kg]', cex.axis = 0.7)
scatter3D(exp(bhd), exp(h), d$derbholzrinde, colkey = F, pch = 16, add = T,
          surf= list(x = x, y =  y, z =  z.werte.DR, border = 'black', facets = NA), zlim = c(0, 4500),
          col = farbe[1:50], ticktype = 'detailed', cex = 2, xlim = c(0, 700)
)
dev.off()

# Vollbaum
farbe <- gray.colors(100)
z.werte.V <- outer(seq(min(exp(bhd)),max(exp(bhd)), length.out = 30), seq(min(exp(h)), max(exp(h)), length.out = 30), VB.plot.boot)
x <- seq(min(exp(0)),700, length.out = 30)
y <- seq(min(exp(h)), max(exp(h)), length.out = 30)
fitt <- VB.plot.boot(exp(bhd),exp(h))

pdf('Tree.pdf', h = 6, w = 10)
op <- par(mar = c(0,2,0,0)+0.6)
perspbox(x = x, y = y, z = z.werte.V, bty = "b2", ticktype = "detailed", xlim = c(0, 700),
         d = 2, plot = F, theta = 20, phi = 15, expand = 0.5, 
         xlab = 'DBH [mm]', ylab = 'Height [cm]', zlab = 'Biomass [kg]', cex.axis = 0.7)
scatter3D(exp(bhd), exp(h), d$vollbaum, colkey = F, pch = 18, add = T, cex = 2.3,
          surf= list(x = x, y =  y, z =  z.werte.V, border = 'black', facets = NA, fit = fitt), zlim = c(0, 4500),
          col = farbe[1:50], ticktype = 'detailed', xlim = c(0, 700)
)
par(op)
dev.off()

############################################
# Plot zur visualisierung der Restrikition #
############################################

VB.plot.boot.L <- function(x,y){
  z <-  para2[1] + para2[2] + para2[3] + para2[4] + para2[5]*x + para2[6]*x + para2[7]*x + para2[8]*x + para2[9]*y + para2[10]*y 
  return(z)
}
Vollbaum <- rowSums(log(d[,c(6 : 9)]))
farbe <- gray.colors(100)
z.werte.VL <- outer(seq(min(bhd),max(bhd), length.out = 30), seq(min(h), max(h), length.out = 30), VB.plot.boot.L)
xl <- seq(min(bhd),max(bhd), length.out = 30)
yl <- seq(min(h), max(h), length.out = 30)

m <- matrix(1, ncol = 2, nrow = 2)

pdf('Restriction.pdf', h = 6, w = 10)
op <- par(mar = c(0,2,0,0)+0.6)
perspbox(x = xl, y = yl, z = z.werte.VL, bty = "b2", ticktype = "detailed", 
         d = 2, plot = F, theta = 20, phi = 20, expand = 0.5, 
         xlab = 'ln DBH [mm]', ylab = 'ln Height [cm]', zlab = 'Sum of ln Biomass', cex.axis = 0.7)
scatter3D(c(xl[4], xl[25], xl[8]), c(yl[7], yl[10], yl[20]), c(z.werte.VL[4,7], z.werte.VL[25,10] , z.werte.VL[8,20]), colkey = F, 
          type = 'h', add = T, pch = 16,
          surf= list(x = xl, y =  yl, z = z.werte.VL, border = 'black', facets = NA), 
          col = 'black', ticktype = 'detailed', cex = 2
)
par(op)
dev.off()

######################
# Table of covariats #
######################

Modelle <- read.csv2('coef.csv')
Modelle <- Modelle[,c(2:6)]
names(Modelle) <- c('Components', 'OLS', 'Restricted OLS', 'NSUR', 'GNLS')
stargazer(Modelle)
