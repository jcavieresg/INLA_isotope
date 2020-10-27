rm(list=ls(all=TRUE))

#====================================================================================================================================
#                                         Proyecto IFOP 2020 - lobos marinos
#====================================================================================================================================

setwd("C:/Users/Usuario/Desktop/Projects/INLA_isotope/data")

# Packages
library(ggplot2)
library(reshape2)
library(RANN)
library(tidyverse)
library(gridExtra)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(abind)
library(ggmap)
library(dplyr)
library(ggpubr)


# load data 
dat = read.csv("data_total.csv", header=T)
dim(dat)
head(dat)
str(dat)

options(scipen=99)

#================================================================================================
# OBSERVACIONES:
# Variable a modelar 'N' (en primera instancia)
#================================================================================================


# Creamos una nueva variable para el total de lobos capturados incidentalmente
dat$N = as.numeric(dat$N)

# Transformar NA a 0 en variable respuesta
dat$N[is.na(dat$N)] = 0

# Vemos el histograma
table(dat$N)
hist(dat$N, breaks = 35)

# Otras variables disponibles
dat$T      # Temperatura
dat$Sal    # Salinidad
dat$Chla   # Clorofila
dat$OD     # Oxigeno disuelto


# Traspasamos a numericas las variables
dat$T = as.numeric(dat$T)
dat$Sal = as.numeric(dat$Sal)
dat$Chla =as.numeric(dat$Chla)
dat$OD = as.numeric(dat$OD)

# Factors la temporada
dat$Season = as.factor(dat$Season)

# # Vamos a imputar los datos de las profundidades por sus mediandas para
# # así no perder información.
# # Imputación de datos para la PROFUNDIAD_PROM
# dat$prof_min[is.na(dat$prof_min)] = median(dat$prof_min, na.rm=TRUE)
# dat$prof_max[is.na(dat$prof_max)] = median(dat$prof_max, na.rm=TRUE)
# 
# 
# # Creamos una variable de profundidad promedio
# dat$prof_med = (dat$prof_max + dat$prof_min) / 2

# Creamos otra variable "TRIM"
# library(varhandle)
# dat$trim = ifelse(dat$mes_l < 4 ,1,ifelse(dat$mes_l < 7, 2, ifelse(dat$mes_l < 10, 3, 4)))
# table(dat$trim)
# 
# # Reemplazamos las celdas NA en las capturas de aves con 0
# dat[, c(21:65)][is.na(dat[, c(21:65)])] = 0
# dat[, c(21:65)]


summary(dat$Lat)  # Celdas con NA
summary(dat$Lon)  # Celdas con NA

# Revisamos si existen NA en las variables relevantes
glimpse(dat)


#=======================================
# Variable 
#=======================================
# De acuerdo a lo conversado con Marcelo, se debería utilizar la variable de mayor
# resolución, por tanto consideramos la variable 'ESPECIE_OBJETIVO_LANCE'
table(dat$T, exclude = NULL)

#=======================================
# Variable 
#=======================================
table(dat$Sal, exclude = NULL)


#=======================================
# Variable 
#=======================================
table(dat$Chla, exclude = NULL)



#=======================================
# Variable 
#=======================================
table(dat$OD, exclude = NULL) # Tenempos 134 NA y los transformaremos a otro nivel
                                      # como S/I (sin información)

 
# #=======================================
# # Variable 
# #=======================================
# table(dat$tsm, exclude = NULL) # Tenempos 602 NA y los imputaremos con las medianas
# 
# dat$tsm[is.na(dat$tsm)] = median(dat$tsm, na.rm=TRUE)  
# summary(dat$tsm)

#================================
#        Análisis
#================================
# Por comodidad vamos a renombrar algunas variables del data.frame

names(dat)[names(dat) == "T"] <- "temp"
names(dat)[names(dat) == "Sal"] <- "sal"
names(dat)[names(dat) == "Chla"] <- "cloro"
names(dat)[names(dat) == "OD"] <- "oxi"


glimpse(dat)



#==========================================================
# Escribimos los datos que vamos a utilizar para modelar
#write.csv(dat ,"dat_fit.csv")




#==========================================================
#                 Codigo con librería mgcv
#==========================================================
# Determinación del mejor sub-conjunto de datos para el Nitrogeno (N)
library(ggplot2)
library(pscl)
library(boot)
library(MASS)
library(dplyr)
library(mgcv)
library(mgcViz)


# Leemos la nueva base
dat_glm = read.csv("dat_fit.csv", header=T)
dim(dat_glm)
head(dat_glm)
glimpse(dat_glm)

# Hacemos una copia para crear modelos con gam
dat_gam = dat_glm

dat_glm = dat_glm[, c(-1, -4, -5, -9)] # Sin variable X, Lat y Lon 
head(dat_glm)

str(dat_glm)

# Declaramos como factor al año (Year)
dat_glm$Year = as.factor(dat_glm$Year)
str(dat_glm)


# Corremos el modelo nulo
full = glm(N~., data = dat_glm, family = gaussian)
step = stepAIC(full, trace = FALSE)
step$anova


backward  = stepAIC(glm(N~.,data = dat_glm, family = gaussian),direction="backward")
forward   = stepAIC(glm(N~.,data = dat_glm, family = gaussian),direction="forward")
both      = stepAIC(glm(N~.,data = dat_glm, family = gaussian),direction="both")


backward$anova
forward$anova
both$anova

# Comparar devianzas
step$deviance
backward$deviance
forward$deviance
both$deviance   # Tienen igual devianza así que podemos ocupar cualquier

# Menos AIC es forward
formula(step)
formula(backward)
formula(forward)
formula(both)

# MEJOR PREDICTOR LINEAL PARA N
# N ~ Year + Season + Dist + Depth + cloro + oxi + temp + sal


# Vemos modelos con GAM
library(mgcv)


m0 <- gam(N ~ Year + Season + Dist + Depth + cloro + oxi + temp + sal, 
          data = dat_gam,
          family = gaussian, method = "REML")


m1 <- gam(N ~ Year + Season + Dist + Depth + cloro + oxi + temp + sal + 
          s(Lon) + s(Lat), 
          data = dat_gam,
          family = gaussian, method = "REML")



m2 <- gam(N ~ Year + Season + Dist + Depth + cloro + oxi + temp + sal + 
          s(Lon, Lat), 
          data = dat_gam, 
          family = gaussian, method = "REML")


m3 <- gam(N ~ Year + Season + Dist + Depth + cloro + oxi + temp + sal + 
          s(Lon,Lat, bs="tp", k=100), 
          data = dat_gam, 
          family = gaussian, method = "REML")


m4 <- gam(N ~ Year + Season + Dist + Depth + cloro + oxi + temp + sal + 
          s(Lon,Lat, bs="gp", k=100), 
          data = dat_gam, 
          family = gaussian, method = "REML")

AIC(m0, m1, m2, m3, m4)
anova(m0, m1, m2, m3, m4, test="Chisq")
gam.check(m2, k.rep=1000)


summary(m2)
coef(m2)


library(mgcViz)
m2 <- getViz(m2)
print(plot(m2) + labs(title = NULL), pages = 1)

print(plot(m2, select = 2:7), pages = 1)

# print(plot(m2), ask = FALSE)
# print(plot(m2, allTerms = TRUE), pages = 1)


vis.gam(x = m2, # GAM object
        view = c("Lon","Lat"), # variables
        plot.type = "persp") # kind of plot

qq(m2, method = "simul1", a.qqpoi = list("shape" = 1), a.ablin = list("linetype" = 2))
qq(m2, rep = 20, showReps = T, CI = "none", a.qqpoi = list("shape" = 19), a.replin = list("alpha" = 0.2))


o <- qq(m3, rep = 10, method = "simul1", CI = "normal", showReps = TRUE,
        ngr = 1e2, a.replin = list(alpha = 0.1), a.qqpoi = list(shape = 19))
o 

gridPrint(o, zoom(o, xlim = c(2, 2.5), ylim = c(2, 2.5)), ncol = 2)

vis.gam(x = m2, # GAM object
        view = c("Lon","Lat"), # variables
        plot.type = "contour", main='Efecto espacial',
        color = "heat") # kind of plot

x11()
plotRGL(sm(m2, 1), fix = c("z" = 0), residuals = TRUE) # Plot 3D




#===================================================================
#                        Modelo con INLA
#===================================================================
library(INLA)
dat_inla = dat_gam

# Defino las coordenadas de estudio
coo = cbind(dat_inla$Lon, dat_inla$Lat)

# Defino la grilla
m1 <- inla.mesh.2d(coo, max.edge = c(0.5, 0.5)) 
m2 <- inla.mesh.2d(coo, max.edge = c(0.5, 0.5), cutoff = 0.1) 
m3 <- inla.mesh.2d(coo, max.edge = c(0.1, 0.5), cutoff = 0.1) 
#m4 <- inla.mesh.2d(coo, max.edge = c(0.1, 0.5), offset = c(0, -0.65)) 
m5 <- inla.mesh.2d(coo, max.edge = c(0.3, 0.5),offset = c(0.03, 0.5)) 
m6 <- inla.mesh.2d(coo, max.edge = c(0.3, 0.5),offset = c(0.03, 0.5), cutoff = 0.1)
m7 <- inla.mesh.2d(coo, max.edge = c(0.3, 0.5), n = 5, offset = c(0.05, 0.1)) 
m8 <- inla.mesh.2d(coo, max.edge = c(0.3, 0.5), n = 7, offset = c(0.01, 0.3)) 
m9 <- inla.mesh.2d(coo, max.edge = c(0.3, 0.5), n = 4, offset = c(0.05, 0.3))

plot(m1)
plot(m2)
plot(m3)
plot(m5)
plot(m6)
plot(m7)
plot(m8)
plot(m9)


c(m1$n, m2$n, m3$n, m5$n, m6$n, m7$n, m8$n, m9$n)


# Boundaries
bound1 <- inla.nonconvex.hull(coo)
bound2 <- inla.nonconvex.hull(coo, convex = 0.5, concave = -0.15)
bound3 <- inla.nonconvex.hull(coo, concave = 0.5)


# Meshes
m10 <- inla.mesh.2d(boundary = bound1, cutoff = 0.05, max.edge = c(0.1, 0.2))
m11 <- inla.mesh.2d(boundary = bound2, cutoff = 0.05, max.edge = c(0.1, 0.2))
m12 <- inla.mesh.2d(boundary = bound3, cutoff = 0.05, max.edge = c(0.1, 0.2))


plot(m10)
plot(m11)
plot(m12)


# Vanos a elegir m3 ya que la triangulación es homogenea y el n no es tan grande
meshbuilder()

## Build boundary information:
## (fmesher supports SpatialPolygons, but this app is not (yet) intelligent enough for that.)
boundary.loc <- SpatialPoints(as.matrix(coo))
boundary <- list(
        inla.nonconvex.hull(coordinates(boundary.loc), 0.1),
        inla.nonconvex.hull(coordinates(boundary.loc), 0.3))

## Build the mesh:
mesh <- inla.mesh.2d(boundary=boundary,
                     max.edge=c(0.05, 0.99),
                     min.angle=c(30, 21),
                     max.n=c(48000, 16000), ## Safeguard against large meshes.
                     max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                     cutoff=0.01, ## Filter away adjacent points.
                     offset=c(0.1, 0.3)) ## Offset for extra boundaries, if needed.

## Plot the mesh:
plot(m3)
points(coo, pch = 19, col = "red")


pcprec <- list(prior = 'pcprec', param = c(1, 0.01))
A <- inla.spde.make.A(m3, loc = coo)

spde <- inla.spde2.pcmatern(mesh = m3,
                            prior.range = c(0.05, 0.01), # P(practic.range < 0.05) = 0.01
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01



stk.dat <- inla.stack(data = list(y = dat_inla$N), 
                                  A = list(A,1),
                                  effects = list(list(s = 1:spde$n.spde), 
                                  data.frame(Intercept = 1, 
                                             Year = dat_inla$Year,
                                             Season = dat_inla$Season,
                                             Dist = dat_inla$Dist,
                                             Depth = dat_inla$Depth,
                                             cloro = dat_inla$cloro,
                                             oxi = dat_inla$oxi,
                                             temp = dat_inla$temp,
                                             sal = dat_inla$sal)),
                                  tag = 'dat') 




# Formulas
formula  = y ~ 0 + Intercept + Year + Season + Dist + Depth + cloro + oxi + temp + sal
formula2 = y ~ 0 + Intercept + Year + Season + f(Dist, model = 'rw1', scale.model = TRUE) + cloro + oxi + temp + sal
formula3 = y ~ 0 + Intercept + Year + Season + Dist + Depth + cloro + oxi + temp + sal + f(s, model = spde)
formula4 = y ~ 0 + Intercept + Year + Season + f(Dist, model = 'rw1', scale.model = TRUE) + cloro + oxi + temp + sal + f(s, model = spde)



# Controles
ini.zb <- c(-1.147,-5.939,5.069)
ini.zb2 <- c(-1.147,-5.939)
ini.zb4 <- c(-1.147,-5.939, 5.069, 5.069)
cres = list(return.marginals.predictor = FALSE, return.marginals.random = FALSE)
cinla <- list(strategy = 'gaussian', int.strategy = 'eb')  # Estrategias de estimaci?n


# Modelo 1
result <- inla(formula, 
               family="gaussian",
               data=inla.stack.data(stk.dat), 
               control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE, link = 1),
               control.results = cres, control.inla = cinla,
               control.mode = list(theta = ini.zb, restart = TRUE),
               verbose=TRUE)


# Modelo 2
result2 <- inla(formula2, 
               family="gaussian",
               data=inla.stack.data(stk.dat), 
               control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE, link = 1),
               control.results = cres, control.inla = cinla,
               control.mode = list(theta = ini.zb2, restart = TRUE),
               verbose=TRUE)


# Modelo 3
result3 <- inla(formula3, 
               family="gaussian",
               data=inla.stack.data(stk.dat), 
               control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE, link = 1),
               control.results = cres, control.inla = cinla,
               control.mode = list(theta = ini.zb, restart = TRUE),
               verbose=TRUE)



# Modelo 4
result4 <- inla(formula4, 
               family="gaussian",
               data=inla.stack.data(stk.dat), 
               control.compute=list(return.marginals = TRUE, dic=TRUE, waic = TRUE, cpo = TRUE),
               control.predictor=list(A=inla.stack.A(stk.dat), compute=TRUE, link = 1),
               control.results = cres, control.inla = cinla,
               control.mode = list(theta = ini.zb4, restart = TRUE),
               verbose=TRUE)


slcpo <- function(m, na.rm = TRUE) {
        - sum(log(m$cpo$cpo), na.rm = na.rm)
}

c(result = slcpo(result), 
  result2 = slcpo(result2),
  result3 = slcpo(result3),
  result4 = slcpo(result4))


result$waic$waic
result2$waic$waic
result3$waic$waic
result4$waic$waic





# Mejor modelo es result4
round(result4$summary.fixed, 2)

# Desviacion estandar y practical range del spatial random field
round(result4$summary.hyperpar[-1, ], 3)



index <- inla.stack.index(stk.dat, tag = "dat")$data

pred_mean = result4$summary.fitted.values[index, "mean"]
pred_ll = result4$summary.fitted.values[index, "0.025quant"]
pred_ul = result4$summary.fitted.values[index, "0.975quant"]

#coords = cbind(dat_fit$LON2, dat_fit$LAT2)

dpm <- rbind(data.frame(Longitud = coo[, 1], Latitud = coo[, 2], value = pred_mean, variable = "media_pred"),
             data.frame(Longitud = coo[, 1], Latitud = coo[, 2], value = pred_ll, variable = "sup_media_pred"),
             data.frame(Longitud = coo[, 1], Latitud = coo[, 2], value = pred_ul, variable = "inf_media_pred"))
dpm$variable <- as.factor(dpm$variable)

dpm$value

ggplot(dpm) + geom_tile(aes(Longitud, Latitud, fill = value)) +
        facet_wrap(~variable, nrow = 1) +
        coord_fixed(ratio = 1) +
        scale_fill_gradient(
                name = "N medido",
                low = "blue", high = "orange"
        ) +
        theme_bw()


summary(pred_mean)
summary(pred_ll)
summary(pred_ul)



# Media posterior para varianza y el rango
r.f <- inla.spde2.result(result4, 's', spde, do.transf=TRUE)   # "i.z" es el indice espacial declarado en data.stack

# # Media marginal para la varianza
# inla.emarginal(function(x) x, r.f$marginals.log.variance.nominal[[1]])
# 
# # Media marginal para el rango
# inla.emarginal(function(x) x, r.f$marginals.range.nominal[[1]])


# Grafico para el intercepto y los efectos fijos
# x11()
# par(mfrow=c(1,2), mar=c(3,3.5,2,2), mgp=c(1.5, .5, 0), las=0)
# plot(result4$summary.random[[1]][,1:2], type='l', xlab='ano', ylab='Effect'); abline(h=0, lty=3)
# plot(result4$summary.random[[2]][,1:2], type='l', xlab='mode', ylab='Effect'); abline(h=0, lty=3)

par(mfcol=c(1,2), mar=c(4,5,3,3), las=1)
plot(result4$marginals.fix[[1]], type='l', xlab=expression(Intercept), ylab="Density", cex.lab=1.6, cex.axis=1.4)
#abline(v=alpha.c, col=4)

plot(result4$marginals.hy[[2]], type='l',xlab=expression(sigma[x]), ylab="Density", cex.lab=1.6, cex.axis=1.4)
abline(v=exp(result4$internal.summary.hyperpar[3,1]), col=4)

plot(result4$marginals.hy[[3]], type='l',xlab=expression(Range (km)), ylab="Density", cex.lab=1, cex.axis=1)
abline(v=exp(result4$internal.summary.hyperpar[1,1]), col=4)

plot(result4$marginals.hy[[4]], type='l',xlab=expression(sigma[i.z]^2), ylab="Density", cex.lab=1, cex.axis=1)
abline(v=exp(result4$internal.summary.hyperpar[2,1]), col=4)



# Prediccion del spatial random field
stepsize <- 4 * 1 / 111

x.range <- diff(range(coo[, 1]))
y.range <- diff(range(coo[, 2]))
nxy <- round(c(x.range, y.range) / stepsize)


projgrid <- inla.mesh.projector(m3, xlim = range(coo[, 1]), 
                                ylim = range(coo[, 2]), dims = nxy)


xmean <- inla.mesh.project(projgrid, result4$summary.random$s$mean)
xsd <- inla.mesh.project(projgrid, result4$summary.random$s$sd)


x11()
par(mfrow=c(1,2), mar=c(5,4,3,7))
image.plot(x=projgrid$x, y=projgrid$y, z=xmean, asp=1,xlab='Longitud', ylab='Latitud', axes=T, cex.axis=0.9, axis.args = list(cex.axis = 0.9))
bnd <- inla.mesh.boundary(m3)
title(main="Media posterior del GMRF")
#plot(norte, add = T, col = 'darkgrey')
# inter <- inla.mesh.interior(m3)
plot(norte, draw.segments=FALSE, main = '', add = T)
points(coo, pch = 19, col = "black")
#lines(inter[[1]], col=1, lwd=2)

#plot(Omega.SP3[[1]], add=T, col='darkgrey')

image.plot(x=projgrid$x, y=projgrid$y, z=xsd, asp=1,xlab='Longitud', ylab='', axes=T, cex.axis=0.9, axis.args = list(cex.axis = 0.9))
plot(norte, draw.segments=FALSE, main = '', add = T)
points(coo, pch = 19, col = "black")
title(main="Desviación estandar posterior del GMRF")
