# Ce fichier contient le code R du projet de Séries temporelles linéaires.
# Ce projet a été réalisé par  Michael Aichoun et Jules Chapon.
# Nous avons utilisé l'indice de la fabrication de produits de voyage en France.


# Tout d'abord, nous commençons par importer les différents packages nécessaires au projet.

require(zoo)
require(tseries)
require(forecast)
library(dplyr)
library(tseries)
library(zoo)
library(ggplot2)
library(fUnitRoots)

## PARTIE I : Les données

# 1.

# Importons notre base de données

datafile <- "https://minio.lab.sspcloud.fr/juleschpn/valeurs_mensuelles.csv"
data_raw <- read.csv(datafile,sep=";") # Importe un fichier .csv dans un objet de classe data.frame

# Nous devons réaliser quelques manipulations pour la rendre utilisable

# Mise en forme

data <- data_raw[4:(nrow(data_raw)-120), c(1,2)] # Enlève les lignes et colonnes inutiles
colnames(data) <- c("Date", "Valeur") # Renomme les colonnes
data$Valeur <- as.numeric(data$Valeur) # Applique le bon format

# Indexation par chronologie croissante

data$Date <- rev(data$Date) # Met les dates dans l'ordre croissant
data$Valeur <- rev(data$Valeur) # Met les valeurs dans le bon ordre
pred <- data[(nrow(data)-1):(nrow(data)), ] # Garde en mémoire les valeurs à prédire
data <- data[1:(nrow(data)-1), ] # Enlève les 4 dernières valeurs pour faire de la prédiction

# Mise en forme des dates

data$Date = NULL
Date <- seq.Date(from = as.Date("2000-01-01"), to = as.Date("2023-01-01"), by = "month")
data <- cbind(Date = Date,data)


# Nous avons alors notre série représentant la fabrication de produits de voyage, de maroquinerie et de sellerie en France.
# Cette série commence en janvier 2000 et s'arrête en décembre 2022.
# Elle propose des valeurs mensuelles.


# 2.

plot(data,type='l')

# Cette série ne semble pas stationnaire, une tendance est visible
# Vérifions cela avec un test de Dickey-Fuller augmenté

# Définition des fonctions

Qtests1 <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

adfTest_valid <- function(series, kmax, adftype){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k," lags: residuals OK? "))
    adf <- adfTest(series, lags=k, type=adftype)
    pvals <- Qtests1(adf@test$lm$residuals, 24, fitdf = length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T)==0) {
      noautocorr <- 1; cat("OK \n")
    } else cat("nope \n")
    k <- k+1
  }
  return(adf)
}

# Nous testons sur notre série

summary(lm(data$Valeur ~ data$Date))
adf <- adfTest_valid(data$Valeur,24,adftype="ct")
adf

# La p-value est très élévée : la série n'est pas stationnaire

# Nous regardons alors les autocorrélogrammes

acf(data$Valeur,24) ; pacf(data$Valeur,24)

# Nous pouvons voir qu'il y a une autocorrélation partielle significative à l'ordre 1.

# Nous allons donc essayer de la rendre stationnaire en la différenciant à l'ordre 1

dif <- diff(data$Valeur,lag = 1)

# Nous effectuons à nouveau les tests

summary(lm(dif ~ Date[-1]))
adf <- adfTest_valid(dif,24,adftype="nc")
adf

# Ainsi la série est bien devenue stationnaire


# 3.

# Nous pouvons alors représenter graphiquement nos deux séries

dif <- as.data.frame(dif)
Dates <- seq.Date(from = as.Date("2000-02-01"), to = as.Date("2023-01-01"), by = "month")
dif <- cbind(Date = Dates,dif)

par(mfrow=c(2,1),mar=c(1, 6, 0.5, 1) + 1)
plot(data,type='l')
plot(dif,type='l')


## Partie II : Modèle ARMA


# 4.

# Pour déterminer les valeurs de p et q, regardons les autocorrélogrammes

acf(dif$dif,24); pacf(dif$dif,24)

# Nous regardons les dernières valeurs pour lesquels l'ACF et le PACF sont significatifs
# Nous avons donc p = 4 et q = 1

p_max <- 4
q_max <- 1

# Nous choisissons donc un modèle ARMA(4,1)

# Nous pouvons estimer les paramètres de ce modèle

arma41 <- arima(dif$dif, c(4,0,1))
arma41

# Nous pouvons alors tester la validité du modèle via des tests de Ljung-Box

Qtests2 <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests2(arma41$residuals, 24, 5) # Tests de LB pour les ordres 1 à 24

# L'absence d'autocorrélation n'est jamais rejetée à 95% jusqu'à 24 retards
# Le modèle est donc valide

# Voyons désormais s'il est bien ajusté

# Fonction de test des significativités individuelles des coefficients

signif <- function(estim){ 
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

signif(arma41) # Tests de siginificativité de l’ARMA(4,1)

# Nous regardons les coefficients des retards les plus élevés AR(4) et MA(1) 
# Ils ne rejettent chacun pas leur nullité à 95% (p-value > 0.05)
# Le modèle ARMA(4,1) est donc mal ajusté

# Nous allons donc tester les différents sous-modèles candidats
# Ce sont les modèles ARMA(1,0), ARMA(2,0), ARMA(3,0), ARMA(0,1), ARMA(1,1), ARMA(2,1) et ARMA(3,1)

# Fonction d’affichage des tests pour la sélection du modèle ARMA

arimafit <- function(estim){
  adjust <- round(signif(estim),3)
  pvals <- Qtests2(estim$residuals,24,length(estim$coef)-1)
  pvals <- matrix(apply(matrix(1:24,nrow=6),2,function(c) round(pvals[c,],3)),nrow=6)
  colnames(pvals) <- rep(c("lag", "pval"),4)
  cat("tests de nullité des coefficients :\n")
  print(adjust)
  cat("\n tests d'absence d'autocorrélation des résidus : \n")
  print(pvals)
}

# ARMA(1,0)

estim <- arima(dif$dif, c(1,0,0)); arimafit(estim)
arma10 <- estim

# Le modèle est bien ajusté mais n'est pas valide

# ARMA(2,0)

estim <- arima(dif$dif, c(2,0,0)); arimafit(estim)
arma20 <- estim

# Le modèle est bien ajusté mais n'est pas valide

# ARMA(3,0)

estim <- arima(dif$dif, c(3,0,0)); arimafit(estim)
arma30 <- estim

# Le modèle est bien ajusté mais n'est pas valide

# ARMA(4,0)

estim <- arima(dif$dif, c(4,0,0)); arimafit(estim)
arma40 <- estim

# Le modèle est valide et bien ajusté
# Nous le gardons comme candidat

# ARMA(0,1)

estim <- arima(dif$dif, c(0,0,1)); arimafit(estim)
arma01 <- estim

# Le modèle est valide et bien ajusté
# Nous le gardons comme candidat

# ARMA(1,1)

estim <- arima(dif$dif, c(4,0,1)); arimafit(estim)
arma11 <- estim

# Le modèle est valide mais n'est pas bien ajusté

# ARMA(2,1)

estim <- arima(dif$dif, c(2,0,1)); arimafit(estim)
arma21 <- estim

# Le modèle est valide mais n'est pas bien ajusté

# ARMA(3,1)

estim <- arima(dif$dif, c(3,0,1)); arimafit(estim)
arma31 <- estim

# Le modèle est valide mais n'est pas bien ajusté

# Nous avons donc deux modèles bien ajustés et valides :ARMA(4,0) et ARMA(0,1)
# Nous regardons les critères d'information associés (AIC et BIC)

models <- c("arma40","arma01"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))

# Nous pouvons voir que notre ARMA(0,1) minimise bien les 2 critères AIC et BIC
# Il est donc le meilleur selon ces critères

# Estimation de notre ARMA(0,1)

arma01

# Création de séries où à chaque colonne sera assignée la prédiction par un modèle

models <- c("arma01")
preds <- zoo(matrix(NA,ncol=1,nrow=2),order.by=(pred$Date))
colnames(preds) <- models

# Prédiction de dif et data par notre modèle

predarma01 <- zoo(predict(get("arma01"),2)$pred, order.by=(pred$Date))
pred1 <- data.frame(Date = pred$Date, arma01 = predarma01)
pred2 <- data.frame(pred1)

for (m in models){
  pred2[1,m] <- pred1[1,m]+tail(data$Valeur,1)[1]
  pred2[2,m] <- pred1[2,m]+pred2[1,m]
}

obs <- data.frame(pred) # 2 dernières observations de la série originelle

y <- merge(pred,pred2, by = "Date") # Affichage des observations et des prédictions
print(y)

# On calcule le RMSE (normalisé par l'écart-type de la série)

apply(y[c("arma01")],2, function(x) sqrt(sum((x-obs$Valeur)^2)/2)/sd(data$Valeur))

# Nous pouvons afficher les résidus de notre modèle

plot(arma01$residuals)
acf(arma01$residuals)
pacf(arma01$residuals)

hist(arma01$residuals)
checkresiduals(arma01)

# Les résidus ne sont pas corrélés, comme nous l'avions montré précédemment
# De plus, ils semblent être gaussiens
# Ils se comportent comme un bruit blanc

# De plus, notre modèle est un ARMA(0,1)
# Il est par conséquent causal