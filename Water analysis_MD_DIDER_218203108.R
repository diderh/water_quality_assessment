#####Water analysis r-code
#MD. Dider Hossain
#Matriculation Number: 218203108

####################################################################
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")#for clear console
library(boot)
library(MASS)
library(rcompanion)
library(ggplot2)
library(gridExtra)
library(FSA)
setwd("C:/Desktop/Water analysis BIPU/Calculation")
wqp <- read.csv("Water quality analysis (ANOVA).csv")
wqp
wqp1 <- data.frame(wqp[,-1])
wqp1


##For Na

#Finding correlated variables
cor(wqp1)
cor(wqp1, use = "complete.obs")

##symnum()
symnum(cor(wqp1, use = "complete.obs"))

##creating the indicator variable

Ind <- function(t){
  x<- dim(length(t))
  x[which(!is.na(t))] = 1
  x[which(is.na(t))] = 0
  return(x)
}

wqp1$I <- Ind(wqp1$Na...mg.L.1.)
wqp1  


###Fitting a linear model Na on Sulphate
y <- wqp1$Na...mg.L.1.
x <- wqp1$Sulphate..mg.L.1.
lm(y ~ x, data = wqp1)
summary(lm(y ~ x, wqp1))

#y= -0.05172 + 0.07286x

###Missing value imputation

for(i in 1:nrow(wqp1)){
  if(wqp1$I[i] == 0)
  {
    wqp1$Na...mg.L.1.[i] = -0.05172 + 0.07286*wqp1$Sulphate..mg.L.1.[i]
  }
}
wqp1 ####imputing data



###############################################################################
##For Ammonia
wqp2 <- wqp1
head(wqp1)
#Finding correlated variables
cor(wqp2)
cor(wqp2, use = "complete.obs")

##symnum()
symnum(cor(wqp2, use = "complete.obs"))

##creating the indicator variable

Ind <- function(t){
  x<- dim(length(t))
  x[which(!is.na(t))] = 1
  x[which(is.na(t))] = 0
  return(x)
}

wqp2$I <- Ind(wqp2$NH4...mg.L.1.)
wqp2  


###Fitting a linear model Na on Sulphate
y <- wqp2$NH4...mg.L.1.
x <- wqp2$KH..mmol.L.1.
lm(y ~ x, data = wqp2)
summary(lm(y ~ x, wqp2))

#y= -1.157 + 1.511 x

###Missing value imputation

for(i in 1:nrow(wqp2)){
  if(wqp2$I[i] == 0)
  {
    wqp2$NH4...mg.L.1.[i] = -1.157 + 1.511*wqp2$KH..mmol.L.1.[i]
  }
}
wqp2 ####imputing data



############################################################################################

### For Calcium
wqp3<- wqp2
head(wqp3)
#Finding correlated variables
cor(wqp3)
cor(wqp3, use = "complete.obs")

##symnum()
symnum(cor(wqp3, use = "complete.obs"))

##creating the indicator variable

Ind <- function(t){
  x<- dim(length(t))
  x[which(!is.na(t))] = 1
  x[which(is.na(t))] = 0
  return(x)
}

wqp3$I <- Ind(wqp3$Ca2...mg.L.1.)
wqp3  


###Fitting a linear model Na on Sulphate
y <- wqp3$Ca2...mg.L.1.
x <- wqp3$EC..ms.m.
lm(y ~ x, data = wqp3)
summary(lm(y ~ x, wqp3))

#y= 25.9830  + 0.2746 x

###Missing value imputation

for(i in 1:nrow(wqp3)){
  if(wqp3$I[i] == 0)
  {
    wqp3$Ca2...mg.L.1.[i] = 25.9830  + 0.2746*wqp3$EC..ms.m.[i]
  }
}
wqp3 ####imputing data




########################################################################
wqp4<- wqp3
head(wqp4)
### For Magnesium
wqp4<- wa2
head(wqp4)
#Finding correlated variables
cor(wqp4)
cor(wqp4, use = "complete.obs")

##symnum()
symnum(cor(wqp4, use = "complete.obs"))

##creating the indicator variable

Ind <- function(t){
  x<- dim(length(t))
  x[which(!is.na(t))] = 1
  x[which(is.na(t))] = 0
  return(x)
}

wqp4$I <- Ind(wqp4$Mg2...mg.L.1.)
wqp4  


###Fitting a linear model Na on Sulphate
y <- wqp4$Mg2...mg.L.1.
x <- wqp4$Ca2...mg.L.1.
lm(y ~ x, data = wqp4)
summary(lm(y ~ x, wqp4))

#y= 2.2239  + 0.1016 x

###Missing value imputation

for(i in 1:nrow(wqp4)){
  if(wqp4$I[i] == 0)
  {
    wqp4$Mg2...mg.L.1.[i] = 2.2239  + 0.1016*wqp4$Ca2...mg.L.1.[i]
  }
}
wqp4 ####imputing data



######################################################################
#For Potassium
wqp5<- wqp4
head(wqp5)
#Finding correlated variables
cor(wqp5)
cor(wqp5, use = "complete.obs")

##symnum()
symnum(cor(wqp5, use = "complete.obs"))

##creating the indicator variable

Ind <- function(t){
  x<- dim(length(t))
  x[which(!is.na(t))] = 1
  x[which(is.na(t))] = 0
  return(x)
}

wqp5$I <- Ind(wqp5$K...mg.L.1.)
wqp5  


###Fitting a linear model Na on Sulphate
y <- wqp5$K...mg.L.1.
x <- wqp5$Na...mg.L.1.
lm(y ~ x, data = wqp5)
summary(lm(y ~ x, wqp5))

#y= 1.030  +  2.552 x

###Missing value imputation

for(i in 1:nrow(wqp5)){
  if(wqp5$I[i] == 0)
  {
    wqp5$K...mg.L.1.[i] = 1.030  + 2.552*wqp5$Na...mg.L.1.[i]
  }
}
wqp5 ####imputing data

### Now take the treatment name

wqp #### initial data
wqp_small = wqp[ , 2:18]
wa_treatment = wqp[,1]
head(wqp_small)
head(wa_treatment)
imputed_data = cbind(wqp5, wa_treatment)
head(imputed_data)
final_data <- imputed_data[,-18]
head(final_data)




##############################################################################################
##To see the difference between median and mean
summary(final_data$BOD5..mg.L.1.)
shapiro.test(final_data$BOD5..mg.L.1.)

##Data should be transformed
#log transformation
BOD <- final_data$BOD5..mg.L.1.
BOD_log <- log10(BOD+1)
plot(density(BOD_log))
summary(BOD_log)
shapiro.test(BOD_log)


## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$ANC..mmol.L.1.))
shapiro.test(final_data$ANC..mmol.L.1.)
plotNormalHistogram(final_data$ANC..mmol.L.1.)

##To see the difference between median and mean
summary(final_data$ANC..mmol.L.1.)
shapiro.test(final_data$ANC..mmol.L.1.)

##Data should be transformed
#log transformation
ANC <- final_data$ANC..mmol.L.1.
ANC_log <- log10(ANC+1)
plot(density(ANC_log))
summary(ANC_log)
shapiro.test(ANC_log)

###Tukey transformation
ANC.tuk <- transformTukey(ANC)

#KH
### Now check normality for individual parameter
##checking the normality of data with graph

## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$KH..mmol.L.1.))
shapiro.test(final_data$KH..mmol.L.1.)
plotNormalHistogram(final_data$KH..mmol.L.1.)

##To see the difference between median and mean
summary(final_data$KH..mmol.L.1.)
shapiro.test(final_data$KH..mmol.L.1.)

##Data should be transformed
#log transformation
KH <- final_data$KH..mmol.L.1.
KH_log <- log10(KH+1)
plot(density(KH_log))
summary(KH_log)
shapiro.test(KH_log)

###Tukey transformation
KH.tuk <- transformTukey(KH)




#Nitrate
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$Nitrate..mg.L.1.))
shapiro.test(final_data$Nitrate..mg.L.1.)
plotNormalHistogram(final_data$Nitrate..mg.L.1.)

##To see the difference between median and mean
summary(final_data$Nitrate..mg.L.1.)
shapiro.test(final_data$Nitrate..mg.L.1.)

##Data should be transformed
#log transformation
Nitrate <- final_data$Nitrate..mg.L.1.
Nitrate_log <- log10(Nitrate+1)
plot(density(Nitrate_log))
summary(Nitrate_log)
shapiro.test(Nitrate_log)

###Tukey transformation
Nitrate.tuk <- transformTukey(Nitrate)



#For sulphate
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$Sulphate..mg.L.1.))
shapiro.test(final_data$Sulphate..mg.L.1.)
plotNormalHistogram(final_data$Sulphate..mg.L.1.)

##To see the difference between median and mean
summary(final_data$Sulphate..mg.L.1.)
shapiro.test(final_data$Sulphate..mg.L.1.)



## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$Na...mg.L.1.))
shapiro.test(final_data$Na...mg.L.1.)
plotNormalHistogram(final_data$Na...mg.L.1.)

##To see the difference between median and mean
summary(final_data$Na...mg.L.1.)
shapiro.test(final_data$Na...mg.L.1.)


#pH
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$pH))
shapiro.test(final_data$pH)
plotNormalHistogram(final_data$pH)

##To see the difference between median and mean
summary(final_data$pH)
shapiro.test(final_data$pH)

#Calcium
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$Ca2...mg.L.1.))
shapiro.test(final_data$Ca2...mg.L.1.)
plotNormalHistogram(final_data$Ca2...mg.L.1.)

##To see the difference between median and mean
summary(final_data$Ca2...mg.L.1.)
shapiro.test(final_data$Ca2...mg.L.1.)



#Magnesium
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$Mg2...mg.L.1.))
shapiro.test(final_data$Mg2...mg.L.1.)
plotNormalHistogram(final_data$Mg2...mg.L.1.)

##To see the difference between median and mean
summary(final_data$Mg2...mg.L.1.)
shapiro.test(final_data$Mg2...mg.L.1.)



#EC
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$EC..ms.m.))
shapiro.test(final_data$EC..ms.m.)
plotNormalHistogram(final_data$EC..ms.m.)

##To see the difference between median and mean
summary(final_data$EC..ms.m.)
shapiro.test(final_data$EC..ms.m.)




#Mean_phosphate
## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$Mean_Phosphate))
shapiro.test(final_data$Mean_Phosphate)
plotNormalHistogram(final_data$Mean_Phosphate)

##To see the difference between median and mean
summary(final_data$Mean_Phosphate)
shapiro.test(final_data$Mean_Phosphate)

#log transformation
phosphate <- final_data$Mean_Phosphate
phosphate_log <- log10(Nitrate+1)
plot(density(phosphate_log))
summary(phosphate_log)
shapiro.test(phosphate_log)

###Tukey transformation
phosphate.tuk <- transformTukey(phosphate)

## Now check normality for individual parameter
##checking the normality of data with graph
par(mfrow = c(2,1))
plot(density(final_data$DISSOLVED.O2....))
shapiro.test(final_data$DISSOLVED.O2....)
plotNormalHistogram(final_data$DISSOLVED.O2....)

##To see the difference between median and mean
summary(final_data$DISSOLVED.O2....)
shapiro.test(final_data$DISSOLVED.O2....)




##################################################################################
###Boxplot
parametrictest <- cbind(final_data,BOD_log, ANC.tuk, KH.tuk, Nitrate.tuk, final_data$Sulphate..mg.L.1., phosphate.tuk)
head(parametrictest)
split.screen(c(2,2))
screen(1)
boxplot(BOD_log ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "BOD (mg/L)", col = "lavender")
screen(2)
boxplot(ANC.tuk ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Acidity (mmol/L)",col = "lavender")
screen(3)
boxplot(DISSOLVED.O2.... ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "DO (%)",col = "lavender")
screen(4)
boxplot(KH.tuk ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Hardness (mmol/L)",col = "lavender")
close.screen(all=TRUE)
split.screen(c(2,2))
screen(1)
boxplot(Nitrate.tuk ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Nitrate (mg/L)",col = "lavender")
screen(2)
boxplot(Sulphate..mg.L.1. ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Sulphate (mg/L)",col = "lavender")
screen(3)
boxplot(Na...mg.L.1. ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Sodium (mg/L)", col = "lavender")
screen(4)
boxplot(Ca2...mg.L.1. ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Calcium (mg/L)",col = "lavender")
close.screen(all=TRUE)
split.screen(c(2,2))
screen(1)
boxplot(Mg2...mg.L.1. ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Magnesium (mg/L)", col = "lavender")
screen(2)
boxplot(pH ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "pH",col = "lavender")
screen(3)
boxplot(EC..ms.m. ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "EC (ms/m)", col = "lavender")
screen(4)
boxplot(phosphate.tuk ~ wa_treatment, data = parametrictest, xlab = "Treatments", ylab = "Phosphate(mg/L)",col = "lavender")
close.screen(all=TRUE)









#############################################################################################################
####ANOVA and ANOVA POST HOC FOR ALL VARIABLES.

parametrictest <- cbind(final_data,BOD_log, ANC.tuk, KH.tuk, Nitrate.tuk, final_data$Sulphate..mg.L.1., phosphate.tuk)
head(parametrictest)
### For Hardness
model1 <- aov(parametrictest$KH.tuk ~ parametrictest$wa_treatment)
summary(model1)
TukeyHSD(model1)
t.test(parametrictest$KH.tuk)

### For Nitrate
model2 <- aov(parametrictest$Nitrate.tuk ~ parametrictest$wa_treatment)
summary(model2)
TukeyHSD(model2)
t.test(parametrictest$Nitrate.tuk)


### For sulphate
model3 <- aov(parametrictest$Sulphate..mg.L.1. ~ parametrictest$wa_treatment)
summary(model3)
TukeyHSD(model3)
t.test(parametrictest$Sulphate..mg.L.1.)

### For sodium
model4 <- aov(parametrictest$Na...mg.L.1. ~ parametrictest$wa_treatment)
summary(model4)
TukeyHSD(model4)
t.test(parametrictest$Na...mg.L.1.)

### For calcium
model5 <- aov(parametrictest$Ca2...mg.L.1. ~ parametrictest$wa_treatment)
summary(model5)
TukeyHSD(model5)
t.test(parametrictest$Ca2...mg.L.1.)

### For Magnesium
model6 <- aov(parametrictest$Mg2...mg.L.1. ~ parametrictest$wa_treatment)
summary(model6)
TukeyHSD(model6)
t.test(parametrictest$Mg2...mg.L.1.)

### For pH
model7 <- aov(parametrictest$pH ~ parametrictest$wa_treatment)
summary(model7)
TukeyHSD(model7)
t.test(parametrictest$pH)


### For EC
model8 <- aov(parametrictest$EC..ms.m. ~ parametrictest$wa_treatment)
summary(model8)
TukeyHSD(model8)
t.test(parametrictest$EC..ms.m.)

### For Phosphate
model9 <- aov(parametrictest$phosphate.tuk ~ parametrictest$wa_treatment)
summary(model9)
TukeyHSD(model9)
t.test(parametrictest$phosphate.tuk)

### For BOD
model10 <- aov(parametrictest$BOD_log ~ parametrictest$wa_treatment)
summary(model10)
TukeyHSD(model10)
t.test(parametrictest$BOD_log)


### For ANC
model11 <- aov(parametrictest$ANC.tuk ~ parametrictest$wa_treatment)
summary(model11)
TukeyHSD(model11)
t.test(parametrictest$ANC.tuk)



### For DO
model12 <- aov(parametrictest$DISSOLVED.O2.... ~ parametrictest$wa_treatment)
summary(model12)
TukeyHSD(model12)
t.test(parametrictest$DISSOLVED.O2....)




############################################################################################

###PCA for mine
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")#for clear console
setwd("C:/Desktop/Water analysis BIPU/Calculation")
mydata <- read.csv("PCA analysis.csv")
summary(mydata)
data<- mydata[,-1]
head(data)

##Principal component analysis
pca <- princomp(data, scores = TRUE, cor = TRUE)
summary(pca)
pca$sdev
attributes(pca)
pca$scores
v <- pca$loadings
pc1 = v[,1] 
print(pc1)
pc2= v[,2] 
print(pc2)
pc3= v[,3] 
print(pc3)
loadingforvariables <- cbind(pc1, pc2, pc3)
head(loadingforvariables)
#####Scree plot of eigenvalues
plot(pca)
screeplot(pca, type="line", main = "Scree plot")

###Biplot of score variables
biplot(pca)


####Data scaling
#PCA tries to get the features with maximum variance and the variance is high for high magnitude features
#Centering is done by subtracting column mean of data from their corresponding collumns this will result in mean deviations
# Scaling is done by dividing the centered columns of data by their standard deviations 
data.scaled <- scale(data,
                     center = TRUE,
                     scale = TRUE)
head(data.scaled)
summary(data.scaled)

####Extraction of principal components (eigenvalues)

#eigen values and eigen vectors determination
#Eigen function used to determine the eigen values and eigen vectors 
#Eigen value shows the variance in data in that direction
#First two components account most of the variation in the data
e = eigen(cov(data))
print(e)

#eigen values and eigen vectors determination for scale data to reduce the data size
e.scaled = eigen(cov(data.scaled))
print(e.scaled)
attributes(e.scaled)
#Proportion of variances

e.scaled$values            # Variance
e.scaled$values/9         # Proportion of variance
cumsum(e.scaled$values/9)  # Cumulative proportion


####Computing principal component vector
#data multiplied by eigenvectors to get the PC vectors
#%*% = R operator for matrix multiplication of scaled data matrix

data.pc = as.matrix(data.scaled) %*% e.scaled$vectors
head(data.pc)
attributes(data.pc)
#or PCA using stats package
require(stats)
pc<- prcomp(x = data,
            center = TRUE,
            scale. = TRUE)
attributes(pc)
pc$sdev
pc$rotation
head(pc$x)
summary(pc)
attributes(pc)
pc$sdev
### Display of the screeplot
#plot of the eigenvalues
#plot for the relative proportion of variances or eigenvalues explained by PCs

plot(e.scaled$values, xlab = "Index", ylab = "Lambda",
     main = "Scree plot",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.8)


#or

screeplot(x=pc, type = "line", main = "Scree plot")
plot(pc)


#### Correlation of the original variables with PCs
# Mean deviations
#mean.dev = matrix of deviations of observations from means

mean.dev <- scale(data,
                  center = TRUE,
                  scale = TRUE)
head(mean.dev) #mean deviation for each variables

##PCs for non scaled data
#m.dev.pc = matrix multiplication of mean deviations and eigen vectors
#%*% = R operator for matrix multiplication 

mean.dev.pc = mean.dev %*% e$vectors
head(mean.dev.pc)

### Correlation between PCs and original variables

r.cor = cor(cbind(mean.dev.pc, data))
print(r.cor)

#correlation of the two most important PCs and variables
# Get values from 10th to 18th variables and of 1st and 2nd row (PCs)

r1 = r.cor[10:18, 1:2] ## Corelation matrix of the first three PC components
print(r1)


# plot for the correlation of the original variables with the PCs

par(pty = "s")

ucircle = cbind(cos((0:360)/180*pi), sin((0:360)/180*pi))
plot(ucircle,
     type = "l", lty = "solid", col = "blue", lwd = 2,
     xlab = "PC1", ylab = "PC2", main = "Water quality parameters",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.8)
abline(h = 0.0, v = 0.0)
text(x = r1,
     label = c("KH", "NO3", "SO4", "Na+", "Ca2+", "Mg2+", "pH", "EC", "Phosphate"),
     cex = 0.9)
#Plotting the principal components
# Plot of the PC1 vs PC2
plot(x = data.pc[, 1], y = data.pc [, 2],
     pch = c (rep(16, 2), rep(17, 2),rep(18, 6), cex = 3),
     col = c(rep("green", 2), rep("red", 2),  rep("blue", 6)),
     xlab = "PC1", ylab = "PC2", main = "First vs Second PC",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.8)
#Add legend
legend(x = "topleft",
       legend = c("River","Effluent", "Mixing"),
       pch = c(16, 17, 18),
       col = c("green", "red", "blue"),
       text.font = 1,
       cex = 1,
       bty = "n",
       x.intersp = 0.5,
       y.intersp = 0.8,
       xpd = FALSE,
       adj = c(0, 0.25))
box()

#Add straight line
abline(v = 0, h = 0, lty = 2, col = "grey25")


#or
biplot(x = pc, choices = 1:2)


#### Plot of the second vs third pcs
plot(x = data.pc[, 2], y = data.pc [, 3],
     pch = c (rep(3, 10), rep(1, 10)),
     col = c(rep("blue", 10), rep("green", 10)),
     xlab = "PC2", ylab = "PC3", main = "Second vs Third PC",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.8)

#or
biplot(x = pc, choices = 2:3)

#### Plot of the second vs third pcs
plot(x = data.pc[, 1], y = data.pc [, 3],
     pch = c (rep(3, 10), rep(1, 10)),
     col = c(rep("blue", 10), rep("green", 10)),
     xlab = "PC1", ylab = "PC3", main = "First vs Third PC",
     cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.8)

#or
biplot(pc, choices = c(1, 3))




#############################################################################
#ggbiplot
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")#for clear console
setwd("C:/Desktop/Water analysis BIPU/Calculation")
data <- read.csv("PCA analysis.csv")
head(data)

####PCA
pc<- prcomp(x = data[,-1],
            center = TRUE,
            scale. = TRUE)
print(pc)
summary(pc)
attributes(pc)
pc$rotation
pc$loading
###Visualizing biplot
require(graphics)
par()
par(pty ="s",
    cex.main = 1.2,
    cex.lab = 1,
    font.main = 2,
    font.lab = 2,
    family = "sans",
    col.main = "gray10",
    col.lab = "gray10",
    fg = "gray10",
    las = 0)
range(pc$x[,1])
range(pc$x[,2])
plot.new()
plot.window(xlim = c(-5, 5),
            ylim = c(-2, 2),
            asp = 1)
axis(side = 1,
     at = c(-5, -2.5, 0, 2.5, 5.0),
     labels = TRUE)
axis(side = 2,
     at = c(-5, -2.5, 0, 2.5, 5.0),
     labels = TRUE)
title(main = "Biplot for PCs of water quality data",
      line = 3,
      adj = 0.5)
title(xlab = paste("PC 1 (",
                   round(summary(pc)$importance[2]*100,
                         digits = 1),
                   "%)",
                   sep = ""),
      ylab = paste("PC 2 (",
                   round(summary(pc)$importance[5]*100,
                         digits = 1),
                   "%)",
                   sep = ""),
      line = 2,
      adj = 0.5)


head(pc$x[,1:2])    

##Add points or object to the plot
points(x = pc$x[,1:2],
       pch = c(rep(16, times = 2),
               rep(17, times =2),
               rep(18, times =6)),
       cex = 1,
       col = c(rep("green", times = 2),
               rep("red", times = 2),
               rep("blue", times = 6)))

#Add ellipse to the plot
library(ellipse)
polygon(ellipse(x = cor(pc$x[5:10, 1], pc$x[5:10, 2]),
                centre = colMeans(pc$x[5:10, 1:2]),
                level = 0.85),
        border = "gray10",
        lty = "solid",
        lwd = 1.0,
        col = adjustcolor("blue", alpha.f = 0.20))

##Add second axis for variable vectors to the plot
## Allow a second plot on the same graph

par(new = TRUE, las = 1)
#xlim for PC1
range(pc$rotation[,1])
range(pc$rotation[,2])

plot.window (xlim = c(-1, 1),
             ylim = c(-1,1),
             asp = 1)
axis(side = 3,
     at = c(-1, -0.5, 0, 0.5, 1),
     labels = TRUE,
     col = "navy",
     col.ticks = NULL,
     lwd = 2,
     col.axis = "navy")

axis(side = 4,
     at = c(-1, -0.5, 0, 0.5, 1),
     labels = TRUE,
     col = "navy",
     col.ticks = NULL,
     lwd = 2,
     col.axis = "navy")

### Adding labels for second axis
mtext((text = "PC 1 rotations"),
      side = 3,
      cex = 1,
      font = 2,
      family = "sans",
      col = "gray10",
      line =2)

mtext((text = "PC 2 rotations"),
      side = 4,
      cex = 1,
      font = 2,
      family = "sans",
      col = "gray10",
      line =2,
      las = 3)

# Draw a box around the plot and straight line

box()

#Add straight line
abline(v = 0, h = 0, lty = 2, col = "grey25")

#Add variables vectors or arrows to the plot
arrows(x0 = 0, x1 = pc$rotation[,1],
       y0 = 0, y1 = pc$rotation[,2],
       col = "navy",
       length = 0.1,
       lwd = 2,
       angle = 30)

#Rotated PC1 PC2 PC3 values for variables

pc$rotation[,1]
pc$rotation[,2]
pc$rotation[,3]

#Add variable labels
text(x = pc$rotation[,1], y = pc$rotation[,2],
     labels = row.names(pc$rotation),
     cex = 0.8,
     font = 1,
     col = "gray10",
     pos = c(4, 3, 2, 1, 3, 1, 2, 3, 4, 1, 2, 4, 3)) 

### Add circle for the correlation of the original variables

ucircle = cbind(cos((0:360)/180*pi), 
                sin((0:360)/180*pi))
polygon(ucircle,
        lty= "solid", border = "gray10", lwd = 1)

#Add legend
legend(x = "topleft",
       legend = c("River","Effluent", "Mixing"),
       pch = c(16, 17, 18),
       col = c("green", "red", "blue"),
       text.font = 2,
       cex = 1.0,
       bty = "n",
       x.intersp = 0.5,
       y.intersp = 0.8,
       xpd = FALSE,
       adj = c(0, 0.25))


##########################################################################################
####Cluster analysis
###Cluster analysis with imputing data1
rm(list = ls(all = TRUE))
graphics.off()
shell("cls")#for clear console
setwd("C:/Desktop/Water analysis BIPU/Calculation")
clan1 <- read.csv("Cluster analysis with imputing data.csv")
head(clan1)
##Scatter plot
plot(pH ~ Na.,clan1)
with(clan1, text(pH ~ Na., labels=Treatment, pos=4, cex=0.7))

####Normalization of data have to be done, so that all the variables have same level playing field
a1<- clan1[, -c(10,10)]
head(a1)
b1<- apply(a1, 2, mean)# here 2 indicate the cloumn
c1<- apply(a1, 2, sd)
d1<- scale(a1,b1,c1)


##Calculating Eucladean distance
#we have 10 observation in our data set, this tables shows how close or how far each one is compared to others
# eucladean distance between river 1 and effluent3 is very high which means that this two sampling locations are very dissimilar in terms of those water quality parameters
#eucladean distance between mixing 4 and mixing 5 is very low which means that this two sampling locations are very similar in terms of those water quality parameters
distance <- dist(d1)
distance
print(distance,digits = 2)

require(graphics)
###Cluster dendogram with complete linkage
hc.c <- hclust(distance) ##, "ward.D2"
plot(hc.c, labels=clan1$Treatment, hang = -1)
rect.hclust(hc.c, k = 3, border = "red")
x<- rect.hclust(hc.c, h=6, which = c(1,2), border = 3:4)
hc.c$height

####Cluster dendogram with average linkage
hc.a <- hclust(distance, method = "average")
plot(hc.a, labels=clan1$Treatment, hang = -1)

####Cluster membership
member.c <- cutree(hc.c,3)
member.a <- cutree(hc.a,3)
table(member.c,member.a)



####Cluster means
aggregate(d1, list(member.c),mean)
aggregate(clan1[, -c(10,10)], list(member.c),mean)


####Silhouete plot
library(cluster)
plot(silhouette(cutree(hc.c,3),distance))



#######Scree plot
wss <- (nrow(d1)-1)*sum(apply(d1,2,var))
for (i in 2:8) wss[i] <- sum(kmeans(d1, centers = i)$withinss)
plot(1:8, wss, type="b", xlab = "Number of Clusters", ylab="Within group SS")

###K-Means clusterung
kc <- kmeans(d1,3)
plot(pH ~ Na.,clan1, col= kc$cluster, cex = 1.25, xlab = "Na (mg/l)", ylab = "pH")
with(clan1, text(pH ~ Na., labels=Treatment, pos=3, cex=0.5))


######################################################################################

