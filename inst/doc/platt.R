### R code from vignette source 'platt.Rnw'

###################################################
### code chunk number 1: platt.Rnw:40-44
###################################################
options(width=75)
set.seed(0)
library(platt)
plattVersion <- packageDescription("platt")$Version


###################################################
### code chunk number 2: platt.Rnw:105-106
###################################################
library(platt)


###################################################
### code chunk number 3: platt.Rnw:115-118
###################################################
library(platt)
data(MMP)
head(MMP)


###################################################
### code chunk number 4: platt.Rnw:124-125
###################################################
plattScalingResult <- plattScaling(MMP$SVM,MMP$target)


###################################################
### code chunk number 5: platt.Rnw:131-132
###################################################
str(plattScalingResult)


###################################################
### code chunk number 6: platt.Rnw:140-151
###################################################
x1 <- MMP$SVM[MMP$target==0]
y1 <- plattScalingResult$pred[MMP$target==0]
plot(x1+rnorm(mean=0,sd=0.01,n=length(x1)),
		y1+rnorm(mean=0,sd=0.01,n=length(x1)),
		main="Platt-Scaling Plot",ylim=c(0,1),pch=".",cex=4, 
		xlab="Original values", ylab="Platt-scaled values")
x2 <- MMP$SVM[MMP$target==1]
y2 <- plattScalingResult$pred[MMP$target==1]
points(x2+rnorm(mean=0,sd=0.01,n=length(x2)),
		y2+rnorm(mean=0,sd=0.01,n=length(x2)),
		pch=".",col="red",cex=4)


###################################################
### code chunk number 7: platt.Rnw:161-164
###################################################
newValues <- c(-1.22,0.51,-0.43, 1.1,-1.01)
newValuesPlattScaled <- predictProb(plattScalingResult,newValues)
(cbind(newValues,newValuesPlattScaled))


###################################################
### code chunk number 8: platt.Rnw:234-242
###################################################
p1 <- plattScaling(MMP$NN,MMP$target)$pred
p2 <- plattScaling(MMP$RF,MMP$target)$pred
p3 <- plattScaling(MMP$SVM,MMP$target)$pred
df <- data.frame(p1,p2,p3)
ensemblePrediction <- ensemble(df,
		class1Prob=length(which(MMP$target==1))/nrow(MMP))

(table(ensemblePrediction>0.5, MMP$target))


