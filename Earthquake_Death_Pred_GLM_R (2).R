library(MASS)
library(boot)

## Load the Earthquake data (with magnitude imputed)
signif<- read.csv("/Users/jyotsnakumar/Documents/BDAP-04/Project/signif_final.csv")
head(signif)
summary((signif))
str(signif)
signif$REGION_CODE<-as.factor(signif$REGION_CODE)

## Create a dataframe where Intensity is not null
##2321 records have data for intensity
nrow(signif[!is.na(signif$INTENSITY),])
signif_intensity<-signif[!is.na(signif$INTENSITY),]
summary(signif_intensity)
nrow(signif_intensity)

## Create a dataframe from Intensity dataset where Death is not null
## There are 970 records with Death
nrow(signif_intensity[!is.na(signif_intensity$DEATHS),])
signif_death<-signif_intensity[!is.na(signif_intensity$DEATHS),]
nrow(signif_death)
summary(signif_death)
str(signif_death)

##Initial analysis

##Mean - Variance Ratio
d<-c(mean(signif_death$DEATHS),var(signif_death$DEATHS))
c(mean=d[1],var=d[2],ratio=d[2]/d[1])


## Average deaths by Countries
with(signif_death, tapply(DEATHS, COUNTRY, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

##Observation
## We see that Mean varies with Countries and that Conditional Variance
## is greater than Conditional Mean

## Average deaths by Earthquake Intensity
with(signif_death, tapply(DEATHS, INTENSITY, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))

##Observation
## We see that Mean varies with the Earthquake Intensity and that it is 
## an important parameter
##Conditional Variance is greater than Conditional Mean

##The above statistics is suggestive of over-dispersion in the data and
#as to why a negative binomial is more appropriate than Poisson 

## Split the data into training  and validation sets.
set.seed(81)
index <- sample(1:nrow(signif_death),800)
train_death <- signif_death[index,]
validate_death <- signif_death[-index,]

summary(train_death)
str(validate_death)
unique(train_death$COUNTRY)

## Create Negative binomial model with dispersion parameter (theta) =1
m1 <- glm(DEATHS ~  INTENSITY+LATITUDE+LONGITUDE+EQ_MAG_MW+DEATHS_DESCRIPTION+COUNTRY+REGION_CODE ,family=negative.binomial(theta=1,link="log"),data=train_death)

## Create Poisson model
m2 <- glm(DEATHS ~  INTENSITY+LATITUDE+LONGITUDE+EQ_MAG_MW+DEATHS_DESCRIPTION+COUNTRY+REGION_CODE , family=poisson(link="log"),data=train_death)

##Check summary of the models 
summary(m1)

summary(m2)

## Comparison between observed and fitted values
data.frame(train_death$DEATHS,m1$fitted.values,m2$fitted.values)

## Finding Pearson's Chi squared for Poisson and by computing
#scalar parameter  phi
pr<-residuals(m2,"pearson")
phi<-sum(pr^2)/m2$df.residual
phi

##Predict both models on the validation set
p_m1<-predict(m1, validate_death, type="response", se.fit=TRUE)
p_m2<-predict(m2, validate_death, type="response", se.fit=TRUE)

validate_death$death_pred_nb <- p_m1$fit
validate_death$death_pred_poiss <- p_m2$fit

##Com
validate_death[,c("DEATHS","death_pred_nb","death_pred_poiss")]

###------Metrics-------###

### Print AIC, BIC ####

sprintf("Negative Binomial:   AIC - %f , BIC - %f",AIC(m1),BIC(m1))
sprintf("Poisson:   AIC - %f , BIC - %f",AIC(m2),BIC(m2))

## Print Log Likelihood##
sprintf("Negative Binomial: %f , Poisson: %f",logLik(m1),logLik(m2))

## Likelihood  Ratio Test##
pchisq(2*(logLik(m1)-logLik(m2)),df=1,lower.tail=FALSE)

## Significant, hence negative binomial which estimates 
##the dispersion parameter is better than Poisson

pchisq(m1$deviance,m1$df.residual)
pchisq(m2$deviance,m2$df.residual)

## Overdispersion test
sprintf("Negative Binom: %f ,Poisson: %f",m1$deviance/m1$df.residual,
        m2$deviance/m2$df.residual)

## Plot observed values with fitted Poisson and Negative binomial distribution
vh<-hist(validate_death$DEATHS,
        prob=TRUE,col="red",border="white",xlim=c(0,900000),breaks=50,
        freq=F,main="Observed Deaths")
vh

v_m1_h<-hist(validate_death$death_pred_nb,
           prob=TRUE,col="blue",border="white",xlim=c(0,100000),breaks=50,
           freq=F,main="Negative binomial")
v_m1_h

v_m2_h<-hist(validate_death$death_pred_poiss,
           prob=TRUE,col="green",border="white",xlim=c(0,100000),breaks=50,
           freq=F,main="Poisson")
v_m2_h

par(mfrow=c(1,1))

png("/Users/jyotsnakumar/Documents/BDAP-04/Project/Distribution_Comparison.png",width=500,height=400)
plot(vh$mids, vh$counts,col="red",type='l',xlab="Deaths",ylab="Frequency",
     main="Model Distribution Comparison on Prediction") 
lines(v_m2_h$mids,v_m2_h$counts,col='green')      
lines(v_m1_h$mids,v_m1_h$counts,col='blue')
legend("topright", legend=c("Observed", "Poisson","Negative-Binomial"), col=c("red", "green","blue"), ncol=1, lty=1)
dev.off()

##QQplot

png("/Users/jyotsnakumar/Documents/BDAP-04/Project/QQplot.png",width=500,height=400)
qqnorm(m1$residuals)

dev.off()

###### Mean-Variance Plot of Quasi-Poisson vs Negative Binomial #####

g <- cut(m1$fitted.values, breaks=quantile(m1$fitted.values,seq(0,100,2)/100))

m <- tapply(train_death$DEATHS, g, mean)

v <- tapply(train_death$DEATHS, g, var)

png("/Users/jyotsnakumar/Documents/BDAP-04/Project/Mean-Variance-plot.png",width=500,height=400)

plot(m, v, xlab="Mean", ylab="Variance",main="Mean-Variance Relationship")

x <- seq(min(m),max(m),400)

lines(x, phi*x, lty="dashed",col='blue')
lines(x, x*(1+x/1.0008),col='red')
legend("topleft", lty=c("dashed","solid"), col=c('blue','red'),
            legend=c("Q. Poisson","Neg. Binom."), inset=0.05)

dev.off()

##Conclusion -
#Negative binomial provides a better fit to the data than Poisson as the 
##plot is quadratic and approaches data
## Negative Binomial has better AIC and BIC as compared to Poisson.
## The Log Likelihood ratio also suggests that  a negative binomial is a
##better fit to the data as the data is over-dispersed.
## Both the methods do not handle extereme values very well.However,
##negative binomial still does a better job in predicting the count of deaths
##as compared to Poisson


