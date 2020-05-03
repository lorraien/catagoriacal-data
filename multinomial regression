#####################################################
# Load packages
library(VGAM)
library( nnet )


#################################
# Multinomial

travel = read.csv("/Users/jing/Desktop/210B/R code//travel.csv")
travel$Travel = rep( 1:4, 210 )
travel = travel[ travel$Mode==1, ]
travel = travel[ ,c("Travel", "Hinc", "Psize" ) ]
##
#####	
##
travel[1:10,]
dim( travel )

library( nnet )
mfit = multinom( Travel ~ Hinc + Psize, data=travel )
mfit
summary( mfit )



summ.mfit = function( model ){
	s = summary( model )
	for( i in 1:length(model$coef) ){
		cat( "\nLevel ", model$lev[i+1],  "vs. Level ", model$lev[1], "\n" )
		coef = s$coefficients[i,]
		rrr = exp( coef )
		se = s$standard.errors[i,]
		zStat = coef / se
		pVal = 2*pnorm( abs(zStat), lower.tail=FALSE )
		ci95.lo = exp( coef - qnorm(.975)*se )
		ci95.hi = exp( coef + qnorm(.975)*se )
		rslt = cbind( rrr, se, zStat, pVal, ci95.lo, ci95.hi )
		print( round( rslt, 3 ) )
	}
}

summ.mfit(mfit)

newdata = data.frame(Hinc=45, Psize=3)

# Keeping travel type as names
travel$Travel2[travel$Travel==1] = "Plane"
travel$Travel2[travel$Travel==2] = "Train"
travel$Travel2[travel$Travel==3] = "Bus"
travel$Travel2[travel$Travel==4] = "Car"

mfit2 = multinom( Travel2 ~ Hinc + Psize, data=travel )
mfit2
summary( mfit2 )


# Will use package called VGAM
library(VGAM)

gator = read.table("/Users/jing/Desktop/210B/R code/gator.txt",header=T)
gator$Size = factor(gator$Size,levels=levels(gator$Size)[2:1])

fit = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~1,data=gator,family=multinomial) # Intercept only
deviance(fit); df.residual(fit); pchisq(deviance(fit), df.residual(fit),lower.tail=F) # GOF

fit0 = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Gender,data=gator,family=multinomial) # Gender
deviance(fit0); df.residual(fit0); pchisq(deviance(fit0), df.residual(fit0),lower.tail=F) # GOF

fit1 = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Size,data=gator,family=multinomial) # Size
deviance(fit1); df.residual(fit1); pchisq(deviance(fit1), df.residual(fit1),lower.tail=F)  # GOF

fit2 = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Lake,data=gator,family=multinomial) # Lake
deviance(fit2); df.residual(fit2); pchisq(deviance(fit2), df.residual(fit2),lower.tail=F)  #GOF

fit3 = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Lake+Size,data=gator,family=multinomial) # Lake + Size
deviance(fit3); df.residual(fit3); pchisq(deviance(fit3), df.residual(fit3),lower.tail=F) # GOF

(diff = deviance(fit2) - deviance(fit3));(df = df.residual(fit2)-df.residual(fit3));pchisq(diff,df,lower.tail=F) # LR test for testing significance of Size

fit4 = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Lake+Size+Gender,data=gator,family=multinomial) # Lake + Size + Gender
deviance(fit4); df.residual(fit4); pchisq(deviance(fit4), df.residual(fit4),lower.tail=F) # GOF

(diff = deviance(fit3) - deviance(fit4));(df = df.residual(fit3)-df.residual(fit4));pchisq(diff,df,lower.tail=F) # Gender not significant


fit5 = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Lake*Size,data=gator,family=multinomial) # Lake+Size + Lake:Size
deviance(fit5); df.residual(fit5); pchisq(deviance(fit5), df.residual(fit5),lower.tail=F) # GOF

(diff = deviance(fit3) - deviance(fit5));(df = df.residual(fit3)-df.residual(fit5));pchisq(diff,df,lower.tail=F) # Lake by Size i/a not significant

summary(fit3)  # summary of the model Lake+Size
t(coef(fit3,matrix=T)) # Model coefficients

junk = expand.grid(Size=levels(gator$Size),Lake=levels(gator$Lake))
(predictions = cbind(junk,predict(fit3,type="response",newdata=junk))) # prediction values for all levels of Lake and Size


# Cumulative Logits/ Proportional Odds

mucositis = as.data.frame(list(grade0 = c(5,8),grade1=c(0,8),grade2=c(5,3),grade3=c(16,2),grade4=c(7,0)))
mucositis = cbind(group = c("control","treated"),mucositis)

muc.vglm = vglm(cbind(grade0,grade1,grade2,grade3,grade4)~group,data=mucositis,family=cumulative(parallel=T)) # Fit the model

round(predict(muc.vglm,type="response"),3) # predicted response probabilties from the model
round(mucositis[,2:6]/apply(mucositis[,2:6],1,sum),3) # observed response probabilities
summary(muc.vglm) # model summary

pchisq(deviance(muc.vglm),df.residual(muc.vglm),lower.tail=F) # not a good fit

muc.vglm2= vglm(cbind(grade0,grade1,grade2,grade3,grade4)+1~group,data=mucositis,family=cumulative(parallel=F)) # fit a non-proportional odds model (added 1 to each cell to avoid convergence problems)

summary(muc.vglm2)  # appears first group is problematic
round(predict(muc.vglm2,type="response"),3) # saturated model is predicting Y=0 and 1 better
round(mucositis[,2:6]/apply(mucositis[,2:6],1,sum),3) # observed response probabilities


# Obtaining confidence intervals
x = c(0,0,0,0,-1,1,0,0)
x%*%coef(muc.vglm2) + c(-1,1)*1.96*sqrt(x%*%vcov(muc.vglm2)%*%x) # 95% CI for \beta_2 - \beta_1 (doesn't cover 0, different)

x = c(0,0,0,0,0,-1,1,0)
x%*%coef(muc.vglm2) + c(-1,1)*1.96*sqrt(x%*%vcov(muc.vglm2)%*%x) # 95% CI for \beta_3 - \beta_2 (covers 0, not different)

x = c(0,0,0,0,0,0,-1,1)
x%*%coef(muc.vglm2) + c(-1,1)*1.96*sqrt(x%*%vcov(muc.vglm2)%*%x) # 95% CI for \beta_4 - \beta_3 (covers 0, not different)

muc.vglm = vglm(cbind(grade0+grade1,grade2,grade3,grade4)~group,data=mucositis,family=cumulative(parallel=T)) # combine grade 0 and grade 1
summary(muc.vglm)
pchisq(deviance(muc.vglm),df.residual(muc.vglm),lower.tail=F) # a good fit
round(predict(muc.vglm,type="response"),3) # predicted response probabilties from the model

muc.vglm2 = vglm(cbind(grade0+grade1,grade2,grade3,grade4)~1,data=mucositis,family=cumulative(parallel=T)) # remove beta

(diff=deviance(muc.vglm2) - deviance(muc.vglm) ) # likelihood ratio statistic of H_0 : beta = 0
(df = df.residual(muc.vglm2) - df.residual(muc.vglm))
pchisq(diff,df,lower.tail=F) # p-value

x = c(-1,1,0,0)
x%*%coef(muc.vglm) + c(-1,1)*1.96*sqrt(x%*%vcov(muc.vglm)%*%x) # 95% CI for \alpha_2 - \alpha_1 (doesn't cover 0, different)

x = c(0,-1,1,0)
x%*%coef(muc.vglm) + c(-1,1)*1.96*sqrt(x%*%vcov(muc.vglm)%*%x) # 95% CI for \alpha_3 - \alpha_2 (doesn't cover 0, different)
# Neither of these cover 0 and so the intercepts are significantly different from one another


# Here is an example of fitting the gator data with the multinom function in package nnet
library(nnet)

fit.multinom = multinom(cbind(Fish,Invertebrate,Reptile,Bird,Other)~Size+Lake,data=gator)

summary(fit.multinom)

library(VGAM)

# Adjacent-Categories Logits Model for nominal response

gator = read.table("gator.txt",header=T)
gator$Size = factor(gator$Size,levels=levels(gator$Size)[2:1])

fit.bcl = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Lake+Size,data=gator,family=multinomial) # Lake + Size

fit.acl = vglm(cbind(Invertebrate,Reptile,Bird,Other,Fish)~Lake+Size,data=gator,family=acat(rev=T)) # consult help(acat) 

deviance(fit.bcl)
deviance(fit.acl) # model fits are equivalent
junk = expand.grid(Size=levels(gator$Size),Lake=levels(gator$Lake))
(pred.bcl = cbind(junk,predict(fit.bcl,type="response",newdata=junk))) # pred. same  
(pred.acl = cbind(junk,predict(fit.acl,type="response",newdata=junk))) # for both models

t(coef(fit.bcl,matrix=T)) 
t(coef(fit.acl,matrix=T)) # coefficients are different, but related:

rev(cumsum(rev(coef(fit.acl,matrix=T)["(Intercept)",])))  # These are the alpha_j for the baseline-category logit model, derived from the adjacent-categories logit model

rev(cumsum(rev(coef(fit.acl,matrix=T)["Lakehancock",]))) # effects for Lake Hancock for bcl model, derived from acl model

rev(cumsum(rev(coef(fit.acl,matrix=T)["Lakeoklawaha",]))) # effects for Lake Oklawaha for bcl model, derived from acl model

rev(cumsum(rev(coef(fit.acl,matrix=T)["Laketrafford",]))) # effect for Lake Trafford for bcl model, derived from acl model

rev(cumsum(rev(coef(fit.acl,matrix=T)["Size<2.3",]))) # effects for small alligators for bcl model, derived from acl model


# Adjacent-Categories Logits Model for ordinal response

mucositis = as.data.frame(list(grade0 = c(5,8),grade1=c(0,8),grade2=c(5,3),grade3=c(16,2),grade4=c(7,0)))
mucositis = cbind(group = c("control","treated"),mucositis)

muc.cl = vglm(cbind(grade0+grade1,grade2,grade3,grade4)~group,data=mucositis,family=cumulative(parallel=T)) # combine grade 0 and grade 1
summary(muc.cl)

muc.acl = vglm(cbind(grade0+grade1,grade2,grade3,grade4)~group,data=mucositis,family=acat(parallel=T,rev=T)) # combine grade 0 and grade 1
summary(muc.acl)

deviance(muc.cl)
deviance(muc.acl) # model fits are similar but not equal

predict(muc.cl,type="response")
predict(muc.acl,type="response") # predicted values are close but not equal


# Job Satisfaction Example

# create the data frame
gender = c("female","male")
income = c("<5k","5k-15k","15k-25k",">25k")
vd = c(1,2,0,0,1,0,0,0)
ls = c(3,3,1,2,1,3,0,1)
ms = c(11,17,8,4,2,5,7,9)
vs = c(2,3,5,2,1,1,3,6)
jobsat = expand.grid(income=income,gender=gender)
jobsat = cbind(jobsat,vd,ls,ms,vs) 

# fit the Y nominal independence model (this is the same regardless if X in ordinal or nominal)
fit.ynom.ind = vglm(cbind(vd,ls,ms,vs)~gender,data=jobsat,family=multinomial)

# fit the Y ordinal independence model (this is the same regardless if X in ordinal or nominal)
fit.yord.ind = vglm(cbind(vd,ls,ms,vs)~gender,data=jobsat,family=cumulative(parallel=T))

# fit the Y ordinal, X nominal alternative model
fit.yord.xnom = vglm(cbind(vd,ls,ms,vs)~gender + income,data=jobsat,family=cumulative(parallel=T))

# fit the Y nominal, X nominal alternative model
fit.ynom.xnom = vglm(cbind(vd,ls,ms,vs)~gender+income,data=jobsat,family=multinomial)

# create a new income variable with scores equal to the midpoints of the categories (Income will now be ordinal with a linear effect)
INCOME = c(3,10,20,35,3,10,20,35)

# fit the Y ordinal, X ordinal alternative model
fit.yord.xord = vglm(cbind(vd,ls,ms,vs)~gender + INCOME,data=jobsat,family=cumulative(parallel=T))

# fit the Y nominal, X ordinal alternative model
fit.ynom.xord = vglm(cbind(vd,ls,ms,vs)~gender+INCOME,data=jobsat,family=multinomial)


# Test for conditional independence in the Y nom, X nom model

(diff = deviance(fit.ynom.ind) - deviance(fit.ynom.xnom))
(df = df.residual(fit.ynom.ind) - df.residual(fit.ynom.xnom))
pchisq(diff,df,lower.tail=F) # Not significant, cannot reject null hypothesis of conditional independence.


# Test for conditional independence in the Y nom, X ord model

(diff = deviance(fit.ynom.ind) - deviance(fit.ynom.xord))
(df = df.residual(fit.ynom.ind) - df.residual(fit.ynom.xord))
pchisq(diff,df,lower.tail=F) # Near significance (note it helps to treat X as ordinal)

# Test for conditional independence in the Y ord, X nom model

(diff = deviance(fit.yord.ind) - deviance(fit.yord.xnom))
(df = df.residual(fit.yord.ind) - df.residual(fit.yord.xnom))
pchisq(diff,df,lower.tail=F) # significant, reject the null of conditional independence (note it helps to treat Y as ordinal)


# Test for conditional independence in the Y ord, X ord model (most powerful model/test)

(diff = deviance(fit.yord.ind) - deviance(fit.yord.xord))
(df = df.residual(fit.yord.ind) - df.residual(fit.yord.xord))
pchisq(diff,df,lower.tail=F) # significant, reject the null of conditional independence 
