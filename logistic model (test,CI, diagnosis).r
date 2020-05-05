# read in the crab data set

hcrabs = read.table("crab.txt",header=T)

crabs.width.factor = cut(hcrabs$Width,breaks=c(0,seq(from=23.25,to=29.25,by=1),Inf)) #make a factor variable out of the continuous variable Width

x = aggregate(hcrabs$Width,by=list(Width=crabs.width.factor),mean)$x  #compute means for Width at each level of the factor crabs.width.factor

bin.response = ifelse(hcrabs$Sat>0,1,0)  # if >=1 satellite set to 1 else to 0

 crabs.width.factor = cut(hcrabs$Width,breaks=c(0,seq(from=23.25,to=29.25,by=1),Inf)) #make a factor variable out of the continuous variable Width
postscript("crabs1.ps",horizontal=F,bg="white",colormodel="rgb")
 plot(hcrabs$Width,bin.response,pch=1,col=4,xlab="Width (cm)",ylab="Prob. of Presence of satellites")  # set up the plot region and add labels
dev.off()


postscript("crabs2.ps",horizontal=F,bg="white",colormodel="rgb")
  plot(hcrabs$Width,bin.response,pch=1,col=4,xlab="Width (cm)",ylab="Prob. of Presence of satellites")  # set up the plot region and add   labels
  prop = aggregate(bin.response,by=list(W=crabs.width.factor),mean)$x # get proportion at each width factor level
  points(x,prop,pch=16) # plot the sample proportions
dev.off()

crabs.fit = glm(bin.response~Width,data=hcrabs,family=binomial()) # fit a logistic regression model with Width as predictor

round(crabs.fit$coeff,3) # coefficients rounded to 3 significant digits
crabs.fit$coeff[1]/crabs.fit$coeff[2] # LD_50
exp(crabs.fit$coeff[2])  # multiplicative increase in odds for 1 cm increase in width
100*(exp(crabs.fit$coeff[2])-1) # percent increase in odds for 1 cm increase in width

new.widths = seq(20,34,by=.1)
y = predict(crabs.fit,type="response",newdata=data.frame(Width=new.widths)) # predict new values for response

postscript("crabs3.ps",horizontal=F,bg="white",colormodel="rgb")
 plot(hcrabs$Width,bin.response,pch=1,col=4,xlab="Width (cm)",ylab="Prob. of Presence of satellites")  # set up the plot region and add labels
 prop = aggregate(bin.response,by=list(W=crabs.width.factor),mean)$x # get proportion at each width factor level
 points(x,prop,pch=16) # plot the sample proportions
 lines(new.widths,y,col=4) # plot the estimated values of pi
dev.off()


#####################################
# Cigarette retrospective studies example

lc = read.table("lungcancer.txt",header=T)

postscript("lc1.ps",horizontal=F,bg="white",colormodel="rgb")
plot(lc$score,lc$cases/(lc$cases+lc$controls),xlab="Smoking Level",ylab="Prob. of Lung Cancer",pch=16)
dev.off()

lc.fit = glm(cbind(cases,controls) ~ score,family=binomial(),data=lc)
summary(lc.fit)
pchisq(summary(lc.fit)$deviance,summary(lc.fit)$df.residual,lower=F) # p-value of goodness-of-fit test based on LR statistic

postscript("lc2.ps",horizontal=F,bg="white",colormodel="rgb")
plot(lc$score,lc$cases/(lc$cases+lc$controls),xlab="Smoking Level",ylab="Prob. of Lung Cancer",pch=16)
lines(seq(0,5,.1),predict(lc.fit,type="response",newdata=data.frame(score=seq(0,5,.1))))
dev.off()

# Removing the base group of no cigarettes
lc.fit2 = glm(cbind(cases,controls) ~ score,family=binomial(),data=lc,subset= score>0)

summary(lc.fit2)

postscript("lc3.ps",horizontal=F,bg="white",colormodel="rgb")
lc.subset = lc$score > 0
plot(lc$score[lc.subset],lc$cases[lc.subset]/(lc$cases[lc.subset]+lc$controls[lc.subset]),xlab="Smoking Level",ylab="Prob. of Lung Cancer",pch=16)
lines(seq(0,5,.1),predict(lc.fit2,type="response",newdata=data.frame(score=seq(0,5,.1))))
dev.off()

pchisq(summary(lc.fit2)$deviance,summary(lc.fit2)$df.residual,lower=F) # p-value of goodness-of-fit test based on LR statistic

pchisq(summary(lc.fit2)$null.deviance-summary(lc.fit2)$deviance,summary(lc.fit2)$df.null-summary(lc.fit2)$df.residual,lower=F) # LR test of H_0 : beta = 0

# Read in the framingham data
framingham = read.table("Framingham.txt", header=T)

# Convert sex to 0 or 1 and create a female variable
framingham$sex = framingham$sex - 1
framingham$female =framingham$sex

# Create a new variable, which categorizes the blood pressure of a subject
framingham$sbpgrp = cut( framingham$sbp,breaks=c(min(framingham$sbp),126,146,166, max(framingham$sbp)),include.lowest=TRUE )

#################################### diagnosis
#####	H-L goodness of fit test
##
binary.gof <- function( fit, ngrp=10, print.table=TRUE ){
	y <- fit$y
	phat <- fitted( fit )
	fittedgrps <- cut( phat, quantile( phat, seq(0,1,by=1/ngrp) ), include.lowest=TRUE )
	n <- aggregate( y, list( fittedgrps ), FUN=length )[,2]
	Obs <- aggregate( y, list( fittedgrps ), FUN=sum )[,2]
	Exp <- aggregate( phat, list( fittedgrps ), FUN=sum )[,2]
	if( print.table==TRUE ){
		cat( "\nFitted Probability Table:\n\n" )
		rslt <- as.data.frame( cbind( 1:ngrp, n, Obs, Exp ) )
		names( rslt )[1] <- "group"
		print( rslt )
	}
	chisqstat <- sum( (Obs - Exp)^2 / ( Exp*(1-Exp/n) ) )
	df <- ngrp-2
	pVal <- pchisq( chisqstat, df, lower.tail=FALSE )
	cat( "\n Hosmer-Lemeshow GOF Test:\n\n" )
	cbind( chisqstat, df, pVal )
}

# Fit full model
full.chd = glm(chdfate~sbp+age, family=binomial(link="logit"),data= framingham)

binary.gof(full.chd)

summary(full.chd)

# ROC Curve
# Will need to install pROC package
# install.packages("pROC)
library(pROC)
roc.curve = roc(framingham$chdfate~fitted(full.chd))
plot(roc.curve)
roc.curve

####################################################
# grouped vs un-grouped data example

snoring = matrix(c(24,35,21,30,1355,603,192,224), nc=2)

snoring.df = data.frame(snoring=gl(4, 1, labels=c("Never", "Occasional", "Nearly every night", "Every night")), disease=gl(2, 4, labels=c("Yes", "No")), counts=as.vector(snoring))
snoring.df = snoring.df[rep(seq_len(nrow(snoring.df)), snoring.df$counts), 1:2]

# snoring is the grouped data and snoring.df is the ungrouped

levels(snoring.df$snoring) = c(0, 2, 4, 5)
y = abs(as.numeric(snoring.df$disease)-2)
x = as.numeric(as.character(snoring.df$snoring))

# This is the ungrouped model
fit.glm1 = glm(y ~ x, family=binomial)

# This is the grouped model
fit.glm2 = glm(snoring ~ c(0, 2, 4, 5), family=binomial)

# Both models are identical in output
summary(fit.glm1)
summary(fit.glm2)

####################################################
# Cochran-Armitage test for trend (similar to linear trend test from before)

malformed = c(48,38,5,1,1)
not_mal = c(17066,14464,788,126,37)
scores = c(0,0.5,1.5,4,7)
x = c(rep(scores,malformed), rep(scores, not_mal)) 
y = c(rep(1,sum(malformed)), rep(0, sum(not_mal)))

z2 = sum(malformed + not_mal)*cor(x,y)^2
z2
pchisq(z2,1,lower.tail=F)
  
# compare to linear logit model

malform.fit = glm(cbind(malformed,not_mal)~scores,family=binomial)
summary(malform.fit)

# AZT and AIDS example

aids = read.table("/Users/jing/Desktop/210B/R code/aids.txt",header=T)
ftable(xtabs(cbind(yes,no)~race + AZT,data=aids))

aids.fit = glm(cbind(yes,no)~AZT + race,family=binomial(),data=aids)
summary(aids.fit)

####################################################
# tests for conditional independence

 
aids.fit2 = update(aids.fit,formula = ~.-AZT) # fit reduce model without AZT 
anova(aids.fit2,aids.fit,test="Chisq") #  compute difference in deviances (LR ratio statistic)


aids.fit3 = update(aids.fit,formula = ~.-race)
anova(aids.fit3,aids.fit,test="Chisq") # LR test of H_0 : \beta_2^Z = 0

####################################################
# GOF

aids.fit$deviance  # get the model fit deviance
sum(residuals(aids.fit,type="deviance")^2) # an alternative way to get the deviance

aids.fit$df.residual # and its d.o.f. (n - p)
pchisq(aids.fit$deviance,aids.fit$df.residual,lower.tail=F)

# alternatively

(X2 = sum(residuals(aids.fit,type="pearson")^2)) # yet another way using X^2 instead of G^2
1-pchisq(X2,aids.fit$df.residual)

anova(aids.fit,test="Chisq") # Adds terms sequentially, one a a time, based on model formula


# stepwise selection

crab =read.table("/Users/jing/Desktop/210B/R code/crab.txt",header=T)
crab$Wt = crab$Wt/1000
crab$C = factor(crab$C) #  nominal to start
crab$S = factor(crab$S) # nominal to start

bin.response = ifelse(crab$Sat>0,1,0) 
crab = cbind(crab,bin=bin.response)

# first forward selection
crab.fit1 = glm(bin~1,family=binomial(),data=crab)
fstep.crab.fit = step(crab.fit1,scope=list(lower=formula(crab.fit1),upper= bin ~ C*S*Width),trace=F,direction="forward")
fstep.crab.fit$anova
summary(fstep.crab.fit)

# next backward selection

crab.fit2 = glm(bin~C*S*Width,family=binomial(),data=crab)
bstep.crab.fit = step(crab.fit2,direction="backward",trace=F)
bstep.crab.fit$anova
summary(bstep.crab.fit)

# finally both directions, starting from full model

step.crab.fit = step(crab.fit2,direction="both",trace=F)
step.crab.fit$anova
summary(step.crab.fit)


#####################################################
# Functions to be used to create confidence intervals in GLM (and also conduct likelihood ratio test and robust variance estimation)


# Ifelse function

ifelse1 =function(test, x, y){ if (test) x else y}

#####################################################
#  Function to exponentiate coefficients and produces CIs for GLMs

glmCI <- function( model, transform=TRUE, robust=FALSE ){
	link <- model$family$link
	coef <- summary( model )$coef[,1]
	se <- ifelse1( robust, robust.se.glm(model)[,2], summary( model )$coef[,2] )
	zvalue <- coef / se
	pvalue <- 2*(1-pnorm(abs(zvalue)))

	if( transform & is.element(link, c("logit","log")) ){
		ci95.lo <- exp( coef - qnorm(.975) * se )
		ci95.hi <- exp( coef + qnorm(.975) * se )
		est <- exp( coef )
	}
	else{
		ci95.lo <- coef - qnorm(.975) * se
		ci95.hi <- coef + qnorm(.975) * se
		est <- coef
	}
	rslt <- round( cbind( est, ci95.lo, ci95.hi, zvalue, pvalue ), 4 )
	colnames( rslt ) <- ifelse1( 	robust, 	
					c("Est", "robust ci95.lo", "robust ci95.hi", "robust z value", "robust Pr(>|z|)"),
					c("Est", "ci95.lo", "ci95.hi", "z value", "Pr(>|z|)") )			
	colnames( rslt )[1] <- ifelse( transform & is.element(link, c("logit","log")), "exp( Est )", "Est" )
	rslt
	}
	
#####################################################
#	Function to estimate linear contrasts of coefficients from a GLM fit

linContr.glm <- function( contr.names, contr.coef=rep(1,length(contr.names)), model, transform=TRUE ){
	beta.hat <- model$coef 
	cov.beta <- vcov( model )

	contr.index <- match( contr.names, dimnames( cov.beta )[[1]] )	
	beta.hat <- beta.hat[ contr.index ]
	cov.beta <- cov.beta[ contr.index,contr.index ]
	est <- contr.coef %*% beta.hat
	se.est <- sqrt( contr.coef %*% cov.beta %*% contr.coef )
	zStat <- est / se.est
	pVal <- 2*pnorm( abs(zStat), lower.tail=FALSE )
	ci95.lo <- est - qnorm(.975)*se.est
	ci95.hi <- est + qnorm(.975)*se.est
	
	link <- model$family$link
	if( transform & is.element(link, c("logit","log")) ){
		ci95.lo <- exp( ci95.lo )
		ci95.hi <- exp( ci95.hi )
		est <- exp( est )
		cat( "\nTest of H_0: exp( " )
		for( i in 1:(length( contr.names )-1) ){
			cat( contr.coef[i], "*", contr.names[i], " + ", sep="" )
			}
		cat( contr.coef[i+1], "*", contr.names[i+1], " ) = 1 :\n\n", sep="" )		
		}
	else{
		cat( "\nTest of H_0: " )
		for( i in 1:(length( contr.names )-1) ){
			cat( contr.coef[i], "*", contr.names[i], " + ", sep="" )
			}
		cat( contr.coef[i+1], "*", contr.names[i+1], " = 0 :\n\n", sep="" )
		}
	rslt <- data.frame( est, se.est, zStat, pVal, ci95.lo, ci95.hi )
	colnames( rslt )[1] <- ifelse( transform && is.element(link, c("logit","log")), "exp( Est )", "Est" )
	round( rslt, 8 )
}

#####################################################
# Function to compute deviance (LR) test p-Value

lrtest <- function( fit1, fit2 ){
	cat( "\nAssumption: Model 1 nested within Model 2\n\n" )
	rslt <- anova( fit1, fit2 )
	rslt <- cbind( rslt, c("", round( pchisq( rslt[2,4], rslt[2,3], lower.tail=FALSE ), 4 ) ) )
	rslt[,2] <- round( rslt[,2], 3 )
	rslt[,4] <- round( rslt[,4], 3 )
	rslt[1,3:4] <- c( "", "" )
	names( rslt )[5] <- "pValue"
	rslt
}

#####	H-L goodness of fit test
##
binary.gof <- function( fit, ngrp=10, print.table=TRUE ){
	y <- fit$y
	phat <- fitted( fit )
	fittedgrps <- cut( phat, quantile( phat, seq(0,1,by=1/ngrp) ), include.lowest=TRUE )
	n <- aggregate( y, list( fittedgrps ), FUN=length )[,2]
	Obs <- aggregate( y, list( fittedgrps ), FUN=sum )[,2]
	Exp <- aggregate( phat, list( fittedgrps ), FUN=sum )[,2]
	if( print.table==TRUE ){
		cat( "\nFitted Probability Table:\n\n" )
		rslt <- as.data.frame( cbind( 1:ngrp, n, Obs, Exp ) )
		names( rslt )[1] <- "group"
		print( rslt )
	}
	chisqstat <- sum( (Obs - Exp)^2 / ( Exp*(1-Exp/n) ) )
	df <- ngrp-2
	pVal <- pchisq( chisqstat, df, lower.tail=FALSE )
	cat( "\n Hosmer-Lemeshow GOF Test:\n\n" )
	cbind( chisqstat, df, pVal )
}

#####	Function to compute robust se for glms
##
robust.se.glm<-function(glm.obj){
	## 	Compute robust (sandwich) variance estimate
	if (is.matrix(glm.obj$x)) 
		xmat<-glm.obj$x
	else {
		mf<-model.frame(glm.obj)
		xmat<-model.matrix(terms(glm.obj),mf)		
	}
	umat <- residuals(glm.obj,"working")*glm.obj$weights*xmat
	modelv<-summary(glm.obj)$cov.unscaled
	robust.cov <- modelv%*%(t(umat)%*%umat)%*%modelv
	
	##	Format the model output with p-values and CIs
	s <- summary( glm.obj) 
	robust.se <- sqrt( diag( robust.cov )) 
	z <- glm.obj$coefficients/robust.se
	p <- 2*pnorm( -abs( z ) ) 
	ci95.lo <- glm.obj$coefficients - qnorm( .975 ) * robust.se
	ci95.hi <- glm.obj$coefficients + qnorm( .975 ) * robust.se
	rslt <- cbind( glm.obj$coefficients, robust.se, ci95.lo, ci95.hi, z, p ) 
	dimnames(rslt)[[2]] <- c( dimnames( s$coefficients )[[2]][1], "Robust SE", "ci95.lo", "ci95.hi", dimnames( s$coefficients )[[2]][3:4] ) 
	rslt 
	}

##############################
# Examples from lecture notes

# Read in data
framingham = read.table("/Users/jing/Desktop/210B/R code/framingham.txt")

# Convert sex to 0 or 1 and create a female variable
framingham$sex = framingham$sex - 1
framingham$female =framingham$sex

# Create a new variable, which categorizes the blood pressure of a subject
framingham$sbpgrp = cut( framingham$sbp,breaks=c(min(framingham$sbp),126,146,166, max(framingham$sbp)),include.lowest=TRUE )

# Fit model with interaction
m1 = glm(chdfate~sbp+age+sbp*age, family=binomial(link="logit"),data= framingham)
summary(m1)

# Get exponentiated confidence intervals for coefficients
glmCI(m1)

# Obtain exponentiated cofidence intervals for linear combinations of parameters
linContr.glm(c("sbp","sbp:age"), c(1,20), model=m1)

linContr.glm(c("sbp","sbp:age"), c(15,15*20), model=m1)


# Testing several coefficients at once
reduced = glm(chdfate~sbp, family=binomial(link="logit"),data= framingham)
full = glm(chdfate~sbp+age+sbp*age, family=binomial(link="logit"),data= framingham)

lrtest(reduced,full)
