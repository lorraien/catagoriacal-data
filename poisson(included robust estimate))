# Run these functions into R

#####################################################
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
	
	
# Poisson Regression section

# Read in crab data and fit some models
crab = read.csv("/Users/jing/Desktop/210B/R code/crab.csv", header=TRUE)

crab.model = glm(Sa~W, family=poisson, data=crab)

crab.model2 = glm(Sa~W+Wt+W*Wt, family=poisson, data=crab)


# Read in ERtemp data and fit models
ERtemp = read.csv("/Users/jing/Desktop/210B/R code/ERtemp.csv", header=TRUE)

ER.model = glm(Admissions~Temperature, family=poisson, data=ERtemp)

linContr.glm(c("(Intercept)", "Temperature"), c(1, 85), model = ER.model)

# Read in std data and fit some models
std = read.csv("/Users/jing/Desktop/210B/R code/stdgrp.csv")

std.model = glm(n.reinfect~condom.always+offset(log(yrsfu)), family=poisson, data=std)

std.model2 = glm(n.reinfect~white+condom.always+edugrp+offset(log(yrsfu)), family=poisson, data=std)

std.model3 = glm(n.reinfect~white+condom.always+relevel(edugrp, ref="[6,11.9]")+offset(log(yrsfu)), family=poisson, data=std)

# LR test to see between full model and reduced model
lrtest(std.model,std.model2)

glmCI(std.model2)

# Read in Gail data, transform latitude, and fit model.
gail = read.table("/Users/jing/Desktop/210B/R code/gail.txt", header=TRUE)

gail$latitude.ord = 1
gail$latitude.ord[ gail$latitude=="Middle  " ] = 2
gail$latitude.ord[ gail$latitude=="Southern" ] = 3

model.gail = glm(formula = inccases ~ factor(ageg) + latitude.ord, family = poisson(link=log), data = gail, offset = log(persyrs))

# Fit Gail model with quasipoisson
model.gail2 = glm(formula = inccases ~ factor(ageg) + latitude.ord, family = quasipoisson(link=log), data = gail, offset = log(persyrs))

# Obtain robust standard errors for model.gail
robust.se.glm(model.gail)


#################################

# Credit card example
creditcard = read.table("/Users/jing/Desktop/210B/R code/creditcard.txt")

lcases= log(creditcard[,2])
creditcard =cbind(creditcard,lcases)
colnames(creditcard)=c("income","cases","CrCards","lcases")

card.model=glm(CrCards~income+offset(lcases),family=poisson,data= creditcard)
summary(card.model)
 	
