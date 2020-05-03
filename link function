# Power as a function of the chi^2 noncentrality parameter

df =1
pchisq(qchisq(.95,df),df,ncp=seq(0,20,.1),lower.tail=F) -> x
df =7
pchisq(qchisq(.95,df),df,ncp=seq(0,20,.1),lower.tail=F) -> y
df =4
pchisq(qchisq(.95,df),df,ncp=seq(0,20,.1),lower.tail=F) -> z

matplot(seq(0,20,.1),cbind(x,z,y),lty=c(1,2,4),type="l",col=1,xlab=expression(paste(chi^2," Noncentrality Parameter, ",lambda,sep="")),ylab="Power",main=expression(paste(Pr({X^2}[nu][lambda] > {chi[nu]^2}(alpha))," as a function of ",lambda," and df",sep="")))
legend(10,.3,c("df = 1","df = 4","df = 7"),lty=c(1,2,4),col=1)


# Using clinical data to fit linear logit and saturated model
infection = read.table("/Users/jing/Desktop/210B/R code/clinical.txt",header=T)

# Saturated model (nominal covariate)
fit = glm(cbind(High,Low)~Clinical,family=binomial(),data=infection)
summary(fit)

# Linear logit model (ordinal covariate)
fit2 = glm(cbind(High,Low)~as.numeric(Clinical),family=binomial(),data=infection)
summary(fit2)


#  Beetle Mortality


ldose = c(1.691,1.724,1.755,1.784,1.811,1.837,1.861,1.884)
nbeetles = c(59,60,62,56,63,59,62,60)
nkilled = c(6,13,18,28,52,53,61,60)

beetle.logit.fit = glm(nkilled/nbeetles ~ ldose,weights=nbeetles,family=binomial())
summary(beetle.logit.fit)

beetle.probit.fit = glm(nkilled/nbeetles ~ ldose,weights=nbeetles,family=binomial(link="probit"))
summary(beetle.probit.fit)

beetle.cloglog.fit = glm(nkilled/nbeetles ~ ldose,weights=nbeetles,family=binomial(link="cloglog"))
summary(beetle.cloglog.fit)


pb.logit = predict(beetle.logit.fit,type="response",newdata=data.frame(ldose=seq(1.6,1.9,.01)))
pb.probit = predict(beetle.probit.fit,type="response",newdata=data.frame(ldose=seq(1.6,1.9,.01)))
pb.cloglog = predict(beetle.cloglog.fit,type="response",newdata=data.frame(ldose=seq(1.6,1.9,.01)))
 

plot(ldose,nkilled/nbeetles,xlab="Log Dose",ylab="Proportion Killed",pch=16)
lines(seq(1.6,1.9,.01),pb.cloglog)
lines(seq(1.6,1.9,.01),pb.probit,lty=2)
lines(seq(1.6,1.9,.01),pb.logit,lty=6)
legend(1.8,.3,c("cloglog","probit","logit"),lty=c(1,2,6),col=1)

