# TWO-WAY TABLES

deg = factor(c("lhs","hs.jc","bs.ms.phd"),levels=c("lhs","hs.jc","bs.ms.phd"))
bel = factor(c("fund","mod","lib"),levels=c("fund","mod","lib"))
relig.belief = cbind(expand.grid(belief=bel,degree=deg),counts=c(178,138,108,570,648,442,138,252,252))
ftable(xtabs(counts~degree+belief,data=relig.belief))

fit.ind = glm(counts~degree+belief,data=relig.belief,family=poisson) # independence model
summary(fit.ind)

deviance(fit.ind) # G^2 statistic
df.residual(fit.ind) # residual degrees of freedom
pchisq(deviance(fit.ind),df.residual(fit.ind),lower.tail=F) # reject null of XY independence ？test model is good fit？

(pearson.ind = sum(residuals(fit.ind,type="pearson")^2)) # pearson X^2
pchisq(pearson.ind,df.residual(fit.ind),lower.tail=F) # reject XY independence

fit.sat =  glm(counts~degree*belief,data=relig.belief,family=poisson)
summary(fit.sat)

# Sample odds ratio n_11 n_22 / n_21 n_12

178*648/570/138  # sample based

exp(0+coef(fit.sat)["degreehs.jc:beliefmod"] + 0 + 0) # model based estimated odds ratio

# Sample odds ratio  n_22 n_33 / n_32 n_23

648*252/252/442 # sample base

exp(coef(fit.sat)["degreehs.jc:beliefmod"] 
   + coef(fit.sat)["degreebs.ms.phd:belieflib"] 
   - coef(fit.sat)["degreebs.ms.phd:beliefmod"] 
   - coef(fit.sat)["degreehs.jc:belieflib"])  # model based estimated odds ratio


# THREE_WAY TABLES

#Example: Alcohol, Cigarette and Marijuana Use among HS Seniors.

tmp = factor(c("yes","no"),levels=c("no","yes"))
counts = c(911,538,44,456,3,43,2,279)
drug.use = cbind(expand.grid(marijuana=tmp,cigarette=tmp,alcohol=tmp),counts)

ftable(xtabs(counts~alcohol+cigarette+marijuana,data=,drug.use))

fit.sat = glm(counts~.^3,data=drug.use,family=poisson()) #saturated model
fitted.sat =fitted(fit.sat) # fitted values equal the observed counts

fit.homog.assoc = update(fit.sat,~.-alcohol:marijuana:cigarette)# homogeneous association model
fitted.ha = round(fitted(fit.homog.assoc),1) # fitted values close to observed

fit.am.cm = update(fit.homog.assoc,~.-alcohol:cigarette) # AC conditional independence
fitted.am.cm = round(fitted(fit.am.cm),1)

fit.ac.am = update(fit.homog.assoc,~.-cigarette:marijuana) # CM conditional independence
fitted.ac.am = round(fitted(fit.ac.am),1)

fit.ac.cm = update(fit.homog.assoc,~.-alcohol:marijuana) # AM conditional independence
fitted.ac.cm = round(fitted(fit.ac.cm),1)

fit.a.cm = update(fit.homog.assoc,~.-alcohol:marijuana-alcohol:cigarette) # A jointly independent of C and M
fitted.a.cm = round(fitted(fit.a.cm),1)

fit.c.am = update(fit.homog.assoc,~.-cigarette:marijuana-cigarette:alcohol) # C jointly independent of A and M
fitted.c.am = round(fitted(fit.c.am),1)

fit.m.ac = update(fit.homog.assoc,~.-alcohol:marijuana-cigarette:marijuana) # M jointly independent of A and C
fitted.m.ac = round(fitted(fit.m.ac),1)

fit.a.c.m = glm(counts~alcohol+cigarette+marijuana,data=drug.use,family=poisson) # Joint independence
fitted.a.c.m = round(fitted(fit.a.c.m),1)

fitted.table = cbind(fitted.sat,fitted.ha,fitted.am.cm,fitted.m.ac,fitted.a.c.m)

fitted.table = cbind(drug.use[,3:1],fitted.table)
fitted.table
