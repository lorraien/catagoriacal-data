# THREE_WAY TABLES for Loglinear models

# GOF Tests for the Drug use example

tmp = factor(c("yes","no"),levels=c("no","yes"))
counts = c(911,538,44,456,3,43,2,279)
drug.use = cbind(expand.grid(marijuana=tmp,cigarette=tmp,alcohol=tmp),counts)

ftable(xtabs(counts~alcohol+cigarette+marijuana,data=,drug.use))

fit.sat = glm(counts~.^3,data=drug.use,family=poisson()) #saturated model

fit.homog.assoc = update(fit.sat,~.-alcohol:marijuana:cigarette)# homogeneous association model

fit.am.cm = update(fit.homog.assoc,~.-alcohol:cigarette) # AC conditional independence

fit.ac.am = update(fit.homog.assoc,~.-cigarette:marijuana) # CM conditional independence

fit.ac.cm = update(fit.homog.assoc,~.-alcohol:marijuana) # AM conditional independence

fit.a.cm = update(fit.homog.assoc,~.-alcohol:marijuana-alcohol:cigarette) # A jointly independent of C and M

fit.c.am = update(fit.homog.assoc,~.-cigarette:marijuana-cigarette:alcohol) # C jointly independent of A and M

fit.m.ac = update(fit.homog.assoc,~.-alcohol:marijuana-cigarette:marijuana) # M jointly independent of A and C

fit.a.c.m = glm(counts~alcohol+cigarette+marijuana,data=drug.use,family=poisson) # Joint independence


model2 = list(fit.sat,fit.homog.assoc,fit.am.cm,fit.ac.am,fit.ac.cm,fit.a.cm,fit.c.am,fit.m.ac,fit.a.c.m)

G2 = round(unlist(lapply(model2,deviance)),1)
X2 = round(unlist(lapply(lapply(lapply(model2,residuals,type="pearson"),"^",2),sum)),1)
AIC = round(unlist(lapply(model2,AIC)),1)
df = round(unlist(lapply(model2,df.residual)),1)
pvalue = rep(0,9)
pvalue = round(pchisq(G2,df,lower.tail=F),3)

model = c("(ACM)","(AC,AM,CM)","(AM,CM)","(AC,AM)","(AC,CM)","(A,CM)","(C,AM)","(M,AC)","(A,C,M)")
GOF.table = as.data.frame(cbind(model,G2,X2,df,pvalue,AIC))

# Inference concerning conditional associations

summary(fit.homog.assoc)

