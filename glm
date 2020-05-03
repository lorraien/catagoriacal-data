# read in the snoring/heart diseaese data set

snore = read.table("snoring.txt",header=T)

# check that is is a data.frame

is.data.frame(snore)

# run glm

snore.glm = glm(cbind(yes,no) ~ snore,data=snore,family=binomial())

# alternatively can do it with wieghts

snore.glm2 = glm(yes/(yes+no) ~ snore,data=snore,family=binomial(),weights=yes+no)

# print out a summary of the fit

summary(snore.glm)

# compute the sample proportions

snore.sample.props = snore$yes/(snore$yes+snore$no)

# extract the coefficients 

snore.glm$coeff # or use the function 'coef'
coef(snore.glm)

# Extract the fitted values from the glm fit (note: for binary data the fitted values are model predicted proportions/probabilities not estimated frequencies)

snore.glm$fitted  #or use the function 'fitted'
fitted(snore.glm)

# fitted and coeff are attributes of a glm object that is returned by the function 'glm'
# any attribute can be extracted by a the name of the object followed by a $ sign followed by the attribute
# To find the attributes of 'snore.glm'

attributes(snore.glm)

# summary also has attributes that can be extracted

attributes(summary(snore.glm))

# plot the sample proportions

plot(snore$snore,snore.sample.props,xlab="Level of Snoring",ylab ="Predicted Proportions",main="Snoring/Heart Disease Example",ylim=c(0,.15))

#add the predicted proportions

lines(snore$snore,fitted(snore.glm),type="o")

#want a smoother fit?  use 'predict.glm'

snore.predict = predict.glm(snore.glm,type="response",newdata=data.frame(snore=seq(0,5,.1)))

lines(seq(0,5,.1),snore.predict,col="red")

#want the linear responses at the 4 levels of snoring

predict.glm(snore.glm,type="link")


# fit a probit model to the data

snore.probit.fit = glm(cbind(yes,no) ~ snore,data=snore,family=binomial(link="probit"))

summary(snore.probit.fit)

# get the coefficients

coef(snore.probit.fit)

# Get the predict proportions, rounded to 3 digits

round(fitted(snore.probit.fit),3)

# add the fitted proportions to the plot

lines(snore$snore,fitted(snore.probit.fit),col="blue",type="o")

# We've talke about three types of residuals in lecture: deviance residuals, Pearson residuals and 
# standarized Pearson residuals.  Each of these are easily obtained from a glm object---

# want deviance residuals from the logit model?

residuals.glm(snore.glm,type="deviance")

# summary also has an attribute named deviance.resid

summary(snore.glm)$deviance.resid

# how about Pearson residuals?

snore.pres = residuals.glm(snore.glm,type="pearson")
snore.pres

# what of standardized residuals? ---Just a little more work

snore.pres/sqrt(1-influence(snore.glm)$hat)

# influence is a function that provides basic quantities used in diagnostics for regression models

influence(snore.glm)

# also check out 'influence.measures' which gives regression deletion diagnostics

influence.measures(snore.glm)

# and diagnostic plots

plot(snore.glm)


#####################################

# Poisson example

# Import crab data

crab = read.table("crab.txt",header=T)

# Fit Poisson model
crab.glm = glm(Sat~Width, family=poisson(link="log"), data=crab)

# Get summary of poisson model
summary(crab.glm)

# Now fit quasi-Poisson model
crab.glmq = glm(Sat~Width, family=quasipoisson(link="log"), data=crab)

# Get summary of poisson model
summary(crab.glmq)
