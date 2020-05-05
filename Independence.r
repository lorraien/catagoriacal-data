# Load needed library
library(vcd)
library(epitools)
library(rmeta)

# X^2 test for independence example

relig.belief = c("fund","mod","lib")
degree = c("<HS","HS or Jr.C.",">= Bachelor")
table3.2 = expand.grid(relig.belief=relig.belief,degree=degree)
count.relig = c(178,138,108,570,648,442,138,252,252)
table3.2 = cbind(table3.2,count=count.relig)
relig.table = xtabs(count ~ degree+relig.belief,table3.2)

relig.chisq = chisq.test(relig.table)
attributes(relig.chisq)
relig.chisq$observed
relig.chisq$expected
relig.chisq$statistic
relig.chisq$p.value  

#alternatively you can use package 'vcd' from CRAN, load it and use the function assoc.stats to get both X^2 and G^2
summary(assocstats(relig.table))

# Computing standarized Pearson residuals
#  chisq.test gives the Pearson residuals

relig.chisq$residuals

n = sum(relig.table)
pi = rowSums(relig.table)/n  
pj = colSums(relig.table)/n
(es = relig.chisq$residuals/sqrt(outer(1-pi,1-pj)))


# Trend test example in lecture 7 using income and job satisfaction

# Create data
income = c("<15k","12k-15k","25k-40k",">40k")
jobstat = c(1,2,3,4)
table3.2 = expand.grid(jobstat=jobstat , income=income)
counts = c(1,3,10,6,2,3,10,7,1,6,14,12,0,1,9,11)
table3.2 = cbind(table3.2,count=counts)
table1 = xtabs(count ~ income+jobstat,table3.2)
table1=as.matrix(table1)

income = c(rep(7.5,20), rep(20,22), rep(32.5,33), rep(60,21))

job=apply(table1, 1, function(x){rep(c(1,2,3,4),x)})

data=cbind(income, as.numeric(unlist(job)))

Msquared = cor(data)[2,1]^2*(nrow(data)-1)
1-pchisq(Msquared,1)

 
#Fisher's exact test

convicted = c("yes","no")
type = c("dizygotic","monzygotic")
counts = c(2,15,10,3)
convict.table = cbind(expand.grid(convicted = convicted,type=type),counts)
convict.xtab = xtabs(counts~ type+convicted,data=convict.table)
fisher.test(convict.xtab,alternative="less")

# Maentel Haenszel example
# Create death penalty data
d.penalty = array( c( 415,38,54,12,17,140,1,5),dim = c(2, 2, 2),dimnames = list( Defendant.Race = c("White","Black"),Verdict = c("No.Death", "Death"),Victim.Race = c("White","Black")))

# Conduct Maentel Haenszel test
mantelhaen.test(d.penalty)


# Test for heterogeneity with the Breslow-Day test
n.black = c( sum( d.penalty[,,1][2,] ), sum( d.penalty[,,2][2,] ) ) 
n.white = c( sum( d.penalty[,,1][1,] ), sum( d.penalty[,,2][1,] ) ) 
b.death = c( sum( d.penalty[,,1][2,2] ), sum( d.penalty[,,2][2,2] ) ) 
w.death = c( sum( d.penalty[,,1][1,2] ), sum( d.penalty[,,2][1,2] ) )

mh.rslt = meta.MH( n.black, n.white, b.death, w.death, names=c("White Victim", "Black Victim" ) )
summary( mh.rslt )
