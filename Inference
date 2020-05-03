# This oddsratio function taken from package "vcd" 
# Only change is the option to compute/not compute with continuity correction
# help can be found after loading the package "vcd"

library("vcd")

oddsratio2 = function (x, stratum = NULL, log = TRUE,correct=TRUE, conf.level = 0.95) 
{
    l <- length(dim(x))
    if (l > 2 && is.null(stratum)) 
        stratum <- 3:l
    if (l - length(stratum) > 2) 
        stop("All but 2 dimensions must be specified as strata.")
    if (l == 2 && dim(x) != c(2, 2)) 
        stop("Not a 2 x 2 - table.")
    if (!is.null(stratum) && dim(x)[-stratum] != c(2, 2)) 
        stop("Need strata of 2 x 2 - tables.")
    lor <- function(y) {
        if (correct)
          y <- y + 0.5
        or <- y[1, 1] * y[2, 2]/y[1, 2]/y[2, 1]
        if (log) 
            log(or)
        else or
    }
    ase <- function(y) sqrt(sum(1/(y + 0.5)))
    if (is.null(stratum)) {
        LOR <- lor(x)
        ASE <- ase(x)
    }
    else {
        LOR <- apply(x, stratum, lor)
        ASE <- apply(x, stratum, ase)
    }
    I <- ASE * qnorm((1 + conf.level)/2)
    Z <- LOR/ASE
    list(Estimate=LOR, ASE = if (log) 
        ASE, lwr = if (log) 
        LOR - I
    else exp(log(LOR) - I), upr = if (log) 
        LOR + I
    else exp(log(LOR) + I), Z = if (log) 
        Z, P = if (log) 
        1 - pnorm(abs(Z)), log = log, class = "oddsratio")
}


#ODD RATIO CI EXAMPLE 
heart.attack = c("Fatal","Non-Fatal")
treatment = c("Placebo","Aspirin")
HA.table = expand.grid(heartAttack=heart.attack,Treat=treatment)
HA.data = c(189,10845,104,10933)
HA.table = cbind(HA.table,count=HA.data)
HA.xtab = xtabs(count ~ Treat + heartAttack,data=HA.table)

oddsratio2(HA.xtab,correct=F)
oddsratio2(HA.xtab,correct=F,log=F)

oddsratio2(HA.xtab,correct=T)
oddsratio2(HA.xtab,correct=T,log=F)


#DIFFERENCE IN PROPORTION CI EXAMPLE

row.mar = margin.table(HA.xtab,1)
HA.prop = HA.xtab[,1]/row.mar
diff.prop = HA.prop[1]-HA.prop[2]
se.diff.prop = sqrt(HA.prop[1]*(1-HA.prop[1])/row.mar[1] + HA.prop[2]*(1-HA.prop[2])/row.mar[2])
(CI.HA.diff = diff.prop + c(-1,1)*qnorm(0.975)*se.diff.prop)
z.wald = diff.prop/se.diff.prop # null hypothesis that the diff in prop = 0.
(p.value = 2*(1-pnorm(abs(z.wald))))  # null hypothesis that the diff in prop = 0.


# RELATIVE RISK CI EXAMPLE

r = HA.prop[1]/HA.prop[2]
log.r = log(r)
se.logr = sqrt((1-HA.prop[1])/(HA.prop[1]*row.mar[1]) + (1-HA.prop[2])/(HA.prop[2]*row.mar[2]))
(CI.logr = log.r + c(-1,1)*qnorm(0.975)*se.logr)
(CI.r = exp(CI.logr))
z.wald = log.r/se.logr # null hypothesis that the relative risk=1.
(p.value = 2*(1-pnorm(abs(z.wald))))  # null hypothesis that the relative risk=1.


