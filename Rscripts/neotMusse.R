Path_base <- "~/Documents/Salamanders/"
library(phytools)
library(BAMMtools)
library(diversitree)

# Read in neoteny scores
setwd(paste(Path_base, "neoteny", sep=""))
neoteny <- read.csv("neoteny.csv", stringsAsFactors=F)
spNames <- apply(neoteny, 1, function(x) paste(x[1], x[2], sep="_", collapse=""))
neoteny <- neoteny[,3]
names(neoteny) <- spNames


par(mfrow=c(1,2))

# Read in tree and BAMM spec/ext run
setwd(paste(Path_base, "bamm/extant_only", sep=""))
extree <- read.tree("modernSal.tre")
exEdata <- getEventData(extree, "ex_event_data.txt", burnin=0.1, nsamples=200)

neot_vec <- neoteny + 2
neot_vec <- neot_vec[exEdata$tip.label]
Cols <- c("black", "gray", "red", "green")
x <- plot(exEdata, breaksmethod="quantile", method="polar", spex="netdiv")
tiplabels(pch=16, col=Cols[neot_vec], cex=0.5)

ex_Sp <- traitDependentBAMM(exEdata, neoteny[exEdata$tip.label], reps=1e3, rate='s')
ex_Ex <- traitDependentBAMM(exEdata, neot_vec, reps=1e3, rate='e')
ex_ND <- traitDependentBAMM(exEdata, neoteny[exEdata$tip.label], reps=1e3, rate='n')

# Read in fossil tree and BAMM spec/ext run
setwd(paste(Path_base, "bamm/output-1", sep=""))
ftree <- read.tree("fullTree.tre")
fEdata <- getEventData(ftree, "full_event_data.txt", burnin=0.1, nsamples=200)

neot_vec <- neoteny + 2
neot_vec <- neot_vec[fEdata$tip.label]
Cols <- c("black", "gray", "red", "green")
plot(fEdata, colorbreaks=x$colorbreaks, method="polar", spex="netdiv")
tiplabels(pch=21, bg=Cols[neot_vec], cex=0.5)

f_Sp <- traitDependentBAMM(fEdata, neot_vec, reps=1e3, rate='s')
f_Ex <- traitDependentBAMM(fEdata, neot_vec, reps=1e3, rate='e')
f_ND <- traitDependentBAMM(fEdata, neot_vec, reps=1e3, rate='n')



# Pull the tip rates from the event data
tipSpec <- exEdata$meanTipLambda
names(tipSpec) <- exEdata$tip.label
tipExt <- exEdata$meanTipMu
names(tipExt) <- exEdata$tip.label

# Find the mean speciation/extinction by neotenic state
meanLambdas <- tapply(tipSpec, neoteny[names(tipSpec)], mean)
meanMus <- tapply(tipExt, neoteny[names(tipExt)], mean)


# set up MuSSE pars that assert an ordered transition regime
mussePars <- c(0.1, 0.15, 0.2,		# lambda 1, 2, 3
				0.03, 0.045, 0.06,	# mu 1, 2, 3
				0.05, 0,				# q12, q13
				0.05, 0.05,			# q21, q23
				0, 0.05)				# q31, q32
musseNames <- c("l1", "l2", "l3",
				"m1", "m2", "m3",
				"q12", "q13",
				"q21", "q23",
				"q31", "q32")

# Build likelihood function and use the spec/ext rates from BAMM to start
lik.m <- make.musse(extree, neoteny[extree$tip.label]+1, k=4)
p <- starting.point.musse(extree, 3)
p2 <- p
p2[1:3] <- meanLambdas
p2[4:6] <- meanMus

# maximum likelihood estimate using starting parameters & BAMM parameters
fit.m <- find.mle(lik.m, p[argnames(lik.m)])
fit.m2 <- find.mle(lik.m, p2[argnames(lik.m)])

prior <- make.prior.exponential( 1 / (2 * (p2[1] - p2[4])))
samples <- mcmc(lik.m, coef(fit.m2), nstep=10000, w=1, prior=prior, print.every=100)

col <- c("blue", "orange", "red")
profiles.plot(samples[2:4], col)



#### QUASSE

lik.q <- make.quasse(extree, liabMeans[extree$tip.label], liabSD[extree$tip.label], sigmoid.x, sigmoid.x)

p <- starting.point.quasse(extree, liabMeans, liabSD)

# ignore drift
lik.q2 <- constrain(lik.q, drift ~ 0)

p.start <- c(p[1], p[1], mean(liabMeans), 1, p[2], p[2], mean(liabMeans), 1, p[3])
names(p.start) <- argnames(lik.q2)

lower <- c(0, 0, min(liabMeans), -Inf, 0, 0, min(liabMeans), -Inf, 0)
fit <- find.mle(lik.q2, p.start, control=list(parscale=0.1, lower=lower))

lik.constant <- constrain(lik.q2, l.y1 ~ l.y0, l.xmid ~ 0, lr ~ 1, m.y1 ~ m.y0, m.xmid ~ 0, mr ~ 1)
fit.constant <- find.mle(lik.constant, p.start[argnames(lik.constant)], control=list(parscale=0.1), lower=0, verbose=0)

anova(fit, constant=fit.constant)
