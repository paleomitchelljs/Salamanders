### Function for placing fossil taxon on an existing phylogeny

# Eqn 2 from Bapst 2013 MEE
# calculate the probability of sampling a fossil
#	lam is speciation
#	mu is extinction
#	psi is preservation
calcPs <- function(lam, mu, psi)	{
	res <- ((lam + mu + psi)^2) - (4 * lam * mu)
	res <- 1 - (((lam + mu + psi) - (sqrt(res)))/(2 * lam))
	return(res)
}

# Eqn 1 from Bapst 2013 MEE
# probability density of some amount of missing (unobserved) evolutionary history (delta, del here)
#	del is the amount of missing history
#	lam is speciation
#	mu is extinction
#	psi is preservation
probDel <- function(del, lam, mu, psi)	{
	Ps <- calcPs(lam, mu, psi)
	Rate <- psi + (Ps * lam)
	out <- dgamma(del, shape=2, rate=Rate)
	return(out)
}


# Function requires: tree, a phylogeny
#					age, time from before-present for fossil
#					Name, tip label for fossil taxon
#					Hs, matrix of node heights (generated if not supplied)
#					taxa, vector of tip labels for the crown-clade the fossil is stem to (can be one tip) OR:
#					node, the node number for the crown clade fossil is stem to.
#					grain, how finely the branch length is estimated. Also minimum size branch is allowed to be.
#					maxAge, if the fossil is stem to the *whole clade*, this is the oldest the fossil could be.
#					rates, a vector of speciation, extinction and preservation rates
placeFossil <- function(tree, age, Name="fossil", Hs=NULL, taxa=NULL, node=NULL, grain=0.01, maxAge=NULL, rates=NULL)	{
	require(phytools)
	if (is.null(rates))	{
		stop("You must provide a speciation, extinction and preservation rate for the tree.")
	}
	if (is.null(taxa) && is.null(node))	{
		stop("Please supply either a node or set of crown-group taxa.")
	}
	if (length(taxa) > 1 && is.null(node))	{
		node <- findMRCA(tree, tips=taxa)
	}
	if (length(taxa) == 1 && is.null(node))	{
		node <- which(tree$tip.label == taxa)
	}
	if (is.null(Hs))	{
		Hs <- nodeHeights(tree)
	}

	# Node the fossil is stem to
	dNode <- node
	
	rownames(Hs) <- tree$edge[,2]
	treeAge <- max(Hs)
	
	# node heights are from-root. 
	# Fossil ages are before-present. 
	# Need to convert.
	t_age <- treeAge - age

	# is the fossil stem to the entire clade?	
	if ((dNode - Ntip(tree)) == 1)	{
		if (is.null(maxAge))	{
			stop("Please supply a maximum age if the fossil is stem to the whole clade.")
		}
		if (maxAge < treeAge)	{
			warning("maxAge less than the age of the tree. Setting maxAge to sum of the two.")
			maxAge <- maxAge + treeAge
		}
		beforeRoot <- maxAge - treeAge
		possAges <- seq(from=grain, to=beforeRoot, by=grain)

		Probs <- sapply(t_age + possAges, probDel, lam=rates[1], mu=rates[2], psi=rates[3])		
		branchPt <- sample(possAges, 1, prob=Probs)
		brLen <- t_age + branchPt
		
		tree2 <- tree
		tree2$root.edge <- beforeRoot
		tree2 <- bind.tip(tree2, Name, edge.length=brLen, where=Ntip(tree)+1, position=branchPt)
	}
	else	{
		# Age of the crown clade (fossil must diverge before this)
		dNode_age <- Hs[as.character(dNode),2]
		
		# Age of the total clade (fossil must diverge after this)
		aNode_age <- Hs[as.character(dNode),1]

		if (t_age < aNode_age)	{
			Msg <- paste("Tip-age of ", Name, " is older than the total-clade...Aborted!")
			stop(Msg)
		}
		
		before_dNode_start <- grain
		
		if (t_age < dNode_age)	{
			before_dNode_end <- t_age - aNode_age
			buffer <- dNode_age - t_age
		}
		else	{
			before_dNode_end <- dNode_age - aNode_age
			buffer <- 0			
		}
		# Smallest possible missing history is set to grain
		possAges <- seq(from=before_dNode_start, to=before_dNode_end, by=grain)
		
		# Use eqns from Bapst 2013 MEE to weight the potential missing histories
		Probs <- sapply(t_age + possAges, probDel, lam=rates[1], mu=rates[2], psi=rates[3])
		
		# Amount of missing history
		branchPt <- sample(possAges, 1, prob=Probs) + buffer
		
		# amt of time before dNode_age
		brLen <- t_age - (dNode_age - branchPt)
	
		tree2 <- bind.tip(tree, Name, edge.length=brLen, where=dNode, position=branchPt)
	}
	if (Ntip(tree2) == Nnode(tree2)+1)	{
		return(tree2)
	}
	else	{
		placeFossil(tree=tree, age=age, Name=Name, Hs=Hs, taxa=taxa, node=dNode, grain=grain, maxAge=maxAge, rates=rates)
	}
}

# Code to test the above function
library(phytools)
tree <- pbtree(b=0.05, d=0, t=5)
Hs <- nodeHeights(tree)
maxAge <- max(Hs) + 1
age <- 0.1
Name="fossil-Tip"
taxa=tree$tip.label[2]
rates=c(0.05,0,0.01)
#test <- placeFossil(tree, age=2.5, taxa=tree$tip.label[c(1,2)], Hs=Heights, maxAge=max(Heights)+1, rates=c(0.05,0,0.01))
#test2 <- placeFossil(test, age=2.5, Name="stem", taxa=test$tip.label[c(1,Ntip(test))], maxAge=max(Heights)+1, rates=c(0.05,0,0.01))
#test3 <- placeFossil(tree, age=0.1, Name="fossil-Tip", taxa=tree$tip.label[2], maxAge=max(Heights)+1, rates=c(0.05,0,0.01))
#par(mfrow=c(1,4))
#plotTree(tree)
#plotTree(test)
#plotTree(test2)
#plotTree(test3)