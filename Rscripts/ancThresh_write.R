require(phytools)
likLiab <- phytools:::likLiab
probMatch <- phytools:::probMatch
to.vector <- phytools:::to.vector
logPrior <- phytools:::logPrior
ancThresh <- function (tree, x, ngen = 1000, sequence = NULL, method = "mcmc", 
    model = c("BM", "OU", "lambda"), control = list(), ...) 
{
    if (!inherits(tree, "phylo")) 
        stop("tree should be an object of class \"phylo\".")
    if (method != "mcmc") 
        stop(paste(c("do not recognize method =", method, ",quitting")))
    model <- model[1]
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    if (is.matrix(x)) {
        X <- x[tree$tip.label, ]
        if (is.null(sequence)) {
            message("**** NOTE: no sequence provided, column order in x")
            seq <- colnames(X)
        }
        else seq <- sequence
    }
    else if (is.vector(x)) {
        x <- x[tree$tip.label]
        if (is.null(sequence)) {
            message("**** NOTE: no sequence provided, using alphabetical or numerical order")
            seq <- sort(levels(as.factor(x)))
        }
        else seq <- sequence
        X <- to.matrix(x, seq)
    }
    X <- X/apply(X, 1, sum)
    X <- X[, seq]
    th <- c(1:length(seq)) - 1
    names(th) <- seq
    x <- to.vector(X)
    l <- sapply(x, function(x) runif(n = 1, min = th[x] - 1, 
        max = th[x]))
    if (model == "OU") 
        alpha <- 0.1 * max(nodeHeights(tree))
    if (model == "lambda") 
        lambda <- 1
    n <- length(tree$tip)
    m <- length(th)
    npar <- if (model == "BM") 
        tree$Nnode + n + m - 2
    else tree$Nnode + n + m - 1
    PrA <- matrix(1/m, tree$Nnode, m, dimnames = list(1:tree$Nnode + 
        n, seq))
    if (!is.null(control$pr.anc)) {
        if (!is.matrix(control$pr.anc)) {
            message("**** NOTE: prior on ancestral states must be in matrix form; using default prior")
            control$pr.anc <- NULL
        }
        else {
            control$pr.anc <- control$pr.anc[, seq, drop = FALSE]
            PrA[rownames(control$pr.anc), ] <- control$pr.anc
            control$pr.anc <- PrA
        }
    }
    con = list(sample = 1000, propliab = 0.5 * max(nodeHeights(tree)), 
        propthresh = 0.05 * max(nodeHeights(tree)), propalpha = 0.1 * 
            max(nodeHeights(tree)), proplambda = 0.01, pr.anc = PrA, 
        pr.th = 0.01, burnin = round(0.2 * ngen), plot = TRUE, 
        print = TRUE, piecol = setNames(palette()[1:length(seq)], 
            seq), tipcol = "input", quiet = FALSE)
    con[(namc <- names(control))] <- control
    con <- con[!sapply(con, is.null)]
    temp <- apply(con$pr.anc, 1, rstate)
    a <- sapply(temp, function(x) runif(n = 1, min = th[x] - 
        1, max = th[x]))
    th[length(th)] <- Inf
    V <- if (model == "BM") 
        vcvPhylo(tree)
    else if (model == "OU") 
        vcvPhylo(tree, model = "OU", alpha = alpha)
    else if (model == "lambda") 
        vcvPhylo(tree, model = "lambda", lambda = lambda)
    if (any(tree$edge.length <= (10 * .Machine$double.eps))) 
        stop("some branch lengths are 0 or nearly zero")
    invV <- solve(V)
    detV <- determinant(V, logarithm = TRUE)$modulus[1]
    lik1 <- likLiab(l, a, V, invV, detV) + log(probMatch(X, l, 
        th, seq))
    A <- matrix(NA, ngen/con$sample + 1, tree$Nnode, dimnames = list(NULL, 
        n + 1:tree$Nnode))
    B <- if (model == "BM") 
        matrix(NA, ngen/con$sample + 1, m + 2, dimnames = list(NULL, 
            c("gen", names(th), "logLik")))
    else if (model == "OU") 
        matrix(NA, ngen/con$sample + 1, m + 3, dimnames = list(NULL, 
            c("gen", names(th), "alpha", "logLik")))
    else if (model == "lambda") 
        matrix(NA, ngen/con$sample + 1, m + 3, dimnames = list(NULL, 
            c("gen", names(th), "lambda", "logLik")))
    C <- matrix(NA, ngen/con$sample + 1, tree$Nnode + n, dimnames = list(NULL, 
        c(tree$tip.label, 1:tree$Nnode + n)))
    A[1, ] <- sapply(a, threshState, thresholds = th)
    B[1, ] <- if (model == "BM") 
        c(0, th, lik1)
    else if (model == "OU") 
        c(0, th, alpha, lik1)
    else if (model == "lambda") 
        c(0, th, lambda, lik1)
    C[1, ] <- c(l[tree$tip.label], a[as.character(1:tree$Nnode + 
        n)])
    message("MCMC starting....")
    logL <- lik1 <- likLiab(l, a, V, invV, detV) + log(probMatch(X, 
        l, th, seq))
    for (i in 1:ngen) {
        lik1 <- logL
        d <- i%%npar
        if (ngen >= 1000) 
            if (i%%1000 == 0) 
                if (con$print) 
                  message(paste("gen", i))
        ap <- a
        lp <- l
        thp <- th
        if (model == "OU") 
            alphap <- alpha
        if (model == "lambda") 
            lambdap <- lambda
        Vp <- V
        invVp <- invV
        detVp <- detV
        if (d <= tree$Nnode && d != 0) {
            ind <- d%%tree$Nnode
            if (ind == 0) 
                ind <- tree$Nnode
            ap[ind] <- a[ind] + rnorm(n = 1, sd = sqrt(con$propliab))
        }
        else {
            if ((d > tree$Nnode && d <= (tree$Nnode + n)) || 
                (npar == (tree$Nnode + n) && d == 0)) {
                if (d == 0) 
                  ind <- n
                else {
                  ind <- (d - tree$Nnode)%%n
                  if (ind == 0) 
                    ind <- n
                }
                lp[ind] <- l[ind] + rnorm(n = 1, sd = sqrt(con$propliab))
            }
            else if (d > (tree$Nnode + n) && d <= (tree$Nnode + 
                n + m - 2) || (npar == (tree$Nnode + n + m - 
                2) && d == 0)) {
                if (d) 
                  ind <- (d - tree$Nnode - n)%%m + 1
                else ind <- m - 1
                thp[ind] <- bounce(th[ind], rnorm(n = 1, sd = sqrt(con$propthresh)), 
                  c(th[ind - 1], th[ind + 1]))
            }
            else {
                if (model == "OU") {
                  alphap <- bounce(alpha, rnorm(n = 1, sd = sqrt(con$propalpha)), 
                    c(0, Inf))
                  Vp <- vcvPhylo(tree, model = "OU", alpha = alphap)
                }
                else if (model == "lambda") {
                  lambdap <- bounce(lambda, rnorm(n = 1, sd = sqrt(con$proplambda)), 
                    c(0, 1))
                  Vp <- vcvPhylo(tree, model = "lambda", lambda = lambdap)
                }
                invVp <- solve(Vp)
                detVp <- determinant(Vp, logarithm = TRUE)$modulus[1]
            }
        }
        lik2 <- likLiab(lp, ap, Vp, invVp, detVp) + log(probMatch(X, 
            lp, thp, seq))
        p.odds <- min(c(1, exp(lik2 + logPrior(sapply(ap, threshState, 
            thresholds = thp), thp, con) - lik1 - logPrior(sapply(a, 
            threshState, thresholds = th), th, con))))
        if (p.odds > runif(n = 1)) {
            a <- ap
            l <- lp
            th <- thp
            V <- Vp
            detV <- detVp
            invV <- invVp
            if (model == "OU") 
                alpha <- alphap
            if (model == "lambda") 
                lambda <- lambdap
            logL <- lik2
        }
        else logL <- lik1
        if (i%%con$sample == 0) {
            A[i/con$sample + 1, ] <- sapply(a, threshState, thresholds = th)
            B[i/con$sample + 1, ] <- if (model == "BM") 
                c(i, th[colnames(B)[1 + 1:m]], logL)
            else if (model == "OU") 
                c(i, th[colnames(B)[1 + 1:m]], alpha, logL)
            else if (model == "lambda") 
                c(i, th[colnames(B)[1 + 1:m]], lambda, logL)
            C[i/con$sample + 1, ] <- c(l[tree$tip.label], a[as.character(1:tree$Nnode + 
                n)])
 	       write.csv(as.data.frame(A), "mcmc_out.csv", quote=F)
 	       write.csv(as.data.frame(B), "param_out.csv", quote=F)
 	       write.csv(as.data.frame(C), "liab_out.csv", quote=F)
        }
    }
    mcmc <- as.data.frame(A)
    param <- as.data.frame(B)
    liab <- as.data.frame(C)
    ace <- matrix(0, tree$Nnode, m, dimnames = list(colnames(A), 
        seq))
    burnin <- which(param[, "gen"] == con$burnin)
    for (i in 1:tree$Nnode) {
        temp <- summary(mcmc[burnin:nrow(mcmc), i])/(nrow(mcmc) - 
            burnin + 1)
        ace[i, names(temp)] <- temp
    }
    if (con$plot) 
        plotThresh(tree, X, list(ace = ace, mcmc = mcmc, par = param, 
            liab = liab), burnin = con$burnin, piecol = con$piecol, 
            tipcol = con$tipcol, ...)
    return(list(ace = ace, mcmc = mcmc, par = param, liab = liab))
}