as.phylo.bammdata <- BAMMtools:::as.phylo.bammdata
NU.branching.times <- BAMMtools:::NU.branching.times
exponentialRate <- BAMMtools:::exponentialRate

getRTT <- function (ephy, tvec = NULL, nslices = 100, 
    node = NULL, nodetype = "include") 
{
    if (!"bammdata" %in% class(ephy)) {
        stop("Object ephy must be of class 'bammdata'\n")
    }
    if (is.null(node)) {
        nodeset <- c(length(ephy$tip.label) + 1, ephy$edge[, 
            2])
    }
    else if (!is.null(node) & nodetype == "include") {
        nodeset <- unlist(sapply(node, function(x) getDesc(ephy, 
            x)$desc_set))
    }
    else if (!is.null(node) & nodetype == "exclude") {
        nodeset <- setdiff(ephy$edge[, 2], unlist(sapply(node, 
            function(x) getDesc(ephy, x)$desc_set)))
        nodeset <- c(length(ephy$tip.label) + 1, nodeset)
    }
    else {
        stop("error in getRateThroughTimeMatrix\n")
    }
    if (is.ultrametric(as.phylo.bammdata(ephy))) {
        bt <- branching.times(as.phylo.bammdata(ephy))
    }
    if (!is.ultrametric(as.phylo.bammdata(ephy))) {
        bt <- NU.branching.times(as.phylo.bammdata(ephy))
    }
    maxpossible <- max(bt[as.character(intersect(nodeset, ephy$edge[, 
        1]))])
     if (!is.null(node)) {
        if (nodetype == "include") {
            nodeset <- nodeset[nodeset != node]
        }
    }
    if (is.null(tvec))	{
    		tvec <- seq(0, max(bt), length.out = nslices)
    }
    tol <- 0.00001
    goodTime <- function(vec, val, tol) {
        (vec[, 2] - val <= tol) & (val - vec[, 3] <= tol)
    }
    getRates <- function(time, es, events, isGoodNode) {
        isGoodTime <- goodTime(es, time, tol = tol)
        if (!(is.null(node))) {
            isGoodNode <- es[, 1] %in% nodeset
        }
        estemp <- es[isGoodTime & isGoodNode, ]
        if (is.vector(estemp)) {
            index <- estemp[4]
        }
        else {
            index <- estemp[, 4]
        }
        tvv <- time - events$time[index]
        rates <- exponentialRate(tvv, events$lam1[index], events$lam2[index])
        return(list(rates, index))
    }
    bySample <- function(counter, ephy) {
        es <- ephy$eventBranchSegs[[counter]]
        events <- ephy$eventData[[counter]]
        isGoodNode <- rep(TRUE, nrow(es))
        ret <- lapply(tvec, function(x) getRates(time = x, es, 
            events, isGoodNode))
        mmRow <- unlist(lapply(ret, function(x) mean(x[[1]])))
        mumatRow <- unlist(lapply(ret, function(x) mean(events$mu1[x[[2]]])))
        return(list(mmRow, mumatRow))
    }
    ret <- lapply(1:length(ephy$eventBranchSegs), function(y) bySample(y, 
        ephy))
    mm <- lapply(ret, function(x) x[[1]])
    mm <- do.call(rbind, mm)
    mumat <- lapply(ret, function(x) x[[2]])
    mumat <- do.call(rbind, mumat)
    obj <- list()
    if (ephy$type == "diversification") {
        obj$lambda <- mm
        obj$mu <- mumat
    }
    if (ephy$type == "trait") {
        obj$beta <- mm
    }
    obj$times <- tvec
    class(obj) <- "bamm-ratematrix"
    if (ephy$type == "diversification") {
        obj$type = "diversification"
    }
    else {
        obj$type = "trait"
    }
    return(obj)
}
