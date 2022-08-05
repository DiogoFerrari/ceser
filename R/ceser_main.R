

## {{{ ancilary functions }}}

check.type <- function(type)
{
    current.implemented = c("HC0", "HC1", "HC2", "HC3", "HC4")
    if (type %!in% c("HC0", "HC1", "HC2", "HC3", "HC4" )) {
        stop("\n\nParameter 'type' must be one of: 'HC0', 'HC1', 'HC2', 'HC3', 'HC4'")
    }
    if (type %!in% current.implemented) {
        stop(paste0("\n\nCurrent implementation of the package only includes the following HC corrections: ", paste0(current.implemented, collapse=", ") ) )
    }
    invisible(NULL)
}
get.residuals  <- function(mod, P, type, n, k)
{
    ## -----------------------------
    ## Selecting levarege correction
    ## ----------------------------- 
    ## leverage correction Davidson and MacKinnon (1993, p. 554)
    ## "hc0"    White SE
    ## "hc1"    Stata's Default
    ## "hc2"    Unbiased under homoskedasticity
    ## "hc3"    Unbiased under homoskedasticity, Default (conservative), Davidson and MacKinnon (1993, p. 554) 
    ## if(is.null(type)) {type <- "HC3"}
    ## type <- match.arg(type, c("HC2", "HC3"))

    e   = stats::residuals(mod)
    hii = diag(P)
    if (!is.null(type)) {
        if (type=="HC0") {e = e}
        if (type=="HC1") {e = e * ( sqrt(n/(n-k)) )}
        if (type=="HC2") {e = e * ( 1/sqrt(1 - hii) )}
        if (type=="HC3") {e = e * ( 1/(1 - hii) )}
        if (type=="HC4") {
            deltai = purrr::map_dbl(.x= hii*(n/k), function(x) min(4, x))
            v      = purrr::map2_dbl(.x=hii, .y=deltai, function(hii, deltai) (1-hii)^deltai)
            e      = e * ( 1/sqrt(v) )
        }
    }
    return(e)
}
check.clusters <- function(mod, cluster)
{
    if (! (class(cluster) %in% c("formula", "character") | is.null(cluster)) ) {
        stop("\n\nThe values of the parameter 'cluster' of the function vcovCESE() must be either a string vector with the names of the variables to cluster, a formula in which the RHS contains the variables to clusters, or 'NULL'. See documentation of vcovCESE().")
    }
    if (class(cluster) == 'character') {cluster = paste0(" ~ ", paste0(cluster, collapse = " + ") )  %>% stats::as.formula(.)}
    if (!is.null(cluster)) {
        tryCatch(
        {
            existing.colnames <- stats::expand.model.frame(mod, cluster, na.expand = FALSE) %>% names
        },
        error = function(e) stop("\n\nCheck the names of the variables used to cluster the standard error in the parameter 'cluster' of the function vcovCESE(). It seems some names are not in the data.")
        )
    }
    invisible(NULL)
}
get.clusters <- function(mod, cluster, n)
{
    if (is.null(cluster)) cluster <- attr(mod, "cluster")    ## select cluster variable. if cluster is not provided (NULL), it can be an attribute of the model. If so, get that info
    if (is.null(cluster)) cluster <- 1L:n                  ## if cluster is not supplied and no attribute have that info, use observarion level cluster
    if(!inherits(cluster, "formula")) {                     ## inherits return TRUE if cluster is a formula
        cluster.formula = paste0("~", paste0(cluster, collapse=" + " ) )  %>% stats::as.formula(.)
    }else{
        cluster.formula = cluster 
    }
    cluster.formula.no.intercept = stats::update(cluster.formula, ~ .- 1)
    
    cluster_tmp <- stats::expand.model.frame(mod, cluster.formula, na.expand = FALSE)
    cluster     <- stats::model.frame(cluster.formula, cluster_tmp, na.action = stats::na.pass)
    
    ## handle omitted or excluded observations
    if((n != nrow(cluster)) && !is.null(mod$na.action) && (class(mod$na.action) %in% c("exclude", "omit"))) {
        cluster <- cluster[-mod$na.action, , drop = FALSE]
    }
    if(nrow(cluster) != n) stop("\n\nThe number of observations in 'cluster' and in the data used in the regression do not match. Check the NA values.")

    return(list(cluster.formula.no.intercept=cluster.formula.no.intercept, cluster=cluster))
}
get.Pg <- function(Pf, ddp, n)
{
    Pg    = matrix(0, nrow=n, ncol=n) # it is P in .do file
    for (i in 1:n)
    {
        for (j in i:n)
        {
            Pg[j, i] = Pf[j, i] * ddp[j, i]
            Pg[i, j] = Pg[j, i]
        }
    }
    return(Pg)
}
get.Zg <- function(Z, ddp, n)
{
    Zg = matrix(0, ncol=n, nrow=n)
    for (i in 1:n)
    {
        for (j in i:n)
        {
            Zg[j, i] = Z[j, i] * ddp[j, i]
            Zg[i, j] = Zg[j, i]
        }
    }
    return(Zg)
}
get.sumZg <- function(Zg, n)
{
    sumZg = 0
    for (i in 1:n)
    {
        j.ini = (i+1)
        if (j.ini<=n) {
            for (j in j.ini:n)
            {
                sumZg = sumZg + Zg[j, i]
            }
        }
    }
    return(sumZg)
}
getQ <- function(Q1, Q2, ng, nclusters)
{
    ## Computing q1Tq1, q2Tq2, and q1Tq2 for each cluster
    q1_2 = 0
    q1q2 = 0
    q2_2 = 0

    t      = 0
    ng.vec = ng %>% dplyr::pull(.)
    for (k in 1:nclusters)
    {
        i.ini = t + 1
        i.end = t + ng.vec[t+1]
        for (i in i.ini:i.end)
        {
            j.ini = i
            j.end = t + ng.vec[t+1]
            for (j in j.ini:j.end)
            {
                tryCatch(
                q1_2 = q1_2 + Q1[j, i]^2
                   ,
                error = function(e) {stop("\n\nError when getting matrix Q. Make sure you estimated the model with a data.frame ordered by the clustering variables. See documentation.\n\n")},
                )
                q1q2 = q1q2 + Q1[j, i]*Q2[j, i]
                q2_2 = q2_2 + Q2[j, i]^2
            }
        }
        t = t + ng.vec[t+1]
    }
    return(list(q1_2=q1_2, q1q2=q1q2, q2_2=q2_2))
}
get.Qe <- function(Q1, Q2, eep, ng, nclusters)
{
    Qe     = matrix(0, nrow=2, ncol=1)
    t      = 0
    ng.vec = ng %>% dplyr::pull(.)
    for (k in 1:nclusters)
    {
        i.ini = t + 1
        i.end = t + ng.vec[t+1]
        for (i in i.ini:i.end)
        {
            j.ini = i
            j.end = t + ng.vec[t+1]
            for (j in j.ini:j.end)
            {
                Qe[1,1] = Qe[1,1] + eep[j,i]*Q1[j,i]
                Qe[2,1] = Qe[2,1] + eep[j,i]*Q2[j,i]
            }
        }
        t = t + ng.vec[t+1]
    }
    return(Qe)
}

## }}}

## {{{ docs }}}

#' Cluster Estimated Standard Errors
#'
#' Cluster Estimated Standard Errors (CESE)
#'
#'
#' @param mod a model object. It can be the output of the functions \code{lm}, \code{glm}, or other regression function that returns compatible objects. Important: See details
#' @param cluster either a string vector with the name of the variables that will be used to cluster the standard errors, or a formula - e.g., ~ rhs, with a summation of the variables that will be used to cluster the standard errors replacing the \code{rhs} -, or a vector, matrix, or data.frame with the clustering data.
#' @param type string with either \code{HC0}, \code{HC1}, \code{HC2}, \code{HC3}, or \code{HC4}. It specifies the type of heteroskedasticity correction to use  (see Davidson and MacKinnon (1993) and Hayes and Cai (2007)).
#'
#' @return The function returns a variance-covariace matrix of the coefficient estimates using the Cluster Estimated Standard Error (CESE) method.
#'
#' @details
#' The data frame must be ordered by the clustering variables before 
#' estimating the model with \code{lm}, etc.
#' 
#' @references
#' Jackson, John (2019) Corrected Standard Errors with Clustered Data. Political Analysis.
#' 
#' Hayes, A. F., & Cai, L., (2007) Using heteroskedasticity-consistent standard error estimators in ols regression: an introduction and software implementation, Behavior research methods, 39(4), 709–722.
#' 
#' Davidson, R., & MacKinnon, J. G., (2004) Econometric theory and methods: Oxford University Press New York.
#' 
#' @examples
#' 
#' mod  = lm(enep ~  enpc + fapres + enpcfapres + proximity + eneg + logmag + logmag_eneg , data=dcese)
#'
#' ## --------------------------------------
#' ## Getting the variance covariance matrix
#' ## -------------------------------------- 
#' ## Original variance-covariance matrix (no clustered std. errors)
#' vcov(mod)
#' 
#' ## Variance-covariance matrix using CRSE (sandwish package)
#' ## sandwich::vcovCL(mod, cluster = ~ country)
#' ## sandwich::vcovCL(mod, cluster = ~ country, type="HC3")
#' 
#' ## Variance-covariance matrix using CESE
#' ceser::vcovCESE(mod, cluster = ~ country)
#' ceser::vcovCESE(mod, cluster = ~ country, type="HC3") # HC3 correction
#'
#' ## ---------
#' ## Summaries
#' ## ---------
#' ## no robust SE 
#' summary(mod)                                                                          
#' 
#' ## summary table using CRSE (sandwich package)
#' ## lmtest::coeftest(mod, vcov = sandwich::vcovCL, cluster = ~ country)                   
#'
#' ## summary using CESE
#' lmtest::coeftest(mod, vcov = ceser::vcovCESE, cluster = ~ country, type='HC3')
#' 
#' 
#' @export

## }}}
vcovCESE <- function(mod, cluster = NULL, type=NULL)
{
    X         = stats::model.matrix(mod) %>% tibble::as_tibble(.) ## data used in the regression
    n         = nrow(X)                                    ## size of the data used in the regression

    ## options(warn=-1)
    ## on.exit(options(warn=0))
    
    ## -------------
    ## error control
    ## -------------
    if (!is.null(type)) {check.type(type)}
    check.clusters(mod, cluster)

    ## ----------------
    ## get cluster info (returns a list)
    ## ----------------
    cluster.info                 = get.clusters(mod, cluster, n)
    cluster                      = cluster.info$cluster
    cluster.formula.no.intercept = cluster.info$cluster.formula.no.intercept
    cluster.formula              = cluster.formula.no.intercept ## change it later if intercept should be used


    ## ------------------
    ## Computing the CESE
    ## ------------------
    x         = X %>% as.matrix()
    k         = ncol(x)
    Sigma.hat = matrix(NA, nrow = k, ncol=k, dimnames=list(names(x),names(x)))
    ng        = cluster %>% dplyr::group_by_all() %>% dplyr::mutate(ng = dplyr::n())  %>% dplyr::ungroup(.)  %>% dplyr::select(ng) 
    nclusters = cluster %>% dplyr::filter(!duplicated(.)) %>% nrow(.)
    err       = matrix(0, nrow=1,ncol=3)

    I    = diag(n)
    iota = stats::model.matrix(cluster.formula, cluster)
    ddp  = iota %*% t(iota)

    xpx  = t(x) %*% x
    xpxi = solve(xpx)
    Pf   = x %*% xpxi %*% t(x)
    M    = matrix(diag(I - Pf) , nrow=n)
    xiix = t(x) %*% ddp %*% x
    Pg   = get.Pg(Pf, ddp, n)

    Z  = x %*% xpxi %*% xiix %*% xpxi %*% t(x)
    Zg = get.Zg(Z, ddp, n)
    sumZg = get.sumZg(Zg, n)

    Q1 = I - Pg
    Q2 = ddp - Q1 - Pg %*% ddp - ddp %*% Pg + Zg
    
    tryCatch(
        q    = getQ(Q1, Q2, ng, nclusters)
       ,
	error = function(e) {stop("\n\nError when getting matrix Q. Make sure you estimated the model with a data.frame ordered by the clustering variables. See documentation.\n\n")},
    )
    q1_2 = q$q1_2
    q1q2 = q$q1q2
    q2_2 = q$q2_2

    ## compute squares and cross-products matrix and its inverse
    z  = matrix(c(q1_2, q1q2, q1q2, q2_2), ncol=2, nrow=2, byrow=F)
    zi = solve(z)

    ## residuals
    e = get.residuals(mod, Pg, type, n, k)
    
    ## cross-products
    eep = e %*% t(e)
    
    ## matrices for accumulating Q'e
    tryCatch(
        Qe = get.Qe(Q1, Q2, eep, ng, nclusters)
       ,
	error = function(e) {stop("\n\nError when getting matrix Qe. Make sure you estimated the model with a data.frame ordered by the clustering variables. See documentation.\n\n")},
    )

    ## compute estimates for sigma and rho (Sigma.hat)
    sighat = zi %*% Qe 

    ## -----------
    ## Corrections
    ## ----------- 
    ## err = 1 if rho > sigma
    err[1,1] =  1*(sighat[2,1] >= sighat[1,1])
    ## make adjustment if err == 1
    sighat[1,1] = (sighat[2,1]< sighat[1,1])*sighat[1,1] + (sighat[2,1] >= sighat[1,1])*(sighat[2,1]+.02) 

    ## compute estimated standard errors from sighat
    sigv       = sighat[2,1]*ddp
    diag(sigv) = sighat[1,1]
    sb = xpxi %*% t(x) %*% sigv %*% x %*%  xpxi

    return(sb)
}

