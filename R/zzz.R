## # required (by devtools) to link the cpp code 
## #' @useDynLib ceser
## #' @importFrom Rcpp sourceCpp
## NULL 

#' @importFrom magrittr %>%



.onAttach<- function(libname, pkgname) 
{
 packageStartupMessage('
 ## ------------------------------------------
 ## Cluster Estimated Standard Errors  (ceser)
 ## ------------------------------------------
 Author(s): Diogo Ferrari and John Jackson
 ')
} 


if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
