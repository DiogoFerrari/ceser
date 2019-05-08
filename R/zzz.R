# required (by devtools) to link the cpp code 
#' @useDynLib ceser
#' @importFrom Rcpp sourceCpp
NULL 

#' @importFrom magrittr %>%



.onAttach<- function(libname, pkgname) 
{
 packageStartupMessage('
 ## -----------------------------------
 ## ceser package
 ## -----------------------------------
 Author(s): Diogo Ferrari, John Jackson
 ')
} 

.onUnload <- function (libpath) {library.dynam.unload("ceser", libpath)}

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))
