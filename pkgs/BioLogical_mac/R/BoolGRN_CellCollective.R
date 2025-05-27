#' @title 70 genetic networks selected from Cell Collective Web
#' @description List, contains all network. Each element is an independent network.
#' The number is coded by the website. For example, \code{c_11863} is the biological network, No.11863.
#' One network is also organized as a \code{List}:
#'  \itemize{
#'   \item \code{AllMember}: \code{CharacterVector}, all genes' name
#'   \item \code{InEdge}: \code{IntegerVector}, Input regulations by order. No input is denoted as \code{NA}. In each vector, first element is highest bit.
#'   \item \code{OutEdge}: \code{IntegerVector}, Output regulations by order. No output is denoted as \code{NA}.
#'   \item \code{BoolFun}: \code{IntegerVector}, Boolean function of each genes (length=2^length(\code{InEdge})). If \code{InEdge} is \code{NA}, it also keeps \code{NA}.
#' }
#' Detail information, please see \href{https://cellcollective.org}{Cell Collective Website}.
#' @usage
#' data(BoolGRN_CellCollective)
#' @format
#' An list of many networks.
"BoolGRN_CellCollective"
