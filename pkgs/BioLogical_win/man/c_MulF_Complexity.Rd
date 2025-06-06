% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{c_MulF_Complexity}
\alias{c_MulF_Complexity}
\title{C++ prototype function for MulVFun_Complexity()}
\usage{
c_MulF_Complexity(avec, k, L)
}
\arguments{
\item{avec, }{IntegerVector, each element in {0,1,...L-1}}

\item{k, }{input number (in-degree)}

\item{L, }{system of L:{v0,v1,...v_{L-1}}}
}
\value{
numeric, the complexity of given function.
}
\description{
C++ prototype function for MulVFun_Complexity()
}
\details{
Please see document of \link{MulVFun_Complexity}
}
