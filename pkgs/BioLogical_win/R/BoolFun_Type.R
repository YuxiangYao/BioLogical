#' @title Determine whether a function belongs to a special type
#' @description This function can determine whether the given function belongs 
#' to the specified category. See the description of parameter \code{bf_type}.
#' @param Bit_Vec a bool (or 0/1-integral) vector, represents a truth table of Boolean function.
#' @param bf_type a character contained in the following:
#'  \itemize{
#'   \item \code{C}: canalized [\href{https://doi.org/10.1073/pnas.0407783101}{Kauffman,2004}]
#'   \item \code{D}: dominating [\href{https://doi.org/10.1016/S0378-4371(98)00260-X}{Challet,1998}]
#'   \item \code{M}: monotone [\href{https://doi.org/10.1070/RM2003v058n05ABEH000667}{Korshunov,2003}]
#'   \item \code{N}: Spin-like [\href{https://doi.org/10.1088/0305-4470/20/11/009}{Derrida,1987}; \href{https://doi.org/10.1073/pnas.1722609115}{Fontclos,2018}]
#'   \item \code{P}: Post^2 (Clique) [\href{https://doi.org/10.1073/pnas.1534782100}{Pogosyan,1997}]
#'   \item \code{S}: Sign [\href{https://doi.org/10.1007/s11538-008-9304-7}{Aracena,2008}]
#'   \item \code{T}: Threshold [\href{https://doi.org/10.1073/pnas.0305937101}{Li,2004}]
#' }
#' @param ShowRegulation Should show the regulationships of T/N-type?
#' @return a Boolean value
#' @export
#' @examples
#' # Is it a Canalized function? [Canalizing/Canalized value is 1/0]
#' # BoolFun_Type(c(runif(8)>0.5,rep(FALSE,each=8)),'C')
#' # Return a TRUE.
#' BoolFun_Type(c(runif(8)>0.5,rep(FALSE,each=8)),'C')
#'
#' # Is it a Post function? [y=x1*x2+x2*x3+x3*x1] ~ c(0,0,0,1,0,1,1,0)
#' # BoolFun_Type(c(0,0,0,1,0,1,1,0),'P')
#' # Return a TRUE.
#' BoolFun_Type(c(0,0,0,1,0,1,1,0),'P')
#'
#' # Is it a Threshold function? [(-x1+x2-x3+1>0)?TRUE:FALSE]
#' # BoolFun_Type(c(1,0,1,1,0,0,1,0),'T')
#' # Return a TRUE.
#' BoolFun_Type(c(1,0,1,1,0,0,1,0),'T')
#' 
BoolFun_Type<-function(Bit_Vec,bf_type=c('C','P','M','S','N','D','T','E'), ShowRegulation=FALSE){
  if(!bf_type%in%c('C','P','M','S','N','D','T','E')){
    stop("Invalid Boolean function types. See HELP.\n");}
  tmp=r_CheckValidBoolFun(Bit_Vec);
  xx=c_BF_isPointed(tmp$bits,tmp$kin,bf_type,ShowRegulation);
  return(xx);
}