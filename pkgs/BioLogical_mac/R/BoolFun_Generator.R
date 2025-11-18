#' @title Generate a specific type of Boolean function
#' @description Generate various specialized types of Boolean functions that 
#' exhibit greater order, efficiency, or redundancy compared to random functions. 
#' Detail definitions are provided in the papers mentioned in \link{BoolFun_Type}.
#' @param bf_type A character included in the following:
#'  \itemize{
#'   \item \code{C}: canalized
#'   \item \code{D}: dominating
#'   \item \code{M}: monotone
#'   \item \code{P}: \eqn{\mathrm{Post}^{2}} (Clique)
#'   \item \code{T}: Threshold
#'   \item \code{R}: Random
#' }
#' @param VarNum An integer between 1 and 15.
#' @param Bias The fraction of "1" in the truth table of a Booelan function. 
#' Note that, except in random cases, the actual bias may differ from the set 
#' value due to the specific function type.
#' @param vars Configure arguments for special type function:
#'  \itemize{
#'  \item an integrate: for \code{c('D','T')}; \code{'1'} or \code{'0'} 
#' for 1/0-type function, other means random set (random 1/0-type).
#'  \item \code{c(int, int)}: for \code{'C'}, the first means layer of 
#' canalization; the second means 1/0-type or random (other number except for 1 and 0)
#'  canalized & canalizing values.
#'  \item \code{c(int, int)}: for \code{'P'}, the first means 1/0-type or 
#' random as \code{c('D','T')} cases; the second determines whether repeated 
#' selection of a state is permitted (non-zero values indicate permission; 
#' zero indicates prohibition).
#'  \item \code{c(int, int)}: for \code{'M'}, the first means 1/0-type or 
#' random as \code{c('D','T')} cases; the second refers to the number of 
#' random seed (between \code{1} and \code{VarNum}).
#' }
#' Detail meaning see References of \link{BoolFun_Type}. 
#' @return An IntegerVector of length \eqn{2^{VarNum}}{2^VarNum}.
#' @export
#' @examples
#' # Generate a random function.
#' # set.seed(2024);
#' # BoolFun_Generator('R', 4L);
#' # Return c(0,1,0,0, 1,0,1,1, 0,1,0,0, 0,0,1,0)
#' set.seed(2024)
#' BoolFun_Generator('R', 4L)
#'
#' # Generate a canalizing function.
#' # set.seed(202401);
#' # BoolFun_Generator('C', 3L);
#' # Return c(1,1,1,1, 1,0,0,0)
#' set.seed(202401)
#' BoolFun_Generator('C', 3L)
#' df_boolfun<-BoolFun_NestedCana(
#'   c(1,1,1,1, 1,0,0,0), PrintOut = TRUE)
#' # if {x_3=0} ==> f(x)=1
#' # else if {x_3=1, x_1=1} ==> f(x)=0
#' # else if {x_3=1, x_1=0, x_2=0} ==> f(x)=1
#' # else if {x_3=1, x_1=0, x_2=1} ==> f(x)=0
#'
#' # Generate a threshold function.
#' # set.seed(202402);
#' # aboolfun=BoolFun_Generator('T', 5L);
#' # aboolfun is c(
#' #   1,1,0,1,0,1,0,0, 
#' #   1,1,1,1,1,1,0,1,
#' #   0,1,0,0,0,0,0,0,
#' #   1,1,0,1,0,1,0,0)
#' # BoolFun_Type(aboolfun, 'T', TRUE);
#' # Return "a_0=1, a_1=-1, a_2=-1, a_3=1, a_4=-1, theta=1", TRUE
#' set.seed(202402)
#' aboolfun=BoolFun_Generator('T', 5L)
#' BoolFun_Type(aboolfun, 'T', TRUE)
#' 
BoolFun_Generator<-function(bf_type=c('R','C','P','M','D','T'), VarNum=3L, Bias=0.50, vars=NULL){
  if(!is.integer(VarNum)||VarNum[1]<=0||VarNum[1]>16){
    stop("Please enter an integer (1~16). If VarNum is a vector, only use the first value.\n");}
  else {
    # Configurate parameters.
    if('C'==bf_type[1]){
      if(is.null(vars)){
        Vars.c=c(1L,9L,rep(0L,2));}
      else if(vars[1]>VarNum){
        stop("Invalid length or other illegal setting of vars. See help.\n");}
      else {
        if(1==length(vars[1])){
          Vars.c=c(as.integer(vars[1]),9L,0L,0L);}
        else {
          Vars.c=c(as.integer(vars[1]),as.integer(vars[2]),0L,0L);}}}
    else if(bf_type[1] %in% c('P','M')){# 0 for 0-type, 1 for 1-type, other or NULL for random
      if(is.null(vars)){
        Vars.c=c(9L,1L,rep(0L,2));}# Second for P_repeated, M_Seed
      else {
        Vars.c=c(as.integer(vars[1]),as.integer(vars[2]),rep(0L,2));}}
    else if(bf_type[1] %in% c('D','T')){# 0 for 0-type, 1 for 1-type, other or NULL for random
      if(is.null(vars)){
        Vars.c=c(9L,rep(0L,3));}
      else {
        Vars.c=c(as.integer(vars[1]),rep(0L,3));}}
    else if('R'==bf_type[1]){
        Vars.c=c(0L,0L,0L,0L);}
    else {
      stop("Illegal function type. Please check help document.\n");}
    resvec=c_BF_Generator(bf_type[1], VarNum, Bias, Vars.c);
    return (resvec);
  }
}

