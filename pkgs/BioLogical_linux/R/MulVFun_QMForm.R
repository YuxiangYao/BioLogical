#' @title Show the Quine-McCluskey form of a multi-valued function
#' @description The function employs our native C++ standard multi-valued 
#' Quine-McCluskey algorithm to obtain prime implicants. Because of algorithms 
#' and system properties, the result is not unique for multi-valued system. See reference 
#' [\href{https://doi.org/10.1007/s00500-007-0175-x}{Petrik2008}].
#' @param aVec An IntegerVector that represents a mapping table of multi-valued function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system.
#' @param VarsName A CharacterVector that represents each name of variable.
#' @param ShowType Which output format should be used? Options are restricted to
#' \code{print} (printed in terminal) and \code{data.frame} (in data.frame form). 
#' @param SourceTable A logical value. Should return the source mapping table? 
#' It is ignored when \code{ShowType} is \code{print} (Default: \code{FALSE}).
#' @return Three different cases:
#' \itemize{
#'    \item \code{NULL}: Only print in terminal (\code{ShowType=="print"})
#'    \item \code{data.frame}: \code{ShowType} is \code{"data.frame"} and 
#' \code{SourceTable} is \code{FALSE}.
#'    \item \code{list}: \code{ShowType} is \code{"data.frame"} and 
#' \code{SourceTable} is \code{TRUE}. [[1]] is the result; [[2]] is the 
#' source funciton.}
#' @export 
#' @examples
#' # Calculate the Quine-McCluskey form of multi-valued function.
#' # Analyze a three-valued multiplexer: 
#' #             | a, iff x=0
#' # f(x,c,b,a)= | b, iff x=1
#' #             | c, iff x=2
#' # Example_3v_Multiplexer; # Show the multiplexer.
#' # MulVFun_QMForm(Example_3v_Multiplexer[,5], k=4, L=3,
#' #   VarsName=colnames(Example_3v_Multiplexer)[4:1]);
#' # 
#' # Return its QMC form:
#' # x  c  b  a f(*)
#' # 2  1 NA NA    1
#' # 2  2 NA NA    2
#' # 1 NA  1 NA    1
#' # 1 NA  2 NA    2
#' # 0 NA NA  1    1
#' # 0 NA NA  2    2
#' MulVFun_QMForm(Example_3v_Multiplexer[,5], k=4, L=3, 
#'    VarsName=colnames(Example_3v_Multiplexer)[4:1])
#' 
MulVFun_QMForm<-function(aVec, k, L, VarsName=NULL, ShowType=c("print","data.frame"), SourceTable=FALSE){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  if(all(a_vec==a_vec[1])){
      cat("It is a constant function.\n");
  } else {
    var.names=LETTERS[k:1];
    if(is.null(VarsName)){;}
    else {
      if(length(VarsName)<as.integer(k)){
        stop("Not enough variable's names.\n");
      }else {
        var.names=VarsName[k:1];}}
    qmc_=c_MulF_QuineMcCluskey(as.integer(a_vec),as.integer(k),as.integer(L));}
  qmc_=as.data.frame(qmc_);
  colnames(qmc_)=c(var.names,"f(*)");
  qmc_[qmc_<0]=NA;# wildcard character '*'-'0' is -6; replaced by NA
  if(ShowType[1]=="data.frame"){
    if(SourceTable){
      sous=r_GenerateTruthTable_M(k,L);
      sous=cbind.data.frame(sous,a_vec);
      colnames(sous)=colnames(qmc_);
      return(list(QMC=qmc_,Source=sous));
    } else {
      return(qmc_);
    }
  }
  else if(ShowType[1]=="print"){
    print(qmc_);
    return(NULL);
  } else {
    stop("Invalid arguments of ShowType");
  }
}