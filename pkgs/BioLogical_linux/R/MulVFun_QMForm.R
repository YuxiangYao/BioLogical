#' @title Show a multi-valued function's Quine-McCluskey form
#' @description The function employs our native C++ standard multi-valued 
#' Quine-McCluskey algorithm to obtain prime implicants. Because of algorithms 
#' and system properties, the result is not unique. See reference 
#' [\href{https://doi.org/10.1007/s00500-007-0175-x}{Petrik,2008}].
#' @param aVec a vector, represents a truth table of valued function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system.
#' @param VarsName a character vector, represents each variable's name
#' @param ShowType \code{c("print","data.frame")}, limited to \code{print} 
#' (Print in terminal) or \code{data.frame} (as a data.frame form).
#' @param SourceTable logical, should return the source mapping table (Default: 
#' \code{FALSE})? If \code{TRUE}, return a list: [[1]] QMC form, [[2]] source 
#' table. It is invalid when \code{ShowType} is \code{"print"}.
#' @return dataframe/list/NULL, show the multi-valued QMC form. 
#' @export
#' @examples
#' # Calculate the Quine-McCluskey form of multi-valued function.
#' # Analyze a three-valued multiplexer: 
#' #             | a, iff x=0
#' # f(x,c,b,a)= | b, iff x=1
#' #             | c, iff x=2
#' # three_valued_Multiplexer; # Show the multiplexer.
#' # MulVFun_QMForm(three_valued_Multiplexer[,5], k=4,L=3,
#' #   VarsName=colnames(three_valued_Multiplexer)[4:1]);# Return the QMC form.
#' # Return: 
#' # x  c  b  a f(*)
#' # 2  1 NA NA    1
#' # 2  2 NA NA    2
#' # 1 NA  1 NA    1
#' # 1 NA  2 NA    2
#' # 0 NA NA  1    1
#' # 0 NA NA  2    2
#' MulVFun_QMForm(three_valued_Multiplexer[,5], k=4,L=3, 
#'    VarsName=colnames(three_valued_Multiplexer)[4:1]);
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
  qmc_[qmc_==-1]=NA;# "-1" replaced by NA instead of wildcard character "*"
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
  } else {
    stop("Invalid arguments of ShowType");
  }
}