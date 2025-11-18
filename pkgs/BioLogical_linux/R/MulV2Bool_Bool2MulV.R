#' @title Boolean to multi-valued and multi-valued to Boolean transformation
#' @param aVec An IntegerVector that represents a mapping table of multi-valued function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system (Default: 3).
#' @param Thres An IntegerVector containing thresholds for a multi-valued system. 
#' Its length is \code{k+1}, with each element ranging from 1 to L-1. If not provided, 
#' the function randomly generate them (Default: \code{NA}).
#' @param Bool2MulV A logical value. Is the original function Boolean? If \code{TRUE} 
#' the transformation involves converting it from Boolean to multi-valued. 
#' Otherwise, from multi-valued to Boolean. (Default: \code{TRUE}).
#' @param MappingTable A logical value. Is the mapping table also returned? 
#' (Default: \code{FALSE}) 
#' @return An IntegerVector (or IntegerMatrix if \code{MappingTable} is \code{TRUE}), 
#' denoting "mapping table" of target system.
#' @export
#' @examples
#' # Transform a Boolean function to a multi-valued one.
#' # The Boolean is f(x,y)=x & y, a possible result see follows:
#' # set.seed(2025);
#' # MulV2Bool_Bool2MulV(c(0,0,0,1), k=2, L=3, 
#' #   c(1,1,1), <--- the two inputs and one output share same threholds:
#' #                  in this case,  0 mapto 0 and 1 mapto {1,2}
#' #   MappingTable = TRUE);
#' # Return
#' #       v1 v0 xx   #  Boolean system
#' #  [1,]  0  0  0   #    0  0  0 
#' #  [2,]  0  1  0   #    0  1  0 
#' #  [3,]  0  2  0   #    0  1  0 
#' #  [4,]  1  0  0   #    1  0  0 
#' #  [5,]  1  1  2   #    1  1  1 
#' #  [6,]  1  2  2   #    1  1  1 
#' #  [7,]  2  0  0   #    1  0  0 
#' #  [8,]  2  1  1   #    1  1  1 
#' #  [9,]  2  2  2   #    1  1  1 
#' set.seed(2025)
#' MulV2Bool_Bool2MulV(c(0,0,0,1), 2, 3, c(1,1,1), MappingTable=TRUE)
#' 
MulV2Bool_Bool2MulV<-function(aVec, k=0, L=3, Thres=NA, Bool2MulV=TRUE, 
  MappingTable=FALSE){
  if(is.na(Thres[1])){
    thres=sample.int((L-1), as.integer(k+1), replace=TRUE);
  } else {
    if((k+1)!=length(Thres)){
      stop("Should give k+1 thresholds of multi-valued system\n");
    }
    if(any(Thres<=0)||any(Thres>=L)){
      stop("Each threshold should be [1,L-1].\n");
    }
    thres=as.integer(Thres);
  }
  if(Bool2MulV){# a Boolen 2 Multi-valued
    a_fun=r_CheckValidBoolFun(aVec)[[1]];
  } else {
    a_fun=r_CheckValidMulVFun(aVec, as.integer(k), as.integer(L));
  }
  xx=c_MulV2Bool_Bool2MulV(a_fun, as.integer(k), as.integer(L), 
    thres, Bool2MulV);
  if(MappingTable){
    if(Bool2MulV){
      xx=cbind(r_GenerateTruthTable_M(k,L), xx);
    } else {
      xx=cbind(r_GenerateTruthTable(k, FALSE), xx);
    }
  }
  return(xx);
}