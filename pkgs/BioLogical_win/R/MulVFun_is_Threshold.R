#' @title Decode the information of a multi-valued linear threshold function
#' @description This function can return a \code{list} records information of
#' whether given function is threshold, weights, baseline and threshold.
#' It is also compatible with Boolean system (L=2).
#' @param aVec an IntegerVector, a mapping table of function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system.
#' @param PrintOut a Boolean value, show output to the terminal.
#' @return a 3-element list (if it is threshold): [[1]] An integer [[2]] describes 
#' weights of each input variable (from low to high), [[3]] records the L-1 boundaries of mapping 
#' values (See examples).
#' @export
#' @examples
#' # Show the information of a multi-valued linear threshold function.
#' # set.seed(2025);
#' # A_MulV_Fun=MulVFun_Generator("T", k=3L, L=3L, MappingTable=TRUE)
#' # View(A_MulV_Fun)# Show the mapping table.
#' # A_MulV_Fun_Thres=MulVFun_is_Threshold(A_MulV_Fun[,4], 3, 3, TRUE);
#' # Return: 
#' # f(~)=2*x_2-1*x_1-2*x_0
#' # Threshold intervals for L values: (-inf,-3], (-3,0], (0,+inf)
#' set.seed(2025);
#' A_MulV_Fun=MulVFun_Generator("T", k=3L, L=3L, MappingTable=TRUE)
#' A_MulV_Fun_Thres=MulVFun_is_Threshold(A_MulV_Fun[,4], 3, 3, TRUE);
#' 
#' # Single input Quaternary funciton:
#' # A_MulV_Fun_Thres=MulVFun_is_Threshold(c(1,1,2,3), 1, 4, TRUE);
#' # Return:
#' # f(~)=1*x_0
#' # Threshold intervals for L values: (-inf,-1], (-1,1], (1,2], (2,+inf)
#' A_MulV_Fun_Thres=MulVFun_is_Threshold(c(1,1,2,3), 1, 4, TRUE);
#' 
MulVFun_is_Threshold<-function(aVec, k, L, PrintOut=FALSE){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  xx=c_M_Threshold(a_vec, as.integer(k), as.integer(L), as.integer(L^k));
  if(xx[[1]]<0){# Not a threshold function!
    if(PrintOut){
      cat("Given function is non-threshold.");}
  } else {
    if(PrintOut){
      wenzi="f(~)=";pp=k;
      while(1){# Fisrt non-zero element noe need "+" symbol
        if(xx[[2]][pp]!=0){
          wenzi=paste0(wenzi,xx[[2]][pp],"*x_",pp-1);
          break;
        }
        pp=pp-1;
      }
      if(2<=(pp)){
        for(ii in c((pp-1):1)){
          if(xx[[2]][ii]>0){
            wenzi=paste0(wenzi,"+",xx[[2]][ii],"*x_",ii-1);
          } else if(xx[[2]][ii]<0){
            wenzi=paste0(wenzi,xx[[2]][ii],"*x_",ii-1);
          } else {# Do none.
            ;# Input weight is zero.
          }
        }
      }
      cat(wenzi);
      wenzi="\nThreshold intervals for L values: (-inf,"
      for(ii in c(1:(L-1))){
        wenzi=paste0(wenzi,xx[[3]][ii],"], (",xx[[3]][ii],",");
      }
      wenzi=paste0(wenzi,"+inf)\n"); 
      cat(wenzi);
    }
    names(xx)=c("IsThreshold","Weight","Boundary");
  }
  return (xx);
}
