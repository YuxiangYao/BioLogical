#' @title Decode the information of a multi-valued domaineted function
#' @description This function can return a \code{list} records information of
#' whether given function is threshold, weights, baselines (deltas).
#' It is also compatible with Boolean system (L=2).
#' @param aVec an IntegerVector, a mapping table of function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system.
#' @param PrintOut a Boolean value, show output to the terminal.
#' @return a 3-element list (if it is threshold): [[1]] An integer [[2]] describes 
#' weights of each input variable (from low to high), [[3]] records the 
#' L baselines of the mapping table (See examples).
#' @export
#' @examples
#' # Show the information of a multi-valued linear threshold function.
#' # Return a list.
#' # set.seed(1002);
#' # A_MulV_Fun=MulVFun_Generator("D", k=2L, L=3L, MappingTable=TRUE);
#' # View(A_MulV_Fun); # Show the mapping table.
#' # A_MulV_Fun_Domaint=MulVFun_is_Domainted(A_MulV_Fun[,3], 2, 3, TRUE);
#' # Terminal output: 
#' # f(~): {w_1=-5, w_0=3 | delta_0=-1, delta_1=0, delta_2=1}
#' set.seed(1002)
#' A_MulV_Fun=MulVFun_Generator("D", k=2L, L=3L, MappingTable=TRUE)
#' A_MulV_Fun_Domaint=MulVFun_is_Domainted(A_MulV_Fun[,3], 2, 3, TRUE)
#' 
MulVFun_is_Domainted<-function(aVec, k, L, PrintOut=FALSE){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  xx=c_M_Domainted(a_vec, as.integer(k), as.integer(L), as.integer(L^k));
  if(xx[[1]]<0){# Not a domainted function!
    if(PrintOut){
      cat("Given function is non-domainted.");}
  } else {
    names(xx)=c("IsDomaint","Weight","Delta");
    if(PrintOut){
      wenzi=paste0("f(~): {w_",k-1,"=",xx[[2]][k],", ");
      if((k-1)>2){
        for(ii in c((k-1):2)){
          wenzi=paste0(wenzi,"w_",ii-1,"=",xx[[2]][ii],", ");
        }
      }
      cat(paste0(wenzi,"w_0=",xx[[2]][1]," | "));
      wenzi=NULL;
      for(ii in c(1:(L-1))){
        wenzi=paste0(wenzi,"delta_",ii-1,"=",xx[[3]][ii],", ");
      }
      wenzi=paste0(wenzi,"delta_",L-1,"=",xx[[3]][L],"}\n");
      cat(wenzi);
    }
  }
  return (xx);
}

