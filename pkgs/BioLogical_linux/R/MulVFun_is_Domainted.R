#' @title Decode the information of a multi-valued domaineted function
#' @description The function analyzes and returns a list containing information
#' on whether the given function is domaineted and corresponding weights and 
#' baselines (deltas).
#' It is also compatible with Boolean system (L=2).
#' @param aVec An IntegerVector that represents a mapping table of multi-valued 
#' function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system.
#' @param PrintOut A logical value. Should the output be displayed in the terminal?
#' @return A list containing three elements: 
#' \itemize{
#' \item \code{[[1]] IsDomaint}: Is domainted? 1 for Yes, 0 for No.
#' \item \code{[[2]] Weight}: The weights of variables, ordered from low to high bits.
#' \item \code{[[3]] Delta}: The baselines of variables, ordered from low 
#' to high bits.}
#' @export
#' @examples
#' # Show the information of a multi-valued linear threshold function.
#' # set.seed(1002);
#' # A_MulV_Fun=MulVFun_Generator("D", k=2L, L=3L, MappingTable=TRUE);
#' # View(A_MulV_Fun); # Show the mapping table.
#' # A_MulV_Fun_Domaint=MulVFun_is_Domainted(A_MulV_Fun[,3], 2, 3, TRUE);
#' # Terminal output: 
#' # f(~): {w_1=-5, w_0=3 | delta_0=-1, delta_1=0, delta_2=1}
#' # Above info means: avg max{i| \sum(w_i*x_i)+delta_i }
#' # A_MulV_Fun_Domaint # Show the information
#' # $IsDomaint
#' # [1] 1
#' # $Weight
#' # [1] 3 -5
#' # $Delta
#' # [1] -1  0  1
#' 
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

