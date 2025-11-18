#' @title Decode the information of a multi-valued signed function
#' @description The function analyzes and returns a list containing information
#' whether given function is signed, and corresponding weight matrix, and baselines.
#' It is also compatible with Boolean system (L=2).
#' @param aVec An IntegerVector that represents a mapping table of multi-valued 
#' function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system.
#' @param PrintOut A logical value. Should output be displayed in the terminal? 
#' (Default: \code{TRUE}).
#' @return A list containing three elements:
#' \itemize{
#' \item \code{[[1]] IsSigned}: Is Signed? 1 for Yes, 0 for No.
#' \item \code{[[2]] Weight}: The weights of all pairwise relationships among
#' varibales; it is a L*(kL) matrix.
#' \item \code{[[3]] Delta}: The baselines of variables, ordered from low 
#' to high bits.}
#' @export
#' @examples
#' # Generate a domainted function and analyze it belong to signed class or not.
#' # set.seed(1001);
#' # A_MulV_Fun=MulVFun_Generator("D", k=3L, L=3L, MappingTable=TRUE);
#' # A_MulV_Fun_Sign=MulVFun_is_Signed(A_MulV_Fun[,4], 3, 3, TRUE);
#' # View(A_MulV_Fun); # Show the mapping table 
#' #       v2 v1 v0 f_out
#' #  [1,]  0  0  0     2
#' #  [2,]  0  0  1     2
#' #  [3,]  0  0  2     1
#' #  [4,]  0  1  0     2
#' #  [5,]  0  1  1     0
#' #  [6,]  0  1  2     0
#' #          ... ... 
#' # Check the "signed" regulatory patterns (User can test other cases)
#' # A_MulV_Fun_Sign[[2]]%*%c(1,0,0, 1,0,0, 1,0,0) ==> (-4,-1,0) ==> f_out(0,0,0)=2;
#' # A_MulV_Fun_Sign[[2]]%*%c(0,0,1, 1,0,0, 1,0,0) ==> (-2,-1,-5) ==> f_out(0,0,2)=1;
#' # A_MulV_Fun_Sign[[2]]%*%c(0,1,0, 0,0,1, 1,0,0) ==> (4,-1,1) ==> f_out(0,2,1)=0;
#' # ... ... 
#' 
#' set.seed(1001)
#' A_MulV_Fun=MulVFun_Generator("D", k=3L, L=3L, MappingTable=TRUE)
#' A_MulV_Fun_Sign=MulVFun_is_Signed(A_MulV_Fun[,4], 3, 3, TRUE)
#' A_MulV_Fun_Sign[[2]]%*%c(1,0,0, 1,0,0, 1,0,0)
#' A_MulV_Fun_Sign[[2]]%*%c(0,0,1, 1,0,0, 1,0,0)
#' A_MulV_Fun_Sign[[2]]%*%c(0,1,0, 0,0,1, 1,0,0)
#' 
MulVFun_is_Signed<-function(aVec, k, L, PrintOut=TRUE){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  xx=c_M_Signed(a_vec, as.integer(k), as.integer(L), as.integer(L^k));
  if(xx[[1]]<0){# Not a domainted function!
    if(PrintOut){
      cat("Given function is non-signed.");}
  } else if(xx[[1]]==0){
    if(PrintOut){
      cat("Given function is unknown due to large undetermined weights.");}
  }else {
    xx[[2]]=t(xx[[2]]);
    tmp.names=paste0("v",c(0:(L-1)))
    rownames(xx[[2]])=tmp.names;
    colnames(xx[[2]])=paste0(rep(paste0("k",c(1:k),"_"),each=L),
      rep(tmp.names,time=k));
    names(xx)=c("IsSigned","Weight","Delta");
    # if(PrintOut){
    #   wenzi=paste0("f(~): {w_",k-1,"=",res[[1]][k],", ");
    #   if((k-1)>2){
    #     for(ii in c((k-1):2)){
    #       wenzi=paste0(wenzi,"w_",ii-1,"=",res[[1]][ii],", ");
    #     }
    #   }
    #   cat(paste0(wenzi,"w_0=",res[[1]][1]," | "));
    #   wenzi=NULL
    #   for(ii in c(1:(L-1))){
    #     wenzi=paste0(wenzi,"delta_",ii-1,"=",res[[2]][ii],", ");
    #   }
    #   wenzi=paste0(wenzi,"delta_",L-1,"=",res[[2]][L],"}\n");
    #   cat(wenzi);
    # }
  }
  return (xx);
}

