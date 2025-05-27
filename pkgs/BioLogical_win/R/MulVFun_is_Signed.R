#' @title Decode the information of a multi-valued signed function
#' @description This function can return a \code{list} records information of
#' whether given function is, the weight matrix, and baselines.
#' @param aVec an IntegerVector, a mapping table of function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system.
#' @param PrintOut a Boolean value, show output abnormal message to the 
#' terminal. (Default: TRUE)
#' @return a 3-element list (if it is signed): [[1]] An integer [[2]] describes 
#' weights of each input variable (L*(kL) matrix), [[3]] records deltas of L 
#' values (See examples).
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
#' # Check the "signed" regulatory pattern (User can test other cases.)
#' # A_MulV_Fun_Sign[[2]]%*%c(1,0,0, 1,0,0, 1,0,0) ==> (-4,-1,0) ==> f_out(0,0,0)=2;
#' # A_MulV_Fun_Sign[[2]]%*%c(0,0,1, 1,0,0, 1,0,0) ==> (-2,-1,-5) ==> f_out(0,0,2)=1;
#' # A_MulV_Fun_Sign[[2]]%*%c(0,1,0, 0,0,1, 1,0,0) ==> (4,-1,1) ==> f_out(0,2,1)=0;
#' # ... ... 
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

