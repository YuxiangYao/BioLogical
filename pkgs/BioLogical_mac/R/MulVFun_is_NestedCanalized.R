#' @title Analyze the structure of a multi-valued nested canalized function
#' @description This function can return a \code{data.frame} of nested canalized 
#' structure. Due to algorithm, same level of canalzing variable will return the 
#' samller variable's ID (See the second example). It is also compatible with 
#' Boolean system, while here recommend \code{BoolFun_NestedCanalized} for its 
#' binary feature.
#' @param aVec an integer vector, represents a truth table of Boolean function.
#' @param k an integer, k-variable of function.
#' @param L an integer, L-level multi-valued system.
#' @param PrintOut a Boolean value, show output to the terminal.
#' @return a 2-element list (if it is (nested) canalized): [[1]] An integer [[2]] A 
#' dataframe descirbes the nested structure of a canalized function. (See examples).
#' @export
#' @examples
#' # Show the canalized structure of a multi-valued nested canalized function.
#' # Return a data.frame (Also show nested "if-else-if" hierarchical relations).
#' # set.seed(2025);
#' # A_MulV_Fun=MulVFun_Generator("C", k=4L, L=3L, CanaDeep=3, 
#' #   CanaVar=c(2,1,3), CanaLayerInfo=list(list(c(2,0),c(0),c(1,2)), 
#' #   list(c(0,1),c(2),c(2,0))), MappingTable=TRUE);
#' # View(A_MulV_Fun);# Show the mapping table.
#' # A_MulV_Fun_Neseted=MulVFun_is_NestedCanalized(A_MulV_Fun[,5], 4, 3, TRUE);
#' # Terminal output: 
#' # if:
#' # x_1=0 ==> f(x)=1
#' # x_1=2 ==> f(x)=0
#' # else if:
#' #  ..., x_0=0 ==> f(x)=2
#' # else if:
#' #  ...,  ..., x_2=1 ==> f(x)=2
#' #  ...,  ..., x_2=2 ==> f(x)=0
#' # View(A_MulV_Fun_Neseted); # Show the nested structure.
#' set.seed(2025);
#' A_MulV_Fun=MulVFun_Generator("C", k=4L, L=3L, CanaDeep=3, 
#'   CanaVar=c(2,1,3), CanaLayerInfo=list(list(c(2,0),c(0),c(1,2)), 
#'   list(c(0,1),c(2),c(2,0))), MappingTable=TRUE);
#' A_MulV_Fun_Neseted=MulVFun_is_NestedCanalized(A_MulV_Fun[,5], 4, 3, TRUE);
#' 
#' # Second example for same level issue.
#' # set.seed(2020);
#' # A_Bool_Fun=MulVFun_Generator("C", k=4L, L=2L, CanaDeep=2, 
#' #   CanaVar=c(2,3), CanaLayerInfo=list(list(c(0),c(1)), list(c(1),c(1))), 
#' #   MappingTable = TRUE);
#' # Please compare the mapping table and printed information.
#' # A_Bool_Fun_Neseted=MulVFun_is_NestedCanalized(A_Bool_Fun[,5],4,2,TRUE);
#' # if:
#' # x_0=1 ==> f(x)=1
#' # else if:
#' #  ..., x_1=1 ==> f(x)=1
#' # else if:
#' #  ...,  ..., x_2=0 ==> f(x)=1
#' # else if:
#' #  ...,  ...,  ..., x_3=0 ==> f(x)=0
#' #  ...,  ...,  ..., x_3=1 ==> f(x)=1
#' # In fact, here set the first layer is f(x_1=0) ==> f(x)=1, 
#' # it is also coupled with f(x_0=0) ==> f(x)=1
#' set.seed(2020)
#' A_Bool_Fun=MulVFun_Generator("C", k=4L, L=2L, CanaDeep=2, CanaVar=c(2,3), 
#'   CanaLayerInfo=list(list(c(0),c(1)), list(c(1),c(1))), MappingTable = TRUE)
#' A_Bool_Fun_Neseted=MulVFun_is_NestedCanalized(A_Bool_Fun[,5],4,2,TRUE)
#' 
MulVFun_is_NestedCanalized<-function(aVec, k, L, PrintOut=FALSE){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  xx=c_M_NestedCanalized(a_vec, as.integer(k), as.integer(L), as.integer(L^k));
  if(xx[[1]]<0){# Not a canalized function!
    if(PrintOut){
      cat("Given function is non-canalized.");}
  } else {
    res=NULL;
    tmp=rep('*',k+1);
    for(ii in c(1:length(xx[[2]]))){# How many variables can be canalized?
      for(jj in c(1:length(xx[[3]][[ii]]))){# IN!! & OUT!!
        tmp[k-xx[[2]][ii]]=xx[[3]][[ii]][jj];
        tmp[k+1]=xx[[4]][[ii]][jj];
        res=rbind(res,tmp);
      }
      tmp[k-xx[[2]][ii]]="Or";
    }
    colnames(res)=c(paste0("x_",c((k-1):0)),"Out");
    rownames(res)=NULL;
    res=as.data.frame(res);
    # Should output detail information?
    if(PrintOut){
      outers=NULL;
      for(ii in c(1:length(xx[[2]]))){# Record the ID of variables.
        if(ii==1){
          outers=paste0(outers,"if:\n"); 
        }
        else {
          outers=paste0(outers,"else if:\n"); 
        }
        for(jj in c(1:length(xx[[3]][[ii]]))){
            outers=paste0(outers,
              paste(rep(" ..., ",ii-1),collapse = ""),
              "x_",xx[[2]][ii],"=",xx[[3]][[ii]][jj]," ==> f(x)=",xx[[4]][[ii]][jj],"\n");
        }
      }
      cat(outers);
    }
    xx=list(IsCana=1L, CanaTopo=res);
  }
  return(xx);
}
