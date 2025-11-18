#' @title Calculate the input/edge effetiveness of a multi-valued function
#' @description The analysis depends on multi-valued Quine-McCluskey method 
#' [\href{https://doi.org/10.1007/s00500-007-0175-x}{Petrik,2008}] to obtain the 
#' prime implicants of input vectors belong to sets "\eqn{f(x)=0}", "\eqn{f(x)=1}",
#' ..., "\eqn{f(x)=i}", ..., "\eqn{f(x)=L-1}", respectively. 
#' @param aVec An IntegerVector that represents a mapping table of multi-valued function.
#' @param k An integer that denotes the number of input variables.
#' @param L An integer that represents the level of multi-valued system (Default: 2).
#' @param Detail A logical value. Should return detail information of each edge?
#' (Default: \code{FALSE})
#' @param Redundancy A logical value. Should return the corresponding redundant 
#' values of all edges rather than the effective ones? (Default: \code{FALSE})
#' @details The results can be obtained through either direct or indirect 
#' methods, which are equivalent in outcome. These computational processes are 
#' transparent to R users and do not impact usability. Advanced developers can 
#' utilize prototype functions to conduct various analyses based on specific 
#' requirements. Relevant concepts and definitions can be found in paper 
#' [\href{https://doi.org/10.1073/pnas.2022598118}{Gates2021}].
#' @return A numeric value (or a NumericVector) denotes the effetive or redundant 
#' edge(s) of a multi-valued function.
#' @export
#' @examples
#' # Calculate the effetive input of a ternary function c(0,1,2, 0,1,2, 0,1,2).
#' # Its can be transformed to 0=f(*0), 1=f(*1), 2=f(*2)
#' # MulVFun_Effectiveness(c(0,1,2, 0,1,2, 0,1,2), k=2, L=3)
#' # Return 1.000 (Only lower bit can control mapping table)
#' MulVFun_Effectiveness(c(0,1,2, 0,1,2, 0,1,2), k=2, L=3)
#' 
#' # Calculate the effetive edges of a ternary function c(0,1,2, 0,1,2, 2,2,2).
#' # Its can be transformed to 0=f(00)|f(10), 1=f(01)|f(11), 2=f(2*)|f(*2)
#' # MulVFun_Effectiveness(c(0,1,2, 0,1,2, 2,2,2), k=2, L=3, Detail=TRUE)
#' # Return 0.72222 0.72222
#' # Explain the values: 0.72222
#' #  00 | 01 | 02 | 10 | 11 | 12 | 20 | 21 | 22 
#' #  00 | 01 | *2 | 10 | 11 | *2 |    |    | *2 
#' #     |    |    |    |    |    | 2* | 2* | 2* 
#' # ( 1 +  1 +  1 +  1 +  1 +  1 +  0 +  0 + 0.5 )/9 = 0.72222 # Lower one.
#' # (1  + 1  + 0  + 1  + 1  + 0  + 1  + 1  + 0.5 )/9 = 0.72222 # Higher one.
#' MulVFun_Effectiveness(c(0,1,2, 0,1,2, 2,2,2), k=2, L=3, Detail=TRUE)
#' 
MulVFun_Effectiveness<-function(aVec, k, L=2, Detail=FALSE, Redundancy=FALSE){
  a_vec=r_CheckValidMulVFun(aVec,k,L);
  if(Detail){# Show detail information
    xx=c_MulF_EffectiveEdges(a_vec,as.integer(k),as.integer(L));
  } else {# Return global information
    xx=c_MulF_Effective(a_vec,as.integer(k),as.integer(L));
  }
  if(Redundancy){# Get opposite value(s)
    if(Detail){
      xx=1.00-xx;
    }
    else {
      xx=k-xx;
    }
  }
  return(xx);
}