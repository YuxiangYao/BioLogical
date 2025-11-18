#' @title Calculate the input/edge effetiveness of a Boolean function
#' @description The function employs the Quine-McCluskey method to derive the 
#' prime implicants of input vectors for which \eqn{\{\vec{x}|f(\vec{x})=0\}}{\{x|f(x)=0\}}
#'  and  \eqn{\{\vec{x}|f(\vec{x})=1\}}{\{x|f(x)=1\}},
#' respectively. 
#' @param Bit_Vec A Boolean vector (or 0/1-valued \code{IntegerVector}) represents 
#' the truth table of a Boolean function.
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
#' edge(s) of a Boolean function. 
#' @export
#' @examples
#' # Calculate the effetive input of a Boolean function c(0,0,0,0,0,0,1,1).
#' # Its can be transformed as: 0=f(*0*), 0=f(0**), and 1=f(11*).
#' # BoolFun_Effectiveness(c(0,0,0,0,0,0,1,1));
#' # Return 1.25000
#' BoolFun_Effectiveness(c(0,0,0,0,0,0,1,1))
#' 
#' # Calculate the effetive edges of a Boolean function c(0,0,0,1,1,1,1,1).
#' # Its can be transformed to 0=f(00*), 0=f(0*0), 1=f(1**) or f(*11)
#' # BoolFun_Effectiveness(c(0,0,0,1,1,1,1,1), Detail=TRUE);
#' # Return 0.3750 0.3750 0.8125 
#' # Explain lowest/middle bit's 0.3750 and highest bit's 0.8125:
#' #  000 | 001 | 010 | 011 | 100 | 101 | 110 | 111  <== 2^3 vectors can be  
#' #  00* | 00* |     |     | 1** | 1** | 1** | 1**      covered by which 
#' #  0*0 |     | 0*0 | *11 |     |     |     | *11      prime implicants?
#' # (0.5 +   0 +   1 +   1 +   0 +   0 +   0 + 0.5 )/8 ==> 0.3750  # Lowest
#' # (0.5 +  1  +  0  +  1  +  0  +  0  +  0  + 0.5 )/8 ==> 0.3750  # Lowest
#' # ( 1  + 1   + 1   + 0   + 1   + 1   + 1   + 0.5 )/8 ==> 0.8125  # Highest bit
#' 
#' BoolFun_Effectiveness(c(0,0,0,1,1,1,1,1), Detail=TRUE)
#' 
BoolFun_Effectiveness<-function(Bit_Vec, Detail=FALSE, Redundancy=FALSE){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  if(Detail){# Show detail information
    xx=c_BF_EffectiveEdges(tmp$bits,tmp$kin);
  } else {# Return global information
    xx=c_BF_Effective(tmp$bits,tmp$kin);
  }
  if(Redundancy){# Get opposite value(s)
    if(Detail){
      xx=1.00-xx;
    }
    else {
      xx=tmp$kin-xx;
    }
  }
  return(xx);
}