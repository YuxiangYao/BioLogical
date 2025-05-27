#' @title Calculate the effetive edges loading of a Boolean function
#' @description The analysis depends on Quine-McCluskey method to obtain the 
#' prime implicants of input vectors belong to "f(x)=1" and "f(x)=0", respectively. 
#' The results within a function can be obtained through direct and indirect 
#' methods that are equivalent. They are invisible to users and do not affect usage.
#' Relevant concepts and definitions are found in the paper
#' [\href{https://doi.org/10.1073/pnas.2022598118}{Gates,2021}].
#' @param Bit_Vec a bool (or 0/1-integral) vector, represents a truth table of Boolean function.
#' @param Detail logical, should return detail of each edge/input information?
#' @param Redundancy logical, should return corresponding redundant values instead of effective values?
#' @return double, Effetive edges of a Boolean function
#' @export
#' @examples
#' # Calculate the effetive edges of a Boolean function c(0,0,0,0,0,0,1,1).
#' # Its can be transformed to 0=f(*0*), 0=f(0**), 1=f(11*)
#' # BoolFun_EffectiveEdges(c(0,0,0,0,0,0,1,1))
#' # Return 1.25000
#' BoolFun_EffectiveEdges(c(0,0,0,0,0,0,1,1))
#' 
#' # Calculate the effetive edges of a Boolean function c(0,0,0,1,1,1,1,1).
#' # Its can be transformed to 0=f(00*), 0=f(0*0), 1=f(1**) or f(*11)
#' # BoolFun_EffectiveEdges(c(0,0,0,1,1,1,1,1), Detail=TRUE)
#' # Return 0.3750 0.3750 0.8125 
#' # Explain lowest/middle bit's 0.3750 and highest bit's 0.8125:
#' #  000 | 001 | 010 | 011 | 100 | 101 | 110 | 111  
#' #  00* | 00* | 0*0 | *11 | 1** | 1** | 1** | 1** 
#' #  0*0 |     |     |     |     |     |     | *11 
#' # (0.5 +   0 +   1 +   1 +   0 +   0 +   0 + 0.5 )/8 ==> 0.3750  # Lowest
#' # ( 1  + 1   + 1   + 0   + 1   + 1   + 1   + 0.5 )/8 ==> 0.8125  # Highest
#' # More detail can see the paper [\href{https://doi.org/10.1073/pnas.2022598118}{Gates,2021}].
#' BoolFun_EffectiveEdges(c(0,0,0,1,1,1,1,1), Detail=TRUE)
BoolFun_EffectiveEdges<-function(Bit_Vec, Detail=FALSE, Redundancy=FALSE){
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