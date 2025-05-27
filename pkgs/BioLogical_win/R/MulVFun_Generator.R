#' @title Generate a special type multi-valued function
#' @description This function can generate following types: Canalization, 
#' Linear threshold-based, Dominated-valued. 
#' @param MF_Type a character contained in the following:
#'  \itemize{
#'   \item \code{C}: Canalization
#'   \item \code{D}: Dominated-valued
#'   \item \code{T}: Linear threshold-based
#'   \item \code{R}: Random
#' }
#' @param k An integer for the number of input-variable.
#' @param L An integer for the level of discrete system.
#' @param Bias a numbericVector. Probabilities of elements in non-controlled slots in mapping tables (see \code{Details}).
#' Default \code{NULL}, means the vector: rep(1.0/L, L).
#' Its length is larger than L-1. All elements are non-minus real numbers
#' Its first L element's normalization serves as the bais's configuration. 
#' @param CanaDeep an integer of 1~k, denoting k-layer canalization for type \code{'C'}. (Default: 1)
#' @param CanaVar an \code{CanaDeep}-length non-repeating integer vector, denoting WHICH variabls should be canalized, [1~k] 
#' @param CanaVarNum an \code{CanaDeep}-length integer vector, denote numbers of each canalizing variable [1~L].
#' If \code{CanaLayerInfo} provided suitably, \code{CanaVarNum} is useless (Default: NULL, 
#' automatically set as sample.int(L, CanaDeep, replace=TRUE) implemented in C++ not rather than the corresponding 
#' native R functions).
#' @param CanaLayerInfo a list, detail information for configuration of canalization. [[1]] for input; [[2]] for output.
#' This argument can should meet specific requirements (see \code{Details}). If not provided (Default: \code{NULL}), 
#' function will randomly configurate them.
#' Detail meaning see References of \link{BoolFun_Type}
#' @param MappingTable logical value. Is the mapping table also returned? (Default: \code{FALSE}) 
#' @details \code{Bias} merely acts only non-controled slots in mapping table. For example, for a funciton that 
#' k=3,L=3, f(x1=0) --> f(~)=2 means for all x1!=0 can be filled with probabilities (bias).
#' Rrequirements \code{CanaLayerInfo}: least two sub-lists denote canalizing and canalized configurating information.
#' The length of each sub-list should be more than \code{k}. In the each sub-list, canalizing and canalized 
#' values should be smaller than \code{L}, and \code{CanaLayerInfo[[1]][[i]]} can't have duplicate elements in vector X
#' Any condition not met wolud throw an error information.
#' One following example explain this concept.
#' @return L^k length integervector
#' @export
#' @examples
#' 
#' # Generate 3-input ternary function. Please note the canalizing/canalized 
#' # values in the following scenarios:
#' set.seed(1234L)
#' Example_01=MulVFun_Generator('C', 3L, 3L, MappingTable=TRUE)
#' Example_02=MulVFun_Generator('T', 3L, 3L, MappingTable=TRUE)
#' Example_03=MulVFun_Generator('D', 3L, 3L, MappingTable=TRUE)
#' # 1=Yes, -1=No (for Canalized)
#' # 1=Yes, 0=Unknown, -1=No (for Threshold and Domainted)
#' MulVFun_is_NestedCanalized(Example_01[,4],3L,3L)[[1]]# 1
#' MulVFun_is_Threshold(Example_01[,4],3L,3L)[[1]]# -1
#' MulVFun_is_Domainted(Example_01[,4],3L,3L)[[1]]# -1
#' MulVFun_is_NestedCanalized(Example_02[,4],3L,3L)[[1]]# -1
#' MulVFun_is_Threshold(Example_02[,4],3L,3L)[[1]]# 1
#' MulVFun_is_Domainted(Example_02[,4],3L,3L)[[1]]# -1
#' MulVFun_is_NestedCanalized(Example_03[,4],3L,3L)[[1]]# -1
#' MulVFun_is_Threshold(Example_03[,4],3L,3L)[[1]]# -1
#' MulVFun_is_Domainted(Example_03[,4],3L,3L)[[1]]# 1
#' 
MulVFun_Generator<-function(MF_Type=c('R','C','D','T'), k=3L, L=3L, Bias=NULL, 
  CanaDeep=1L, CanaVar=NULL, CanaVarNum=NULL, CanaLayerInfo=NULL,
  MappingTable=FALSE){
  kk=as.integer(k); LL=as.integer(L);
  if(kk<1||LL>10||L^k>10000){# !is.integer(k)||!is.integer(L)||
    stop("Please enter integers for 'k' and 'L'. The length of (L^k) should be smaller than 10,000.\n");}
  if(is.null(Bias[1])||is.na(Bias[1])){
    bias=rep(1.0/as.numeric(L), LL);
  } else {
    if(length(Bias)<LL){# Check bias.
      stop("Is following condtion 'length(Bias)>=L' met?");
    } else {
      if((is.numeric(Bias)||is.integer(Bias))&&(!any(is.na(Bias)))&&all(Bias>=0)){
        bias=Bias[1:LL];
        if(sum(bias)<1e-7){
          warning("Summaiton of Bias[1:L] is zero. Here reset Bias as rep(1.0/L, L)",
            call. = FALSE, immediate.=TRUE);
          bias=rep(1.0/as.numeric(L), L);
        } else {
          bias=bias/sum(bias);
        }
      } else {
        stop("'Bias' should be a numeric/intger non-minus vector!");
      }
    }
  }
  # Configurate parameters.
  cana.var=as.integer(c(NA,1L));
  cana.var.num=as.integer(c(NA,1L));
  c.layer.info.1=c.layer.info.2=list(a=NA,b=NA);
  cana_free=TRUE;
  if('C'==MF_Type[1]){
    # Check CanaDeep
    if(CanaDeep>k||CanaDeep<1){
      stop("Please check: '1<=CanaDeep<=k'.");
    }
    # Check CanaVar
    if(is.null(CanaVar)){
      cana.var=as.integer(c(-1L,1L));# "-1" means in CPP start ID is "0".
    } else {
      if(length(CanaVar)!=CanaDeep||!(all(CanaVar>=0)&&all(CanaVar<k))||any(duplicated(CanaVar))){
        stop("Please check 'CanaVar' is appropriate: non-repeating, [0,k-1], length(CanaVar)==CanaDeep!!");
      }
      else {
        cana.var=as.integer(CanaVar)-1L;# R id to CPP id.
      }
    }
    # Check CanaVarNum
    if(is.null(CanaVarNum)){
      cana.var.num=as.integer(c(-1L,1L));
    } else {
      if(length(CanaVarNum)!=CanaDeep||!(all(CanaVarNum>0)&&all(CanaVarNum<=LL))){
        stop("Please check 'CanaVarNum' is appropriate [1,L], length(CanaVarNum)==CanaDeep!!");
      }
      else {
        cana.var.num=as.integer(CanaVarNum);
      }
    }
    # CanaLayerInfo
    if(is.null(CanaLayerInfo)){# Not provide by user, will randomly generate via function.
        c.layer.info.1=list(a=as.integer(NA),b=0L);
    } else {
      if(length(CanaLayerInfo)<2||length(CanaLayerInfo[[1]])<CanaDeep||length(CanaLayerInfo[[2]])<CanaDeep){
        stop("'CanaLayerInfo' not meet the condition. Please check help document.");
      } else {
        len1=sapply(CanaLayerInfo[[1]], function(x) (length(x)));# Cana_In: variable's canalizing type number. 
        len2=sapply(CanaLayerInfo[[2]], function(x) (length(x)));# Cana_Out: variable's canalizing type number.
        num1=sapply(CanaLayerInfo[[1]], function(x) (all(x<LL)&&all(x>=0)));# Each type should be [0,L-1]
        num2=sapply(CanaLayerInfo[[2]], function(x) (all(x<LL)&&all(x>=0)));# Each type should be [0,L-1]
        dup1=sapply(CanaLayerInfo[[1]], function(x) (any(duplicated(x))));
        tt1=all(len1>0)&&all(len1<=LL)&&all(len2==len1)&&all(num1)&&all(num2)&&(!any(dup1));
        if(!tt1){
          stop("'CanaLayerInfo' not meet the condition. Please check help document.");
        } else {
          c.layer.info.1=CanaLayerInfo[[1]];
          c.layer.info.2=CanaLayerInfo[[2]];
          cana_free=FALSE;
        }
      }
    }
  } else if(MF_Type[1]%in%c('R','D','T')){
      ;# Do nothing
  } else {
    stop("Illegal function type. Please check help document.\n");}

  resvec=c_MulVF_Generator(MF_Type[1],kk,LL,
    CanaDeep, cana.var, cana.var.num, c.layer.info.1, c.layer.info.2, bias, cana_free);
  if(MappingTable){
    resvec=cbind(r_GenerateTruthTable_M(kk,LL), resvec);
    colnames(resvec)[kk+1]="f_out";
  }
  return (resvec);
}
