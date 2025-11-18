#' @title Generate a specific type of multi-valued function
#' @description This function can generate following types: Canalization, 
#' Linear threshold-based, Dominated-valued. 
#' @param MF_Type A character included in the following:
#'  \itemize{
#'   \item \code{C}: Canalization
#'   \item \code{D}: Dominated-valued
#'   \item \code{T}: Linear threshold-based
#'   \item \code{R}: Random
#' }
#' @param k An integer that denotes the number of input variables (Default: 3).
#' @param L An integer that represents the level of multi-valued system (Default: 3).
#' @param Bias A NumericVector representing the probabilities of elements in 
#' non-controlled slots within mapping tables (See \code{Details}).
#' By default, \code{NULL}, which corresponds to the vector \code{rep(1.0/L, L)}.
#' If provided, the vector must have a length greater than L-1 and contain only
#' non-negative real numbers. The normalization of its first L elements determines
#' the baseline configuration.
#' @param CanaDeep An integer ranging from 1 to k, indicating k-layer canalization
#' for type \code{'C'}. (Default: 1)
#' @param CanaVar A \code{CanaDeep}-length non-repeating IntegerVector, indicating
#' which variables should be canalized, with values ranging from 1 to k.
#' @param CanaVarNum An integer vector of length \code{CanaDeep} specifying the 
#' number of each canalizing variable (ranging from 1 to L). If \code{CanaLayerInfo}
#' is appropriately provided, \code{CanaVarNum} is ignored. (Default: \code{NULL}, 
#' automatically set using \code{sample.int(L, CanaDeep, replace=TRUE)} implemented
#' in C++, rather than the native R function).
#' @param CanaLayerInfo A list containing detail information for canalization:
#' [[1]] for input (canalizing patterns); [[2]] for output (canalized patterns).
#' This argument must meet specific requirements (see \code{Details}). 
#' If not provided (Default: \code{NULL}), the function will automatically 
#' generate a random configuration. For detailed meaning, refer to \link{BoolFun_Type}.
#' @param MappingTable A logical value. Is the mapping table also returned? 
#' (Default: \code{FALSE}) 
#' @details The parameter Bias applies to non-controlled slots in the mapping table. 
#' For instance, a canalized function (k=3 and L=3) where \eqn{f(x_1=0)=2}{f(x1=0)=2}, indicates 
#' that only inputs with \eqn{x_1\not=0}{x1!=0} can be assigned probabilistic 
#' values according to the specified bias. 
#' 
#' The \code{CanaLayerInfo} must satisfy the following requirements: 
#' it must contain at least two sub-lists, representing the canalizing and 
#' canalized configurations, respectively. Each sub-list must have a length 
#' exceeding \code{k}. Within each sub-list, all canalizing and canalized values
#' must be strictly less than \code{L}. Furthermore, the elements in 
#' \code{CanaLayerInfo[[1]][[i]]} must be unique, which represent different 
#' canalizing values (See examples of \link{MulVFun_is_NestedCana}). 
#' Any violation of these conditions will lead to an error message and execution
#' halts.
#' @return An integervector of length \eqn{L^k}{L^k}.
#' @export
#' @examples
#' # Generate 3-input ternary function. Please note the canalizing/canalized 
#' # values in the following scenarios:
#' set.seed(1234L)
#' Example_01=MulVFun_Generator('C', 3L, 3L, MappingTable=TRUE)
#' Example_02=MulVFun_Generator('T', 3L, 3L, MappingTable=TRUE)
#' Example_03=MulVFun_Generator('D', 3L, 3L, MappingTable=TRUE)
#' # 1=Yes, -1=No (for Canalized)
#' # 1=Yes, 0=Unknown, -1=No (for Threshold and Domainted)
#' 
#' MulVFun_is_NestedCana(Example_01[,4],3L,3L)[[1]]
#' # 1
#' MulVFun_is_Threshold(Example_01[,4],3L,3L)[[1]]
#' # -1
#' MulVFun_is_Domainted(Example_01[,4],3L,3L)[[1]]
#' # -1
#' MulVFun_is_NestedCana(Example_02[,4],3L,3L)[[1]]
#' # -1
#' MulVFun_is_Threshold(Example_02[,4],3L,3L)[[1]]
#' # 1
#' MulVFun_is_Domainted(Example_02[,4],3L,3L)[[1]]
#' # -1
#' MulVFun_is_NestedCana(Example_03[,4],3L,3L)[[1]]
#' # -1
#' MulVFun_is_Threshold(Example_03[,4],3L,3L)[[1]]
#' # -1
#' MulVFun_is_Domainted(Example_03[,4],3L,3L)[[1]]
#' # 1
#' 
MulVFun_Generator<-function(MF_Type=c('R','C','D','T'), k=3L, L=3L, Bias=NULL, 
  CanaDeep=1L, CanaVar=NULL, CanaVarNum=NULL, CanaLayerInfo=NULL,
  MappingTable=FALSE){
  k_i=as.integer(k);
  L_i=as.integer(L);
  if(abs(k_i-k)>1e-5||abs(L_i-L)>1e-5||k_i<1||L_i>10||L_i^k_i>5000000){
    stop("Please enter integers for 'k' and 'L'. The length of (L^k) should be smaller than 5,000,000.\n");}
  if(is.null(Bias[1])||is.na(Bias[1])){
    bias=rep(1.0/as.numeric(L_i), L_i);
  } else {
    if(length(Bias)<L_i){# Check bias.
      stop("Is following condtion 'length(Bias)>=L' met?");
    } else {
      if((is.numeric(Bias)||is.integer(Bias))&&(!any(is.na(Bias)))&&all(Bias>=0)){
        bias=Bias[1:L];
        if(sum(bias)<1e-7){
          warning("Summaiton of Bias[1:L] is zero. Here reset Bias as rep(1.0/L, L)",
            call. = FALSE, immediate.=TRUE);
          bias=rep(1.0/as.numeric(L_i), L_i);
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
    if(CanaDeep>k_i||CanaDeep<1){
      stop("Please check: '1<=CanaDeep<=k'.");
    }
    # Check CanaVar
    if(is.null(CanaVar)){
      cana.var=as.integer(c(-1L,1L));# "-1" means in CPP start ID is "0".
    } else {
      if(length(CanaVar)!=CanaDeep||!(all(CanaVar>=0)&&all(CanaVar<k_i))||any(duplicated(CanaVar))){
        stop("Please check 'CanaVar' is appropriate: non-repeating, [0,k-1], length(CanaVar)==CanaDeep!!");
      }
      else {
        cana.var=as.integer(CanaVar);# R id to CPP id.
      }
    }
    # Check CanaVarNum
    if(is.null(CanaVarNum)){
      cana.var.num=as.integer(c(-1L,1L));
    } else {
      if(length(CanaVarNum)!=CanaDeep||!(all(CanaVarNum>0)&&all(CanaVarNum<=L_i))){
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
        num1=sapply(CanaLayerInfo[[1]], function(x) (all(x<L)&&all(x>=0)));# Each type should be [0,L-1]
        num2=sapply(CanaLayerInfo[[2]], function(x) (all(x<L)&&all(x>=0)));# Each type should be [0,L-1]
        dup1=sapply(CanaLayerInfo[[1]], function(x) (any(duplicated(x))));
        tt1=all(len1>0)&&all(len1<=L_i)&&all(len2==len1)&&all(num1)&&all(num2)&&(!any(dup1));
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

  resvec=c_MulVF_Generator(MF_Type[1],k_i,L_i,
    CanaDeep, cana.var, cana.var.num, c.layer.info.1, c.layer.info.2, bias, cana_free);
  if(MappingTable){
    resvec=cbind(r_GenerateTruthTable_M(k_i,L_i), resvec);
    colnames(resvec)[k+1]="f_out";
  }
  return (resvec);
}
