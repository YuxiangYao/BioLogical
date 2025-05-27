#' @title Analyze damage spread within system
#' @description Here employs the concept of Derrida Curve to measure the damage's 
#' effect. Its main process involves: randomly generate two system states 
#' with a specified Hamming distance; two states spontaneously evolve based on 
#' same update rules; observe final their Hamming distance after certain steps. 
#' Detail concepts can been found in 
#' [\href{https://doi.org/10.1016/j.physa.2004.05.018}{Kauffman,2004}]. 
#' @param Size integer, size of system
#' @param SimStep integer, steps of simulating system
#' @param Init_Dist float, initial normalized Hamming distance of two vector
#' @param Init_1_Ratio floatvector, proportion of each value (same length as L)
#' @param OBF_Type a character, ordered Boolean function type (OBF)
#' @param OBF_iPara1 integer, configuration parameter for OBF
#' @param OBF_iPara2 integer, configuration parameter for OBF (Not necessary for some types of OBF)
#' @param OBF_Ratio float, proportion of ordered Boolean function within system
#' @param RBF_Bias floatvector, biases of random function (each is (0,1), sum=1)
#' @param Net_Type a character, system topological type
#' @param Net_fPara float, topological configured parameters
#' @param NumSys integer, number of discrete value
#' @param UpdateRule integer, update rule: 1 syn, 2 asy, 3 quick-asy
#' @details
#' Ensure that all parameters are properly set. \code{Size}, \code{SimStep}, \code{Init_Dist}, \code{Init_1_Ratio} are dynamic parameters for simualtion. 
#' To avoid transient states, recommend \code{SimStep} >= \code{Size}. \code{Init_Dist} has little effect on the final distance, 
#' but \code{Init_1_Ratio} can casue difference in some cases (when OBF is D-type). \code{OBF_Type}, \code{OBF_iPara1}, \code{OBF_iPara2}, \code{RBF_Bias}, \code{OBF_Ratio}, 
#' configure ordered functions. See their detail roles in \link{BoolFun_Generator}. \code{Net_Type} and \code{Net_fPara} determine the topological type.
#'  \itemize{
#'   \item \code{K}: Kauffman model, \code{Net_fPara} is in-degree.
#'   \item \code{E}: Erdos-Renyi graph, \code{Net_fPara} is average degree.
#'   \item \code{R}: Regular random graph, \code{Net_fPara} is connecting number.
#'   \item \code{L}: Lattic sqaure, \code{Net_fPara} is type (4, Square; 3, triangle; 6, hexagon).
#'   \item \code{N}: Null model (Not enabled here).
#' }
#' @return float, normalized Hamming distance of two vector at final time
#' @export
#' @examples
#' # set.seed(20250101L)
#' # DNS_DamageSpread(Size=1000L, SimStep=1000L, OBF_Type='C',
#' #   OBF_iPara1=1, OBF_iPara2=-1, OBF_Ratio=0.5)
#' # Return 0.335
#' set.seed(20250101L)
#' DNS_DamageSpread(Size=1000L, SimStep=1000L, OBF_Type='C',
#'   OBF_iPara1=1, OBF_iPara2=-1, OBF_Ratio=0.5)
DNS_DamageSpread<-function(Size=1000L, SimStep=1000L, Init_Dist=0.1, Init_1_Ratio=c(0.5,0.5),
  OBF_Type='R', OBF_iPara1=1L, OBF_iPara2=1L, OBF_Ratio=0.1, RBF_Bias=c(0.5,0.5),
  Net_Type='K', Net_fPara=4.00,
  NumSys=2L, UpdateRule=1L){
  # Check parameter.
  f_par=c(Size, SimStep, Init_Dist, Init_1_Ratio, OBF_iPara1, OBF_iPara2, OBF_Ratio, RBF_Bias, Net_fPara);
  c_par=c(OBF_Type,Net_Type);
  if((!all(is.numeric(f_par)))||(!all(is.character(c_par)))){
    stop("Invalid inputs. Please check the help documentation.\n");}
  if((0>Size)||(0>SimStep)||(0>Init_Dist||Init_Dist>1)||
    any(0>Init_1_Ratio)||any(Init_1_Ratio>1)||
    (0>OBF_Ratio||OBF_Ratio>1)||
    any(RBF_Bias<0)||any(RBF_Bias>1)||(NumSys>2&&length(RBF_Bias)!=NumSys)||
    (0>Net_fPara||Net_fPara>12)||# OBF_iPara1, OBF_iPara2 checked in other fun.
    !(OBF_Type%in%c('R','C','P','M','D','T'))||
    !(Net_Type%in%c('K','E','R','L'))){
    stop("Invalid inputs. Please check the help documentation.\n");}
  # Simulate
  xx=c_Derrida_Simualtion(as.integer(Size), as.integer(NumSys),
    as.integer(SimStep), OBF_Type[1], RBF_Bias, OBF_Ratio, 
    Init_Dist, Init_1_Ratio, Net_Type[1], Net_fPara, 
    as.integer(OBF_iPara1), as.integer(OBF_iPara2), as.integer(UpdateRule));
  return (xx);
}
