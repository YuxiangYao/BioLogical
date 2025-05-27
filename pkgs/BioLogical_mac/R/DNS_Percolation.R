#' @title Analyze percolation within system
#' @description 
#' Here employs maximal stable components (MSC) of Boolean networks embedded in square 
#' lattice as measure of percolations. A random state of a given Boolean network
#' finally fall in one attractor with stable or oscillatory nodes. Here observe 
#' the MSC can form a continuous path as the occurrence of percolation. Please note
#' that this definition sources from [\href{10.1073/pnas.1534782100}{Shmulevich,2003}]. 
#' There are also other definitions of "percolation".
#' @param Size integer, size of system
#' @param SimStep integer, steps of simulating system
#' @param ObsWin integer, Observe window, ensure the stable cluster.
#' @param OutPutState bool, should output all states? (Default: \code{FALSE})
#' @param OBF_Type a character, ordered Boolean type (OBF)
#' @param OBF_iPara1 integer, configuration parameter for OBF
#' @param OBF_iPara2 integer, configuration parameter for OBF (Not necessariy for some types of OBF)
#' @param OBF_Ratio float, proportion of ordered Boolean function within system
#' @param RBF_Bias floatvector, biases of random function (each is (0,1), sum=1)
#' @param Net_fPara float, topological configured parameters
#' @param Init_1_Ratio floatvector, proportion of each value (same length as L)
#' @param NumSys integer, number of discrete value
#' @param LatType integer, 4:square, 6: hexangular, 3:triangular
#' @param UpdateRule integer, update rule: 1 syn, 2 asy, 3 quick-asy
#' @details 
#' Ensure that all parameters are properly set. Some parameters are fixed due to specific tasks. \code{Size}, \code{SimStep}, 
#' \code{ObsWin} are dynamic parameters for simualtion. To avoid transient states, recommend \code{SimStep} >= \code{Size}. 
#' \code{OBF_Type}, \code{OBF_iPara1}, \code{OBF_iPara2}, \code{RBF_Bias}, \code{OBF_Ratio}, configure ordered functions. 
#' See their detail roles in \link{BoolFun_Generator}. \code{Net_fPara} control lattice type: (\code{4}, Square; \code{3}, triangle; \code{6}, hexagon); here \code{Size} is an even squared number. 
#' @return List, [[1]] NumericVector[2], [1] max stable cluster fraction, 
#' non-zero means percolation happens; [[2]] IntegerVector[\code{Size}], 
#' System stable/unstable state cases (if \code{OutPutState} is \code{TRUE}). 
#' [[3]] if those stable nodes belong to MAX cluster.
#' @export
#' @examples
#' # Test percolation random and canalized functions (70%)
#' # set.seed(20250101L)
#' # DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#' #   OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.7)
#' # Return [[1]] [1] 0.6456 0.6628
#' set.seed(20250101L)
#' DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#'   OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.7)
#' 
#' # Test percolation random and canalized functions (30%)
#' # set.seed(20250102L)
#' # DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#' #  OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.3)
#' # Return [[1]] [1] 0.0000 0.3448
#' set.seed(20250102L)
#' DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#'   OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.3)
#' 
DNS_Percolation<-function(Size=1000L, SimStep=1000L, ObsWin=1000L, OutPutState=FALSE,
  OBF_Type='R', OBF_iPara1=1L, OBF_iPara2=1L, OBF_Ratio=0.1, RBF_Bias=c(0.5,0.5),
  Net_fPara=4.00, Init_1_Ratio=c(0.5,0.5), NumSys=2L, LatType=4L, UpdateRule=1L){
  # Check parameter.
  f_par=c(Size, SimStep, ObsWin, OBF_iPara1, OBF_iPara2, OBF_Ratio, RBF_Bias, Net_fPara);
  c_par=c(OBF_Type);
  if((!all(is.numeric(f_par)))||(!all(is.character(c_par)))){
    stop("Invalid inputs. Please check the help documentation.\n");}
  if((0>Size)||(0>SimStep)||(0.1*Size>ObsWin)||
    any(0>Init_1_Ratio)||any(Init_1_Ratio>1)||
    any(RBF_Bias<0)||any(RBF_Bias>1)||(NumSys>2&&length(RBF_Bias)!=NumSys)||
    (0>Net_fPara||Net_fPara>12)||# OBF_iPara1, OBF_iPara2 checked in other fun.
    !(OBF_Type%in%c('R','C','P','M','D','T'))||
    !(LatType%in%c(3,4,6)) ){
      stop("Invalid inputs. Please check the help documentation.\n");}
  # Simulate 
  xx=c_Percolation_Simualtion(
    as.integer(Size), as.integer(NumSys), as.integer(SimStep), 
    as.integer(LatType), as.integer(ObsWin), 
    OBF_Type[1], RBF_Bias, OBF_Ratio,
    Init_1_Ratio, 'L',# Fixed parameters.
    Net_fPara, as.integer(OBF_iPara1), as.integer(OBF_iPara2), 
    OutPutState, as.integer(UpdateRule));
  if(!OutPutState){# null
    xx[[2]]=NULL;
    xx[[3]]=NULL;}
  xx[[1]]=xx[[1]]/as.numeric(Size);# Convert to pro
  return (xx);
}
