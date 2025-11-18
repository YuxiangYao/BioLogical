#' @title Analyze percolation within system
#' @description 
#' Here employs maximal stable components (MSC) of Boolean networks embedded in square 
#' lattice as measure of percolations. A random state of a given Boolean network
#' finally fall in an attractor with stable or oscillatory pattern. Here observe 
#' the MSC can form a continuous path as the occurrence of percolation. Please note
#' that this definition sources from [\href{10.1073/pnas.1534782100}{Shmulevich2003}]. 
#' There are also other definitions of so-called "percolation".
#' @param Size an integer, size of system.
#' @param SimStep an integer, steps of simulating system.
#' @param ObsWin an integer, the size of observing window to ensure the stable cluster.
#' @param OutPutState a logical value, should output all states? (Default: \code{FALSE})
#' @param OBF_Type a character, ordered Boolean type (OBF)
#' @param OBF_iPara1 an integer, configuration parameter for OBF.
#' @param OBF_iPara2 an integer, configuration parameter for OBF (Not necessary for some types of OBF).
#' @param OBF_Ratio a numeric value, proportion of ordered Boolean function within system.
#' @param RBF_Bias a NumvericVector, biases of random function (each is (0,1), sum=1).
#' @param Net_fPara a numeric value, topological configured parameters.
#' @param Init_1_Ratio a NumvericVector, proportion of each value (its length is \code{NumSys}).
#' @param NumSys an integer, number of discrete value.
#' @param LatType an integer, 4:square, 6: hexangular, 3:triangular
#' @param UpdateRule an integer, update rule: 1=synchronous, 2=asynchronous, 
#' 3=quick-asynchronous.
#' @details 
#' Ensure that all parameters are properly set. Some parameters are fixed due to specific tasks. 
#' \code{Size}, \code{SimStep}, \code{ObsWin} are dynamic parameters for simualtion. 
#' To avoid transient states, recommend \code{SimStep} >= \code{Size}. 
#' \code{OBF_Type}, \code{OBF_iPara1}, \code{OBF_iPara2}, \code{RBF_Bias}, 
#' \code{OBF_Ratio}, configure functions. See their roles in \link{BoolFun_Generator}. 
#' \code{Net_fPara} control lattice type: (\code{4}, Square; \code{3}, triangle;
#'  \code{6}, hexagon); here \code{Size} is an even squared number. 
#' @return A three-element list:
#' \itemize{
#'   \item \code{[[1]] NumericVector[2]}: [1] max stable cluster
#' fraction (non-zero means percolation happens); [2] stable node fraction.
#'   \item \code{[[2]] IntegerVector[Size]}: Stable or unstable states of 
#' all nodes. Value "1", "2", "3" represent the node being stable at state 0, 
#' stable at state 1, and unstable, respectively. (if \code{OutPutState} is 
#' \code{TRUE}). 
#' \item \code{[[3]] IntegerVector[Size]}: record if each stable nodes 
#' belong the maximal stable components (MSC).
#' }
#' @export
#' @examples
#' # Test percolation random and canalized functions (70%)
#' # set.seed(20250101L);
#' # DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#' #   OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.7);
#' # Return [[1]] [1] 0.6456 0.6628 (Percolation)
#' set.seed(20250101L)
#' DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#'   OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.7)
#' 
#' # Test percolation random and canalized functions (30%)
#' # set.seed(20250102L);
#' # DNS_Percolation(2500L, 2500L, 2000L, OBF_Type='C', 
#' #  OBF_iPara1=2, OBF_iPara2=-1, OBF_Ratio=0.3);
#' # Return [[1]] [1] 0.0000 0.3448 (No percolation)
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
