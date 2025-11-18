#' @title Simulating a transient process (Automatically generate network)
#' @description Same as \link{DNS_DamageSpread}, but record the transient process.
#' @param Size an integer, size of system.
#' @param SimStep an integer, steps of simulating system.
#' @param Init_Dist a numeric value, initial normalized Hamming distance of two vector.
#' @param Init_1_Ratio a NumvericVector, proportion of each value (its length is \code{NumSys}).
#' @param OBF_Type a character, ordered Boolean function type (OBF).
#' @param OBF_iPara1 an integer, configuration parameter for OBF.
#' @param OBF_iPara2 an integer, configuration parameter for OBF (Not necessary for some types of OBF).
#' @param OBF_Ratio a numeric value, proportion of ordered Boolean function within system.
#' @param RBF_Bias a NumvericVector, biases of random function (each is (0,1), sum=1).
#' @param Net_Type a character, system topological type.
#' @param Net_fPara a numeric value, topological configured parameters.
#' @param NumSys an integer, number of discrete value.
#' @param UpdateRule an integer, update rule: 1=synchronous, 2=asynchronous, 
#' 3=quick-asynchronous.
#' @details See \link{DNS_DamageSpread}.
#' @return A list containing two elements: 
#' \itemize{
#' \item \code{[[1]] Transient states}: A matrix whose rows and columns represent 
#'  time points and nodes, respectively.
#' \item \code{[[2]] Network}: The structure is consistent with that in 
#' \link{BoolGRN_CellCollective}
#' }
#' @export
#' 
Transient_Process_G<-function(Size=1000L, SimStep=1000L, Init_Dist=0.1, Init_1_Ratio=c(0.5,0.5),
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
  xx=cin_OnlyTempState(as.integer(Size), as.integer(NumSys),
    as.integer(SimStep), OBF_Type[1], RBF_Bias, OBF_Ratio, 
    Init_Dist, Init_1_Ratio, Net_Type[1], Net_fPara, 
    as.integer(OBF_iPara1), as.integer(OBF_iPara2), as.integer(UpdateRule));
  xx[[1]]=matrix(xx[[1]],as.integer(SimStep),as.integer(Size),byrow=TRUE);
  return (xx);
}
