#' @title Analyze potential engaged nodes/scaling law within Boolean/multi-valued system
#' @description This analysis only emphasizes the static feature of the network.
#' The function recursively analyzes engaged nodes that contribute to systematic 
#' final dynamic feature. The not "engaged" nodes including terminal nodes 
#' (out-degree is 0) and useless nodes (all out-links are useless). "Useless":
#' {A,B,C} ==> X, and f(X)=A*B, so C is useless for X. If all out-edges/links of 
#' C meets this scenario, node C belongs to useless class. Detail information and 
#' concepts can see the paper [\href{https://doi.org/10.1103/PhysRevLett.90.068702}{Socolar,2003}].
#' @param Size integer, size of system
#' @param OBF_Type a character, ordered Boolean type (OBF)
#' @param OBF_iPara1 integer, configuration parameter for OBF
#' @param OBF_iPara2 integer, configuration parameter for OBF (Not necessariy for some types of OBF)
#' @param OBF_Ratio float, proportion of ordered Boolean function within system
#' @param RBF_Bias floatvector, biases of random function (each is (0,1), sum=1)
#' @param Net_Type a character, system topological type
#' @param Net_fPara float, topological configured parameters
#' @param Controller integerVector, denote which nodes should be manually controlled. (Default: \code{NULL})
#' @param ControllerValue Boolean vector, controlled nodes' values. (Default: \code{NULL})
#' @param NodeAttri logic, should return node's attribute? (Default: \code{FALSE})
#' @param ResidualNet logic, should return residual network? (Default: \code{FALSE})
#' @param NumSys integer, number of discrete value
#' @details
#' Ensure that all parameters are properly set. \code{Size} is dynamic parameters for simualtion. \code{OBF_Type}, \code{OBF_iPara1}, \code{OBF_iPara2}, \code{RBF_Bias}, \code{OBF_Ratio}, configure ordered functions. See their detail roles in \link{BoolFun_Generator}. \code{Net_Type} and \code{Net_fPara} determine the topological type.
#'  \itemize{
#'   \item \code{K}: Kauffman model, \code{Net_fPara} is in-degree.
#'   \item \code{E}: Erdos-Renyi graph, \code{Net_fPara} is average degree.
#'   \item \code{R}: Regular random graph, \code{Net_fPara} is connecting number.
#'   \item \code{L}: Lattic sqaure, \code{Net_fPara} is type (\code{4}, Square; \code{3}, triangle; \code{6}, hexagon).
#'   \item \code{N}: Null model (Not enabled here).
#' } 
#' 
#' \code{Controller} denotes manual setting of controlling nodes. Node's indexes or names are both acceptable. Note
#' that function does not check index's correctness. Please ensure the correct indexes. \code{ControllerValue} denotes 
#' corresponding values of controlled genes. If not provided (\code{NULL}), they are generated randomly. If provided,
#' only the [1,N_Control] are utilized; shorter than N_Control return an error message and stop. One gene is
#' both included in \code{Controller} and \code{Exponents}, the configurated value of the former covers the latter one. 
#' 
#' \code{Residual networks} still contain all nodes within original systems for comparison. Stable nodes, useless 
#' nodes, and corresponding egdes would be removed and cut-off. Remaing nodes, edges, and Boolan functions
#' are engaged in terminal dynamic behaviors.
#' 
#' @return List[[5]]: [[1]] overall information (See example) 
#' [[2]] stable and unstable node's detail information each possible state denoted as (0101...)2=> a int32
#' [[3]] nodes are stable, useless, engaged ?
#' [[4]] Possible Value information: a brief vector denotes stable(1), two-possible(2), ....
#' [[5]] the residual network structure
#' @export
#' @examples
#' # Test a (k=2, p=0.5) Kauffman model. It is critical 2p(1-p)K=1.
#' # set.seed(20250101L);
#' # DNS_Engaged(10000L, RBF_Bias=0.5, Net_Type='K', Net_fPara=2.00);
#' # Return $[[Overview]] (other are NA)
#' #  Stable    Useless    Engaged   External Controlled
#' #    7641       1713        646          0          0 
#' set.seed(20250101L)
#' DNS_Engaged(10000L, RBF_Bias=0.5, Net_Type='K', Net_fPara=2.00)
#' 
#' # Test a (k=3, q=0.788) Kauffman model. It is critical 2p(1-p)K~1.
#' # set.seed(20250102L);
#' # DNS_Engaged(10000L, RBF_Bias=0.788, Net_Type='K', Net_fPara=3.00);
#' # Return $[[Overview]] (other are NA)
#' #  Stable    Useless    Engaged   External Controlled
#' #    7991       1171        838          0          0 
#' set.seed(20250102L)
#' DNS_Engaged(10000L, RBF_Bias=0.788, Net_Type='K', Net_fPara=3.00)
#' 
#' # Test a (k=3, q=0.5) Kauffman model. It is chaotic, 2p(1-p)K>1.
#' # set.seed(20250103L);
#' # DNS_Engaged(10000L, RBF_Bias=0.5, Net_Type='K', Net_fPara=3.00);
#' # Return $[[Overview]] (other are NA)
#' #   Stable    Useless    Engaged   External Controlled 
#' #      118        732       9150          0          0
#' set.seed(20250103L)
#' DNS_Engaged(10000L, RBF_Bias=0.5, Net_Type='K', Net_fPara=3.00)
DNS_Engaged<-function(Size=1000L,
  OBF_Type='R', OBF_iPara1=1L, OBF_iPara2=1L, OBF_Ratio=0.1, RBF_Bias=c(0.5,0.5),
  Net_Type='K', Net_fPara=4.00, Controller=NULL, ControllerValue=NULL,
  NodeAttri=FALSE, ResidualNet=FALSE, NumSys=2L){
  # Check parameter.
  f_par=c(Size, OBF_iPara1, OBF_iPara2, OBF_Ratio, RBF_Bias, Net_fPara);
  c_par=c(OBF_Type,Net_Type);
  if((!all(is.numeric(f_par)))||(!all(is.character(c_par)))){
    stop("Invalid inputs. Please check the help documentation.\n");}
  if((0>Size)||(0>OBF_Ratio||OBF_Ratio>1)||
  any(RBF_Bias<0)||any(RBF_Bias>1)||(NumSys>2&&length(RBF_Bias)!=NumSys)||
    (0>Net_fPara||Net_fPara>12)||# OBF_iPara1, OBF_iPara2 checked in other fun.
    !(OBF_Type %in% c('R','C','P','M','D','T'))||
    !(Net_Type %in% c('K','E','R','L'))){
    stop("Invalid inputs. Please check the help documentation.\n");}
  # Check existing some manual controlled nodes >>>
  con.id=Controller;
  con.vals=ControllerValue;
  if(!is.null(Controller)){# Has controllers.
    con.id=Controller-1;# Note the differences in the index between C++ and R.
    if(is.null(ControllerValue)){# Not provided, random setting.
      con.vals[Controller]=runif(length(Controller))>0.5;}
    else {# Provided.
      if(length(Controller)<=length(ControllerValue)){# Enough long
        con.vals[Controller]=as.logical(ControllerValue)[1:length(Controller)];}
      else {
        stop("Length of provided 'ControllerValue' is insufficient.");}}}
  else {# Non controlling con.id[0] serves as label.
    con.id=-666L;con.vals=c(-666,-666);}
  # Execute analysis >>>
  xx=c_ScalingLaw_Simualtion(as.integer(Size), as.integer(NumSys),
    OBF_Type[1], RBF_Bias, OBF_Ratio, Net_Type[1], Net_fPara, 
    as.integer(OBF_iPara1), as.integer(OBF_iPara2), as.integer(con.id), 
    as.integer(con.vals), as.integer(NodeAttri), as.integer(ResidualNet));
  # Set return list's name.
  names(xx)=c("Overview","Node_S01U","Node_SUE","ResidualNetwork","StableNodeInfo");
  names(xx[[1]])=c("Stable","Useless","Engaged","External","Controlled");
  if(ResidualNet&&(xx[[1]][3]>0)){# Return ResNet.
    xx[[5]][[1]]=paste0("N",formatC(c(1:Size)-1, width=floor(log10(Size-1))+1, flag="0"));
    names(xx[[5]][[2]])=xx[[5]][[1]];
    names(xx[[5]][[3]])=xx[[5]][[1]];
    names(xx[[5]][[4]])=xx[[5]][[1]];}
  return (xx);
}
