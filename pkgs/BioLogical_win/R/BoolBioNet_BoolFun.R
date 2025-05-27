#' @title Explore features of Boolean functions in biological/genetic networks
#' @description This function can calculate and return the categories, sensitivity, effectiveness, and complexity of each Boolean function in a given Boolean network. 
#' @param OneBioNet a 4-element formatted list, record Boolean network information.
#' @param Type Which statistical method is used (\code{rough} or \code{detail})
#' @details Codes can identify function types, including canalized, threshold, clique,
#' dominating, effective, monotone, signed, spin-like. The definations of various funciton can be seen 
#' the document of \link{BoolFun_Type}. Please note that single regulatory and external factor are 
#' ignored from statistics (only 2 <= in-degree <= 16, other cases return \code{NA}). When \code{Type} is 
#' \code{rough}, function returns a \code{data.frame} that contains statistics of each "in-degree" cluster. 
#' When \code{Type} is \code{detail}, the returned data.frame contains each case of Boolean function.
#' @return \code{data.frame}, counts the number of various types of functions (\code{Type} is \code{rough}).
#' \code{data.frame}, records the specific cases of each function within systems (\code{Type} is \code{detail}).
#' @export
#' @examples
#' # Analysis of one real Boolean genetic network.
#' # BoolBioNet_BoolFun(BoolGRN_CellCollective$c_1778,"rough");
#' # Return a data.frame:
#' #                 Number Canalized Clique Monotone Signed Threshold
#' # InDegree_1           1        NA     NA       NA     NA        NA
#' # InDegree_2           3         3      3        1      3         3
#' # InDegree_3           3         3      2        0      3         3
#' #                             ... ...
#' BoolBioNet_BoolFun(BoolGRN_CellCollective$c_1778,"rough")
#' 
#' # BoolBioNet_BoolFun(BoolGRN_CellCollective$c_1778,"detail");
#' # Return a data.frame
#' #            Input Canalized Clique Monotone Signed Threshold Sensitivity Effectiveness Complexity
#' # ADD            5      TRUE  FALSE    FALSE   TRUE      TRUE   1.3125000      1.802083        3.0
#' # ATM            3      TRUE  FALSE    FALSE   TRUE      TRUE   1.2500000      1.562500        2.0
#' # ATR            4      TRUE  FALSE    FALSE   TRUE      TRUE   1.2500000      1.625000        2.5
#' #                                ... ...
#' BoolBioNet_BoolFun(BoolGRN_CellCollective$c_1778,"detail")
#' 
BoolBioNet_BoolFun<-function(OneBioNet,Type=c("rough","detail")){
  # Boolean networks not distinguish biologically genetic or manually generated systems.
  # RealBioNet should be standard four-list format data that records information of bool-nets.
  # Check rationality of Boolean networks ...
  if(!BoolBioNet_Checker(OneBioNet)){
    stop("Invalid Boolean network! Analysis stops here.\n");}
  n.bf=length(OneBioNet[[1]]);
  t.bf=OneBioNet[[4]];
  if("rough"==Type[1]){
    Res=matrix(0,18,6);# Number, Ca, Th, Mo, Sign, Spin-like
    colnames(Res)=c("Number", "Canalized","Clique","Monotone","Signed","Threshold");
    Res=as.data.frame(Res);
    for(ii in c(1:n.bf)){
      if(!is.na(t.bf[[ii]][1])){# Have input 
        bn_tmp=OneBioNet[[2]][[ii]];
        k.var=length(bn_tmp);# Get in-degree
        if((1<k.var)&&(k.var<17)){# Input larger than 2
          Res[k.var,1]=Res[k.var,1]+1;
          bnbn=as.logical(OneBioNet[[4]][[ii]]);
          Res[k.var,2]=Res[k.var,2]+c_BF_isPointed(bnbn,k.var,'C',FALSE);
          Res[k.var,3]=Res[k.var,3]+c_BF_isPointed(bnbn,k.var,'P',FALSE);
          Res[k.var,4]=Res[k.var,4]+c_BF_isPointed(bnbn,k.var,'M',FALSE);
          Res[k.var,5]=Res[k.var,5]+c_BF_isPointed(bnbn,k.var,'S',FALSE);
          Res[k.var,6]=Res[k.var,6]+c_BF_isPointed(bnbn,k.var,'T',FALSE);}
        else if(k.var>=17){# Larger than 16.
          Res[17,1]=Res[17,1]+1;}
        else {# Single input.
          Res[1,1]=Res[1,1]+1;}}
      else {# Null inputs
        if(length(OneBioNet[[2]][[ii]])>1){# Some times, may exist network information.
          Res[17,1]=Res[17,1]+1;}
        else {# Real external control nodes.
          Res[18,1]=Res[18,1]+1;}}}
    # Formatting returned data.frame.
    Res[1,c(2:6)]=NA; Res[c(17:18),c(2:6)]=NA;
    rownames(Res)=c(paste0("InDegree_",c(1:16)),"Larger16","External_Factor");
  } else if("detail"==Type[1]){
    Res=matrix(NA,n.bf,9);# in-degree, {Ca, Th, Mo, Sign, Sp}, {Sensi, Effec, Compl}.
    rownames(Res)=OneBioNet[[1]];
    colnames(Res)=c("Input", "Canalized","Clique","Monotone","Signed","Threshold",
      "Sensitivity","Effectiveness","Complexity");
    Res=as.data.frame(Res);
    for(ii in c(1:n.bf)){
      if(!is.na(t.bf[[ii]][1])){# Have input 
        bn_tmp=OneBioNet[[2]][[ii]];
        k.var=length(bn_tmp);
        if(k.var>1){# Input larger than 2
          Res[ii,1]=k.var;
          bnbn=as.logical(OneBioNet[[4]][[ii]]);
          # Boolean function types.
          Res[ii,2:6]=c(c_BF_isPointed(bnbn,k.var,'C',FALSE), 
            c_BF_isPointed(bnbn,k.var,'P',FALSE), c_BF_isPointed(bnbn,k.var,'M',FALSE), 
            c_BF_isPointed(bnbn,k.var,'S',FALSE), c_BF_isPointed(bnbn,k.var,'T',FALSE));
          # Boolean function feeatures (Due to CP-complex, here only consider <=10 indegree cases.)
          if(k.var<=10){
            Res[ii,7:9]=c(c_BF_Sensitivity(bnbn,k.var),c_BF_Effective(bnbn,k.var),
              c_BF_Complexity(bnbn,k.var));
          } else {
            Res[ii,7:9]=c(c_BF_Sensitivity(bnbn,k.var),NA,NA);}}
        else {
           Res[ii,1]=1;}}
      else Res[ii,1]=0;}
  } else {
    stop("Incorrect type name of analysis (\"rough\",\"detail\"). See help.\n");}
  return (Res);
}
