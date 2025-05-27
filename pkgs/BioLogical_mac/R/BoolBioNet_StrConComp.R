#' @title Caluclate network's strong connect component
#' @description This function briefly analyze the strong connect component for 
#' comparsion of core dynamic anslysis.
#' @param OneNet a 4-element formatted list, record Boolean network information.
#' @return List, each element denotes one strong connect component.
#' @export
#' @examples
#' # Analyze the strong connect componets of a network.
#' # Here show an example of net[c_2202] in Cell Collective.
#' # BoolBioNet_StrConComp(BoolGRN_CellCollective$c_2202)
#' # Return a List
#' # [[1]]
#' # [1] "Exocytosis"
#' # [[2]]
#' # [1] "Packaging_Proteins"
#' # [[3]]
#' # [1] "Protein_Phosphatase_1" "DARPP32"  "Calcineurin" 
#' # [4] "Calcium"               "Glutamate_Receptor"   
#' # ... ... 
#' BoolBioNet_StrConComp(BoolGRN_CellCollective$c_2202)
BoolBioNet_StrConComp<-function(OneNet){
  Size=length(OneNet[[1]]);
  ind=otd=rep(0L,Size);
  con.vals=rep(666L,Size);# Here is real controlled values.
  NonIputer=NULL;# Label recording.
  for(ii in c(1:Size)){
    if(is.na(OneNet[[2]][[ii]][1])||is.null(OneNet[[2]][[ii]][1])){
      ind[ii]=0L;
      NonIputer=c(NonIputer,ii);}# Insert Exponents
    else {
      ind[ii]=length(OneNet[[2]][[ii]]);}
    if(is.na(OneNet[[3]][[ii]][1])||is.null(OneNet[[3]][[ii]][1])){
      otd[ii]=0L;}
    else {
      otd[ii]=length(OneNet[[3]][[ii]]);}}
  # Conver to node name.
  tmp=c_StrongConnectComponent(OneNet,ind,otd);
  for(ii in c(1:length(tmp))){
    tmp[[ii]]=OneNet[[1]][tmp[[ii]]+1];
  }
  return (tmp);
}