#' @title Convert a multi-valued function as polynomial like forms
#' @description Unlike the Boolean system, we adopt a definition similar to 
#' that of domaineted function; please refer to \link{MulVFun_is_Domainted}.
#' @param aVec An IntegerVector that represents a mapping table of multi-valued 
#' function.
#' @param k An integer that denotes the number of input variables (Default: 2).
#' @param L An integer that represents the level of multi-valued system (Default: 3).
#' @return A list containing three elements:
#' \itemize{
#' \item \code{[[1]] sat}: Is satisfiability? \code{TRUE} or \code{FALSE}
#' \item \code{[[2]] Wights}: The weights and thresholds. 
#' \item \code{[[3]] HighOrder}: The highest order of terms contained in the 
#' function (See examples).}
#' @export 
#' @examples
#' # Convert the a multi-valued function into a polynomial like expression.
#' # In fact, the result is same as MulVFun_is_Domainted().
#' # MulVFun_Polynomial(c(0,1,2,1,1,2,2,1,2), k=2L, L=3L);
#' # Return 
#' # [[1]] sat:TRUE, 
#' # [[2]] Weight:
#' #    x_1         x_2 threshold_1 threshold_2 threshold_3 
#' #      2           1          -2           0           0
#' # [[3]] HighOrder: 1
#' 
#' MulVFun_Polynomial(c(0,1,2,1,1,2,2,1,2), k=2L, L=3L)
#' 
MulVFun_Polynomial<-function(aVec, k=2, L=3){
  in_deg=as.integer(k);
  l_sys=as.integer(L);
  tmp=r_CheckValidMulVFun(a_Vec=aVec, k=in_deg, L=l_sys);
  ttt=r_GenerateTruthTable_M(inDegree=in_deg, Lsystem=l_sys);
  ttt=ttt[,c(in_deg:1)];# Low bit should place the first.
  # Return high-ordered results (Here uses a MOD operator)
  sub_1=function(vec,TruthTable){return(apply(TruthTable[,vec], 1, 
    function(xx,L){return(sum(xx)%%L);}, L=l_sys));}# <--- use MOD (%%)
  # In various order levels.
  sub_2=function(TruthTable,InDeg,nChoose){# Obtain i-order coupled truth tables.
    if(1L==nChoose){
      return (TruthTable);
    } else {
      index=t(combn(InDeg,nChoose));
      return (apply(index,1,sub_1,TruthTable=TruthTable));}}
  sub_3=function(InDeg,nChoose){# Set variable names.
    if(1L==nChoose){
      VarName=paste0("x_",c(1:InDeg));
    } else {
      ids=t(combn(InDeg,nChoose));
      ids=apply(ids,c(1,2),function(x){paste0("x_",x);});
      VarName=apply(ids,1,paste,collapse="");}
    return (VarName);}
  # Repeated check the order of ordered polynomial vector.
  Res=list(sat=FALSE);
  whole.index=NULL;
  for(ii in c(1:in_deg)){
    tmp.index=sub_2(ttt,in_deg,ii);
    whole.index=cbind(whole.index,tmp.index);
    tmp.analyze=c_MulVFun2Polynomial(whole.index, tmp, in_deg, l_sys);
    if(tmp.analyze[[1]]){
      Wights_Thres=c(tmp.analyze[[2]], tmp.analyze[[3]]);
      VarNames=NULL;
      for(jj in c(1:ii)){
        VarNames=c(VarNames,sub_3(in_deg, jj));}
      VarNames=c(VarNames,paste0("threshold_",(c(1:l_sys)-0)));
      names(Wights_Thres)=VarNames
      Res$sat=TRUE;
      Res$Wights=Wights_Thres;
      Res$HighOrder=ii;
      break;}}
  return (Res);
}