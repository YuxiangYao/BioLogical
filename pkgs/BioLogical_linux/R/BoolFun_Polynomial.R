#' @title Convert a Boolean function into threshold-based or polynomial forms
#' @description Boolean functions are typically expressed as logical expressions
#' involving variables and logical operators such as AND, OR, and NOT. When 
#' represented in arithmetic form, a Boolean function corresponds to a unique 
#' polynomial, including higher-order terms. It should be noted that the term 
#' \emph{threshold-based} used here does not conform to the standard definition,
#' which usually refers to models restricted to first-order terms.
#' @param Bit_Vec A Boolean vector (or 0/1-valued \code{IntegerVector}) represents the 
#' truth table of a Boolean function.
#' @param SpinLikeForm A logical value. Should the Boolean function be converted
#' into a spin-like form? (Default: \code{FALSE})
#' @param PolyForm A logical value. Should be converted according to a strict 
#' polynomial form.
#' @details Boolean functions are typically represented in two standard forms:
#' the binary form (0/1) and the spin-like form (+1/-1), which are linearly 
#' interconvertible. \code{SpinLikeForm} allows users to obtain representations 
#' in the spin-like form, while \code{PolyForm} provides an option to output 
#' the strict polynomial representation of a Boolean function. For example, for 
#' the expression "\eqn{f=a|b}", the polynomial form is "\eqn{f=a+b+ab}", and 
#' the threshold form is "\eqn{f=(a+b)>0?}". See the provided examples for details.
#' In the output, variables (\eqn{x_1}, \eqn{x_2}, ... , \eqn{x_k}, for a 
#' k-input function) represent the input bits ordered from the low bit to the high bit.
#' @return A three-element formatted list:
#' \itemize{
#' \item [[1]] Is Boolean satisfiability? (\code{TRUE} or \code{FALSE})
#' \item [[2]] The weight value of each variable and variable combination.
#' \item [[3]] The highest order of terms contained in the function.
#' }
#' @export 
#' @examples
#' # Convert the a Boolean function into a polynomial expression:
#' # Case1: Boolean threshold
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),F,F);
#' # Case2: Spin-like threshold
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),T,F);
#' # Case3: Boolean strict polynomial form
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),F,T);
#' # Case4: Spin-like strict polynomial form
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),T,T);
#' # Return [[1]] sat:TRUE, [[2]] (see below), [[3]] HighOrder: 2 or 3
#' # Detail Weight of four types:
#' #         x_1    x_2    x_3  x_1x_2  x_1x_3  x_2x_3 x_1x_2x_3  thershold 
#' # Case1     1      1      1      -2      -2      -2      NULL          0 
#' # Case2    -2     -2     -2      -2      -2      -2      NULL          0 
#' # Case3     1      1      1      -2      -2      -2         3          0 
#' # Case4 -0.25  -0.25  -0.25   -0.25   -0.25   -0.25      0.75      -0.25
#' 
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), FALSE, FALSE) # Case1
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), TRUE,  FALSE) # Case2
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), FALSE, TRUE ) # Case3
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), TRUE,  TRUE ) # Case4
#' 
BoolFun_Polynomial<-function(Bit_Vec, SpinLikeForm=FALSE, PolyForm=FALSE){
  tmp=r_CheckValidBoolFun(Bit_Vec);
  if(SpinLikeForm){
    bitmap=2*(tmp$bits)-1;# {-1,+1} form
  } else {
    bitmap=tmp$bits;}# The {0,1} form 
  in_deg=tmp$kin;
  ttt=r_GenerateTruthTable(inDegree=in_deg,SpinLike=SpinLikeForm);
  ttt=ttt[,c(in_deg:1)];# Low bit should place the first.
  sub_1=function(vec,TruthTable){return (apply(TruthTable[,vec],1,prod));}
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
  Res=list(sat=FALSE);whole.index=NULL;
  PolyFormNum=as.integer(2*SpinLikeForm+PolyForm);
  for(ii in c(1:in_deg)){
    tmp.index=sub_2(ttt,in_deg,ii);
    whole.index=cbind(whole.index,tmp.index);#cat(whole.index,"\n");
    tmp.analyze=c_BoolFun2Polynomial(whole.index, bitmap, PolyFormNum);
    if(tmp.analyze[[1]]){
      Wights=tmp.analyze[[2]];
      VarNames=NULL;
      for(jj in c(1:ii)){
        VarNames=c(VarNames,sub_3(tmp$kin,jj));}
      Wights=c(Wights,tmp.analyze[[3]]);
      VarNames=c(VarNames,"thershold");
      names(Wights)=VarNames;
      Res$sat=TRUE;
      Res$Wights=Wights;
      Res$HighOrder=ii;
      break;}}
  return (Res);
}