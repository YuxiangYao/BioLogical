#' @title Convert a Boolean function as threshold-based/polynomial forms
#' @description Usually, Boolean functions are presented as logical expressions
#' including variables and symbols (AND, OR, NOT). When limited to arithmetic 
#' denotations, a Boolean function can exist only one strict polynomial forms and 
#' infinite threshold-based interpretations. Sometimes, it can be involved with 
#' higher-order terms. Please note that this "threshold-based" is not standard 
#' definition of threshold-based one that only contains ist-order terms.
#' @param Bit_Vec a bool (or 0/1-integral) vector, represents a truth table of Boolean function.
#' @param SpinLikeForm logical, convert the Boolean function into a spinlike form based or not?
#' @param PolyForm logical, should be strict form of polyforms.
#' @details Boolean functions are generally represented in standard form (0/1) and 
#' spin-like form (+1/-1). They can be linearly transformed into each other. 
#' \code{SpinLikeForm} offers users with an option to output representations in 
#' spin-like form. \code{PolyForm} offers an option to output strict polynomial 
#' form of Boolean function. For instance, a|b, poly-form: a(1-b)+(1-a)b+ab=a+b+ab; 
#' threshold-form: theta(a+b)>0. See examples. In the return, x_1, x_2, ... , x_k 
#' (k-input) present low to high bits.
#' @return a List, [[1]] is.SAT ture or not? [[2]] Weights of coupled variable. [[3]] Highest order number.
#' @export 
#' @examples
#' # Convert the a Boolean function into a polynomial expression (4 cases).
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),F,F);# Boolean threshold (case1)
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),T,F);# Spin-like threshold (case2)
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),F,T);# Boolean strict poly-forms (case3)
#' # BoolFun_Polynomial(c(0,1,1,0,1,0,0,0),T,T);# Spin-like strict poly-forms (case4)
#' # Return [[1]] sat:TRUE, [[2]] (see below), [[3]] HighOrder: 2 or 3
#' # Detail Weight of four types:
#' #         x_1       x_2       x_3    x_1x_2    x_1x_3    x_2x_3 x_1x_2x_3  thershold 
#' # case1     1         1         1        -2        -2        -2      NULL          0 
#' # case2    -2        -2        -2        -2        -2        -2      NULL          0 
#' # case3     1         1         1        -2        -2        -2         3          0 
#' # case4 -0.25     -0.25     -0.25     -0.25     -0.25     -0.25      0.75      -0.25
#' 
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), FALSE, FALSE) # Case1
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), TRUE, FALSE) # Case2
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), FALSE, TRUE) # Case3
#' BoolFun_Polynomial(c(0,1,1,0,1,0,0,0), TRUE, TRUE) # Case4
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