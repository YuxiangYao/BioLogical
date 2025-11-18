#' @title Decode the structure of a nested canalized Boolean function
#' @description A nested canalized Boolean function typically exhibits an "if-else
#' if..." regulatory structure. This function decodes and visually represents the
#' pattern to illustrate its logical hierarchy intuitively.
#' @param Bit_Vec A Boolean vector (or 0/1-valued \code{IntegerVector}) represents
#' the truth table of a Boolean function.
#' @param PrintOut A logical value. Should output be displayed in the terminal? 
#' (Default: \code{FALSE}).
#' @return A matrix that descirbes the nested structure of a canalized function.
#' @export
#' @examples
#' # Show the canalized structure of a nested canalized function.
#' # set.seed(202501);
#' # aboolfun<-BoolFun_Generator("C", VarNum=7L, vars=c(4L,9L));
#' # df_boolfun<-BoolFun_NestedCana(aboolfun, PrintOut=TRUE);
#' # Return a matrix and display the nested relationships.
#' # if {x_5=1} ==> f(x)=0
#' # else if {x_5=0, x_6=1} ==> f(x)=0
#' # else if {x_5=0, x_6=0, x_7=1} ==> f(x)=0
#' # else if {x_5=0, x_6=0, x_7=0, x_3=1} ==> f(x)=1
#' # else if {x_5=0, x_6=0, x_7=0, x_3=0, x_4=1} ==> f(x)=0
#' # df_boolfun is:
#' #    x_1 x_2 x_3 x_4 x_5 x_6 x_7 Out
#' #  1  NA  NA  NA  NA   1  NA  NA   0
#' #  2  NA  NA  NA  NA   0   1  NA   0
#' #  3  NA  NA  NA  NA   0   0   1   0
#' #  4  NA  NA   1  NA   0   0   0   1
#' #  5  NA  NA   0   1   0   0   0   0
#' 
#' set.seed(202501)
#' aboolfun<-BoolFun_Generator("C", VarNum=7L, vars=c(4L,9L))
#' df_boolfun<-BoolFun_NestedCana(aboolfun, PrintOut=TRUE)
#' 
BoolFun_NestedCana<-function(Bit_Vec, PrintOut=FALSE){
  if(!is.numeric(Bit_Vec)&&!is.integer(Bit_Vec)&&!is.logical(Bit_Vec)){
    stop("Invalid mapping table input. The table should be numeric, integer or logical.\n");}
  else {
    bit_vec=as.logical(Bit_Vec);}
  if(length(bit_vec)<=0||log2(length(bit_vec))>16){
    stop("Invalid length of bit_vec. The number of variables is less than 16.\n");}
  else {
    tmp=log2(length(bit_vec));
    if(abs(tmp-round(tmp))>1e-7){
      stop("Invalid length of bit_vec. The length should be 2^k.\n");}
    vars=as.integer(tmp);
    if(!c_BF_isPointed(as.logical(bit_vec),vars,"C",FALSE)){
      stop("This function is not a canalized function\n");}
    xx=c_B_NestedCanalized(as.logical(bit_vec),vars);
    n_item=length(xx[[1]]);
    yy=matrix(NA,n_item,vars+1);
    yy[,vars+1]=xx[[3]];
    for(ii in c(1:n_item)){
      yy[ii,xx[[1]][ii]+1]=xx[[2]][ii];# ii-layer, which variable denoted!
      if(ii<n_item){
        yy[(ii+1):n_item,xx[[1]][ii]+1]=1L-xx[[2]][ii];}}
    colnames(yy)=c(paste0("x_",c(1:vars)),"Out");
    yy=as.data.frame(yy);
    if(PrintOut){
      outers=NULL;
      reorder.id=unique(xx[[1]]+1);
      maptab=yy[,reorder.id];
      fff_bn=yy[,vars+1];
      if(1==length(reorder.id)){# Only one variable (single/double canalized values) 
        outers=paste0(outers,"if {x_",reorder.id,"=",maptab[1],"} ==> f(x)=",fff_bn[1],"\n");
        if(1<length(maptab)){
          outers=paste0(outers,"else if {x_",reorder.id,"=",maptab[2],"} ==> f(x)=",fff_bn[2],"\n");}
      } else {
        for(ii in c(1:nrow(yy))){
          id=!is.na(maptab[ii,]);
          tmp=paste(paste0(colnames(maptab)[id],"=",maptab[ii,id]), collapse=", ");
          if(ii>1){
            outers=paste0(outers,"else ");}
          outers=paste0(outers,"if {",tmp,"} ==> f(x)=",fff_bn[ii],"\n");}}
      cat(outers);}
    return (yy);}
}
