#' @title Load a Boolean genetic network
#' @description Load genetic network information from external files.
#' @param NetName a valid file/address of target GenNet.
#' @param NetType Which input format? Limited to \code{BoolExpression}, 
#' \code{TruthTable}, \code{Threshold}.
#' \itemize{
#' \item \code{BoolExpression}: all genes in a single file are shown as Boolean expression like,
#' \preformatted{
#'   X_1 = X_3 AND X_4
#'   X_2 = NOT (X_1 OR X_3)
#'   ... ...
#' }
#' \item \code{TruthTable}: all genes stored in a specified directory like (Based on organizational forms of Cell Collective Web),
#' \preformatted{
#'   in NetName/X_1.csv (for gene X_1)
#'      X_3, X_5, X_1
#'       0 ,  0 ,  0
#'       0 ,  1 ,  1
#'       1 ,  0 ,  1
#'       1 ,  1 ,  0
#'   in NetName/X_2.csv (for gene X_2)
#'      X_4, X_2
#'       0 ,  0
#'       1 ,  1
#'       ... ...
#' }
#' \item \code{Threshold}: represent all regulations 1 for activation, 2 for inhibition.
#' \preformatted{
#'      X_1  X_3  1
#'      X_4  X_2  2
#'      X_5  X_2  1
#'       ... ...
#' }
#' }
#' The detail files and organizational structures can see \href{https://cellcollective.org}{Cell Collective Website}
#' @return a list contains member names of system, topological structure,
#' and Boolean function.
#' @export
#'
LoadBoolBioNetwork<-function(NetName,NetType=c("BoolExpression","TruthTable","Threshold")){
  BitTable<-function(VarNum){
    lens=bitwShiftL(1,VarNum)-1;
    # tabs=NULL;
    # for(ii in c(0:lens)){
    #   tmp=as.integer(as.vector(intToBits(ii)));
    #   tabs=rbind(tabs,tmp[VarNum:1]);}
    tabs=c_FrameTruthTable(VarNum);
    rownames(tabs)=c(0:lens);
    colnames(tabs)=paste0("v",(c(VarNum:1)-1));
    return(tabs);
  }
  if("TruthTable"==NetType[1]){
    files=dir(NetName);
    independen.factor=paste0(NetName,"/external_components.ALL.txt");
    if(file.info(independen.factor)$size>0){# [external_components.ALL.txt] is empty (no independent variables)
      exponter=as.vector(as.matrix(read.table(paste0(NetName,
        "/external_components.ALL.txt"),sep=",")));
      files=setdiff(files,"external_components.ALL.txt");
      nodes=unique(c(gsub(".csv","",files),exponter));}
    else {# [external_components.ALL.txt] is not empty (exist independent variables)
      files=setdiff(files,"external_components.ALL.txt");
      nodes=c(gsub(".csv","",files));}
    # Deal network information
    nodes=cbind.data.frame(gensys=nodes,id=as.integer(c(1:length(nodes))-1));
    rownames(nodes)=nodes[,1];
    sys.node=unique(nodes[,1]);
    empty.list=vector("list",length = length(sys.node));
    names(empty.list)=sys.node;
    inedge=empty.list;# Which input?
    otedge=empty.list;# Output which?
    boolfn=empty.list;# Pointed which?
    for(ii in files){
      id_c=gsub(".csv","",ii);
      bnf=read.csv(paste0(NetName,"/",ii),header = TRUE);
      tmp=colnames(bnf);
      loc=ncol(bnf);
      boolfn[[id_c]]=bnf[,loc];
      tmp=tmp[c(1:(loc-1))];
      inedge[[id_c]]=nodes[tmp,2];
      for(jj in tmp){
        otedge[[jj]]=c(otedge[[jj]],nodes[id_c,2]);}}
    for(ii in c(1:length(sys.node))){
      if(is.null(inedge[[ii]]))inedge[[ii]]=NA;
      if(is.null(otedge[[ii]]))otedge[[ii]]=NA;
      if(is.null(boolfn[[ii]]))boolfn[[ii]]=NA;}
  }
  else if("BoolExpression"==NetType[1]){# Folder form of TTT.
    regus=read.table(NetName,sep="@",header = F);# Some split-symbol is \t instead of space.
    nodes=NULL;
    for(ii in c(1:nrow(regus))){
      chrx=gsub("\\("," \\( ",regus[ii,],perl = TRUE);
      chrx=gsub("\\)"," \\) ",chrx,perl = TRUE);
      chrx=gsub("\t"," ",chrx);# Remove symbol TAB
      chrx=gsub("/", "_", chrx, perl = TRUE);
      chrx=gsub("\\.", "_", chrx, perl = TRUE);
      chrx=gsub("-", "_", chrx, perl = TRUE);
      nodes=c(nodes, unlist(strsplit(chrx," ")[[1]]) );
      # nodes=c(nodes, unlist(strsplit(regus[ii,]," ")[[1]]) );
    }
    nodes=setdiff(unique(nodes),c("=","NOT","not","AND","and","OR","or","(",")",""));
    nodes=cbind.data.frame(gensys=nodes,id=c(1:length(nodes))-1);
    rownames(nodes)=nodes[,1];
    sys.node=unique(nodes[,1]);
    empty.list1=vector("list",length = length(sys.node));# NULL;
    empty.list2=as.list(rep(NA,length(sys.node)))# NA
    names(empty.list1)=names(empty.list2)=sys.node;
    inedge=empty.list2;# Which input?
    otedge=empty.list1;# Output which?
    boolfn=empty.list2;# Pointed which?
    for(ii in c(1:nrow(regus))){
      # bnf=BoolExpTTT(regus[ii,]);
      bnf=BoolFun_Expression2MapTable(regus[ii,],TRUE);
      tmp=colnames(bnf);
      loc=length(tmp)
      id=tmp[loc];
      boolfn[[id]]=bnf[,loc];
      tmp=tmp[c(1:(loc-1))];
      inedge[[id]]=as.integer(nodes[tmp,2]);
      for(jj in tmp){
        otedge[[jj]]=c(otedge[[jj]],nodes[id,2]);}}
    for(ii in sys.node){
      if(is.null(otedge[[ii]])){
        otedge[[ii]]=NA;}
      else {
        otedge[[ii]]=as.integer(otedge[[ii]]);}}
  }
  else if("Threshold"==NetType[1]){# Threshold type.
    regus=read.table(NetName,sep="\t",header = T);
    nodes=unique(c(regus[,1],regus[,2]));
    nodes=cbind.data.frame(gensys=nodes,id=c(1:length(nodes))-1);
    rownames(nodes)=nodes[,1];
    sys.node=unique(nodes[,1]);
    empty.list=list();
    for(ii in sys.node){empty.list[[ii]]=NA;}
    inedge=empty.list;# Which input?
    otedge=empty.list;# Output which?
    boolfn=empty.list;# Pointed which?
    for(ii in sys.node){
      index=(regus[,2]==ii);
      if(0==sum(index)){
        inedge[[ii]]=NA;
      } else {
        inedge[[ii]]=as.integer(nodes[regus[index,1],2]);}
      if(0==sum(regus[,1]==ii)){
        otedge[[ii]]=NA;
      } else {
        otedge[[ii]]=as.integer(nodes[regus[regus[,1]==ii,2],2]);}
      tmp=regus[index,3];
      if(length(tmp)>0){
        tmp[2==tmp]=-1;
        boolfn[[ii]]=as.integer(BitTable(length(tmp))%*%tmp>0);
      }else {
        boolfn[[ii]]=NA;}}
  }
  else {
    stop("Invalid Boolean network type.\n Please choose 'BoolExpression', 'TruthTable', 'Thrshold'.\n");
  }
  if(!((length(sys.node)==length(inedge))&&
    (length(sys.node)==length(otedge))&&
    (length(sys.node)==length(boolfn)))){
    stop("Loaded Boolean network may exist abnormal information. 
      Please carefully check the original files to avoid unusual characters in factor/file names, 
      including '+', '-','(',')', etc.");}
  return(list(AllMember=sys.node,InEdge=inedge,OutEdge=otedge,BoolFun=boolfn));
}
