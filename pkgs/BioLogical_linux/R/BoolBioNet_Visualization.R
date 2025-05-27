#' @title Print-out/Visualize networks
#' @description Convert network into other standardized format or visualization 
#' [text.vector, data.frame, to.terminal or igraph.object (if available)].
#' @param RealBioNet a 4-element formatted list, record biological network information
#' @param Outtype Character, how the information will be presented? 
#' \code{print}: Display directly to the terminal; \code{text}: a string that 
#' save all information as \code{print} scenario; \code{dataframe}: a data.frame
#' that record "source-target" pairs; \code{igraph}: a igraph-object that visualize
#' the network (if package \code{igraph} is available).
#' @return Specific formatted data of networks
#' @export
#' @examples
#' # Print the network of c_1969 in Cell Collective set.
#' # BoolBioNet_Visualization(CellCollective$c_1969,"print");
#' # Print Out to terminal
#' #  "PI3K-> PTEN-> [Akt] ->Cdc42_Rac1"
#' #  "cFOS-> cJUN-> [AP1] ->uPAR"
#' #  "JNK-> p38-> [ATF2] ->CyclinD ->PTGS2"
#' #               ... ...
#' BoolBioNet_Visualization(BoolGRN_CellCollective$c_1969,"print")
#'
#' # BoolBioNet_Visualization(BoolGRN_CellCollective$c_1969,"text");
#' # Return a text vector
#' #                                Akt
#' # "PI3K-> PTEN-> [Akt] ->Cdc42_Rac1"
#' #                          AP1
#' # "cFOS-> cJUN-> [AP1] ->uPAR"
#' #          ... ...
#' BoolBioNet_Visualization(BoolGRN_CellCollective$c_1969,"text")
#'
#' # BoolBioNet_Visualization(BoolGRN_CellCollective$c_1969,"dataframe");
#' # Return a data.frame
#' #         source        target
#' # 1         PI3K           Akt
#' # 2         PTEN           Akt
#' # 3         cFOS           AP1
#' #          ... ...
#' BoolBioNet_Visualization(BoolGRN_CellCollective$c_1969,"dataframe")
#' 
BoolBioNet_Visualization<-function(RealBioNet, Outtype=c("print","text","dataframe","igraph")){
  genenames=RealBioNet$AllMember;
  n.size=length(genenames);
  # Output to terminal.
  if("print"==Outtype[1]){
    for(ii in c(1:n.size)){
      regulator=NULL;
      if(!is.na(RealBioNet$InEdge[[ii]][1])){
        for(jj in c(1:length(RealBioNet$InEdge[[ii]]))){
          regulator=paste0(regulator,genenames[RealBioNet$InEdge[[ii]][jj]+1],"-> ");}}
      regulator=paste0(regulator,"[",genenames[ii],"]");
      if(!is.na(RealBioNet$OutEdge[[ii]][1])){
        for(jj in c(1:length(RealBioNet$OutEdge[[ii]]))){
          regulator=paste0(regulator," ->",genenames[RealBioNet$OutEdge[[ii]][jj]+1]);}}
      regulator=paste0(regulator,"\n");
      cat(regulator);}}
  # Convert into a text vector.
  else if("text"==Outtype[1]){
    Returner=rep(NA,n.size);
    for(ii in c(1:n.size)){
      regulator=NULL;
      if(!is.na(RealBioNet$InEdge[[ii]][1])){
        for(jj in c(1:length(RealBioNet$InEdge[[ii]]))){
          regulator=paste0(regulator,genenames[RealBioNet$InEdge[[ii]][jj]+1],"-> ");}}
      regulator=paste0(regulator,"[",genenames[ii],"]");
      if(!is.na(RealBioNet$OutEdge[[ii]][1])){
        for(jj in c(1:length(RealBioNet$OutEdge[[ii]]))){
          regulator=paste0(regulator," ->",genenames[RealBioNet$OutEdge[[ii]][jj]+1]);}}
      Returner[ii]=regulator;}
    names(Returner)=genenames;
    return(Returner);}
  # Data frame or igraph object
  else if(Outtype[1]%in%c("dataframe","igraph")){
    node.source=NULL;
    node.target=NULL;
    ExternalFactor=rep(FALSE,n.size);
    TerminalFactor=rep(FALSE,n.size);
    IsoformsFactor=rep(FALSE,n.size);
    # Repetitive checking
    for(ii in c(1:n.size)){
      tmp=RealBioNet$InEdge[[ii]];# Only consider in-degree cases.
      if(!is.na(tmp[1])){# Has inputs.
        node.source=c(node.source,genenames[tmp+1]);
        node.target=c(node.target,rep(genenames[ii],length(tmp)));
        if(is.na(RealBioNet$OutEdge[[ii]][1])){# No output.
          TerminalFactor[ii]=TRUE;}}# XXX -> node -> NULL
      else {# No inputs
        #ExternalFactor[ii]=TRUE;
        if(is.na(RealBioNet$OutEdge[[ii]][1])){# No output.
          IsoformsFactor[ii]=TRUE;}# NULL -> node -> NULL
        else {# Has output.
          ExternalFactor[ii]=TRUE;}}}# NULL -> node -> XXX
    # Output formatted data.
    if("dataframe"==Outtype[1]){
      Returner=cbind.data.frame(source=node.source,target=node.target);}
    else {
      if(!requireNamespace("igraph", quietly=TRUE)){
        stop("The 'igraph' package is required but not installed.\n");}
      edges=cbind.data.frame(source=node.source,target=node.target);
      nodes=cbind.data.frame(name=genenames,Nulls=ExternalFactor);
      Returner=igraph::graph_from_data_frame(d=edges,directed=TRUE,vertices=nodes);
      colsss=rep("#acacac",n.size);
      colsss[ExternalFactor]="#EE312E";# External factors.
      colsss[TerminalFactor]="#2156A6";# Terminal nodes.
      colsss[IsoformsFactor]="#FFB600";# Isoformed nodes (for residual/reduced network).
      igraph::V(Returner)$size=8;
      igraph::V(Returner)$color=colsss;
      #igraph::V(Returner)$alpha = 0.8;
      igraph::E(Returner)$width=1.5;
      igraph::E(Returner)$color="#080808";
      igraph::E(Returner)$arrow.type=1;
      igraph::E(Returner)$arrow.size=0.2;}
    return(Returner);}
  else {stop("Unknown or invalid data format.");}
  # Recommand:
  # plot(Gve,vertex.label.color = "#000000",vertex.label.cex = 0.55,
  #   layout=igraph::layout_as_tree/igraph::layout_with_fr);
}
