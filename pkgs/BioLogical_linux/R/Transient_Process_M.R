#' @title Simulating roughly attractor analysis of multi-valued model
#' @description Like \link{Transient_Process_A}, but use multi-valued model. 
#' @param RealBioNet_M A four-element formatted list that records information of a multi-valued threshold network.
#' @param SampleNum An integer pointing the number of randomly sampling initial states.
#' @param Times an integer, steps of simulating system.
#' @param NumSys an integer, number of discrete value.
#' @details See the reference [\href{https://doi.org/10.1186/1752-0509-1-4}{Schaub2007}]
#' @return A matrix where each row represents the final state of a stochastic 
#' simulation, and each column represents a gene.
#' @export
#' 
Transient_Process_M<-function(RealBioNet_M, SampleNum=100L,
  Times=1000L, NumSys=2L){
  # Prepare for transmitting parameters to C-prototype function.
  RealBioNet=RealBioNet_M;
  Size=length(RealBioNet[[1]]);
  ind=otd=rep(0L,Size);
  All_Nega=rep(0L,Size);
  # Check no input nodes.
  for(ii in c(1:Size)){
    if(is.na(RealBioNet[[2]][[ii]][1])||is.null(RealBioNet[[2]][[ii]][1])){
      ind[ii]=0L;
      RealBioNet[[4]][[ii]]=0;
      }# Insert Exponents
    else {
      ind[ii]=length(RealBioNet[[2]][[ii]]);}
    if(is.na(RealBioNet[[3]][[ii]][1])||is.null(RealBioNet[[3]][[ii]][1])){
      otd[ii]=0L;}
    else {
      otd[ii]=length(RealBioNet[[3]][[ii]]);}
    All_Nega[ii]=all(RealBioNet[[4]][[ii]]<0);
  }
  aaa=rep(-1,2);# Useless vector
  bbb=rep(0,Size);# Useless vector
  # Execute analysis >>>
  xx=c_AttractorSim_MulVal(RealBioNet, as.integer(SampleNum),
    aaa, bbb, ind, otd, All_Nega, as.integer(Times), as.integer(NumSys));
  xx=t(xx);
  colnames(xx)=RealBioNet[[1]];
  return (xx);
}
