args <- commandArgs(trailingOnly=TRUE);
param1 <- as.integer(args[1]);# Bias of value-0
param2 <- as.integer(args[2]);# Bias of value-1
param3 <- as.integer(args[3]);# Bias of value-2
param4 <- as.integer(args[4]);# Random seed

res=matrix(NA,1000,3);
set.seed(param4+20000000L);
for(ii in c(1:1000)){
  BioLogical::DNS_Engaged(
    Size =10000L,
    Net_fPara = 2,Net_Type = "K",
    OBF_Type = "R",
    OBF_Ratio = 0,    
    RBF_Bias = c(param1,param2,param3)/300.00,
    NumSys = 3L,
    ResidualNet = FALSE)[[1]]->yyx0
 res[ii,]=yyx0[1:3]
}

mingzi=paste(c("slk3",param1,param2,param3,"f.txt"),collapse = "_");


wenben=paste(
  c(param1,param2,param3,colMeans(res)),collapse = ",")

writeLines(wenben, con = mingzi)
