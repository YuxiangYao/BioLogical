# Fig2c:

args <- commandArgs(trailingOnly = TRUE)

param1 <- as.numeric(args[1])# Number system
param2 <- as.numeric(args[2])# Canalized Deep
param3 <- as.numeric(args[3])# Canalized Variable in each layer
#param4 <- as.numeric(args[4])# Canalized ratio [1~50], bin:0.02
param4 <- as.numeric(args[4])# Rand seed

library(BioLogical)
set.seed(param4);
res=matrix(NA,51,6);
seed_1=0;
for(ii in c(0:50)){
  tmp=rep(NA,1000);
  for(jj in c(1:1000)){
    tmp[jj]=DNS_DamageSpread(
      Size = 1000L,
      SimStep = 2000L,
      Init_Dist = 0.1,
      Init_1_Ratio = rep(1.0,param1)/param1,
      OBF_Type = "C",
      OBF_iPara1 = param2,
      OBF_iPara2 = param3,
      OBF_Ratio = 0.02*ii,
      RBF_Bias = rep(1.0,param1)/param1,
      NumSys = param1);
    seed_1=seed_1+1;
  }
  res[(ii+1),]=
    c(param1,param2,param3,0.02*ii,mean(tmp),round(sd(tmp),6));
}

# outputs=c(param1,param2,param3,param4,param5,mean(tmp),sd(tmp))
finename=paste(args,collapse = "_");
finename=paste0("dss_",finename,".csv")
# writeLines(paste0(paste(as.character(outputs),collapse =","),"\n"), finename)

write.table(res,finename,sep=",",quote = FALSE,
            row.names = FALSE,col.names = FALSE)
