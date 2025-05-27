library(BioLogical)
set.seed(30201L);# 
ttt=MulVFun_Generator("C",k=3,L=3, 
                  CanaDeep = 3, CanaVarNum =c(1,1,1), 
                  MappingTable = TRUE)

MulVFun_Sensitivity(ttt[,4],k=3,L=3)
MulVFun_EffectiveEdges(ttt[,4],k=3,L=3,Detail = TRUE)
MulVFun_Complexity(ttt[,4],k=3,L=3)

MulVFun_is_NestedCanalized(ttt[,4],3L,3L)

MulVFun_QMForm(ttt[,4],k=3L,L=3L,)



MulVFun_is_Domainted(ttt[,4],3L,3L)
MulVFun_is_Threshold(ttt[,4],3L,3L)
MulVFun_is_Signed(ttt[,4],3L,3L)
