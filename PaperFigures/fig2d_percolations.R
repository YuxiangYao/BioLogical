# Fig2d: percolation
library(BioLogical);
set.seed(54321L);
res=DNS_Percolation(Size = 2500L,2500L,1000L,
                           OutPutState = TRUE,
                       OBF_Type = "C",
                       OBF_iPara1 = 3L,
                       OBF_iPara2 = 2L,
                       OBF_Ratio = 0.85,
                       LatType=4,Net_fPara = 4,
                       RBF_Bias = c(1,1,1)/3.0,
                       NumSys = 3L)

stasss=letters[(res[[2]]+(res[[3]]*3+1))+1]
FigLatticePlot(stasss,50L,50L,
               Type=4L,
               # ColorFills = c("#fefefe","#f8e0a0","#f4b643","#f39c12",
               #                "#73B1C9","#2B7E9C","#1B5B7E"),
               #ColorFills = c("#fefefe","#E8D5B5","#B3987B","#82725B",
               #                "#A8B86B","#799956","#456A2F"),
               ColorFills = c("#fefefe","#FDC9C1","#FA8C7E","#ED1C24",
                              
                              "#B0CEE3","#5796C4","#006699"),
               Linewidth = 0.25)
library(ggplot2)
ggsave("./fig2d.pdf",width=4,height = 4,units="in")
