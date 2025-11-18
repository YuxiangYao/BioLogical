# Fig2c: Derrida damage speard..

model3=rbind(
  read.csv("./simulation_analysis/dss_3_1_1_250422.csv",header = F),
  read.csv("./simulation_analysis/dss_3_2_1_250423.csv",header = F),
  read.csv("./simulation_analysis/dss_3_3_1_250424.csv",header = F),
  read.csv("./simulation_analysis/dss_3_1_2_250425.csv",header = F))

model4=rbind(
  read.csv("./simulation_analysis/dss_4_1_1_250426.csv",header = F),
  read.csv("./simulation_analysis/dss_4_2_1_250427.csv",header = F),
  read.csv("./simulation_analysis/dss_4_3_1_250428.csv",header = F),
  read.csv("./simulation_analysis/dss_4_1_2_250429.csv",header = F))

model3=cbind.data.frame(model3[,4:6],t=rep(letters[1:4],each=51));
model4=cbind.data.frame(model4[,4:6],t=rep(letters[1:4],each=51));
colnames(model3)=colnames(model4)=c('x','y','p','t');

library(ggplot2)

ggplot(data=model3,mapping = aes(x=x,y=y,col=t))+
  geom_ribbon(aes(ymin = y-p, ymax = y+p,fill=t), 
              alpha = 0.30, linewidth = 0) +
  geom_line(linewidth=1.35)+
  scale_color_manual(values = c("#F5B800","#FF7BAC","#3DB8B8","#b19e8b"))+
  scale_fill_manual(values = c("#F5B800","#FF7BAC","#3DB8B8","#b19e8b"))+
  theme_minimal()

ggsave(filename = "./fig2ca.pdf",width = 4.5,height = 2.3,units = "in")

ggplot(data=model4,mapping = aes(x=x,y=y,col=t))+
  geom_ribbon(aes(ymin = y-p, ymax = y+p,fill=t), 
              alpha = 0.30, linewidth = 0) +
  geom_line(linewidth=1.35)+
  scale_color_manual(values = c("#F5B800","#FF7BAC","#3DB8B8","#b19e8b"))+
  scale_fill_manual(values = c("#F5B800","#FF7BAC","#3DB8B8","#b19e8b"))+
  theme_minimal()
ggsave(filename = "./fig2cB.pdf",width = 4.5,height = 2.3,units = "in")

# 7.5*3 inch
