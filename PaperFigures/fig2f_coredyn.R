# Figure 2f: core dynamic samples
library(BioLogical,dplyr)
nets=BoolGRN_CellCollective;
set.cc=matrix(NA,length(nets),5);
for(ii in c(1:length(nets))){
  if(length(nets[[ii]][[1]])>20){
    set.cc[ii,]=BoolBioNet_CoreDyn(nets[[ii]])[[1]];
  }
}

nets=BoolGRN_KadelkaSet;
set.ka=matrix(NA,length(nets),5);
for(ii in c(1:length(nets))){
  if(length(nets[[ii]][[1]])>20){
    set.ka[ii,]=BoolBioNet_CoreDyn(nets[[ii]])[[1]];
  }
}


nets=BoolGRN_ThresholdModel;
set.th=matrix(NA,length(nets),5);
for(ii in c(1,2,3,4,5,6,7,8)){
  if(length(nets[[ii]][[1]])>20){
    set.th[ii,]=BoolBioNet_CoreDyn(nets[[ii]],Times = 4L)[[1]];
  }
}

coredyn=rbind(
  cbind.data.frame(set.cc,type="Cc"),
  cbind.data.frame(set.ka,type="Ka"),
  cbind.data.frame(set.th,type="Tr"))
coredyn=coredyn[!is.na(coredyn[,1]),c(1:3,6)];
coredyn=cbind.data.frame(n=paste0("N",sprintf("%02d", c(1:nrow(coredyn)))),
                         coredyn);
coredyn=cbind.data.frame(coredyn$n,rowSums(coredyn[,2:3]),coredyn[,4:5])
colnames(coredyn)=c("n","u","e","t");
reshape2::melt(coredyn)->tmp
library(dplyr)
tmp= tmp %>% group_by(n) %>%
  mutate(total = sum(value), proportion = value / total) %>%
  ungroup()
tmp=cbind.data.frame(tmp,col=NA,ttcol="#000000")
tmp$col[tmp$t=="Cc"&tmp$variable!="e"]="#66c2a5";# 
tmp$col[tmp$t=="Cc"&tmp$variable=="e"]="#003F5C";
tmp$col[tmp$t=="Ka"&tmp$variable!="e"]="#F88BA9";
tmp$col[tmp$t=="Ka"&tmp$variable=="e"]="#D45087";
tmp$col[tmp$t=="Tr"&tmp$variable!="e"]="#F7F056";
tmp$col[tmp$t=="Tr"&tmp$variable=="e"]="#B5813B";
tmp$ttcol[tmp$variable=="e"]="#ffffff"

library(ggplot2);
ggplot(tmp, aes(x = n, y = proportion, fill = col)) +
  geom_bar(stat = "identity", position = "stack") +
  geom_text(aes(label = value), 
            position = position_stack(vjust = 0.5), angle=90,
            color = tmp$ttcol) +
  labs(title = "core dynamic component",
       x = "Group",
       y = "Proportion") +
  theme_minimal() +
  theme(legend.position = "bottom",
        panel.grid=element_blank() ) +
  scale_fill_identity(guide = "legend",labels=LETTERS[1:9])


# 10*4.2 inch
