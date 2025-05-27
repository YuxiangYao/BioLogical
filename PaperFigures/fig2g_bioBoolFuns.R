# Fig2g: Obtain sensitivities & effective edges: 

library(boolGenDyn)
nets=BoolGRN_CellCollective;
set.cc=NULL;
for(ii in c(1:length(nets))){
  tmp.ind=nets[[ii]][[2]];
  tmp.bns=nets[[ii]][[4]];
  for(jj in c(1:length(tmp.ind))){
    tmp_ind=length(tmp.ind[[jj]]);
    if(11>tmp_ind&&tmp_ind>2){# Only check larger than 2!
      cat(ii,"~",jj,"~",tmp_ind,"\n")
      set.cc=rbind(set.cc,c(ii,jj,tmp_ind,
                  c_BF_Sensitivity(tmp.bns[[jj]],tmp_ind),
                  c_BF_Effective(tmp.bns[[jj]],tmp_ind) ));
    }
  }
}


nets=BoolGRN_KadelkaSet;
set.ka=NULL;
for(ii in c(1:length(nets))){
  tmp.ind=nets[[ii]][[2]];
  tmp.bns=nets[[ii]][[4]];
  for(jj in c(1:length(tmp.ind))){
    tmp_ind=length(tmp.ind[[jj]]);
    if(11>tmp_ind&&tmp_ind>2){# Only check larger than 2!
      cat(ii,"~",jj,"~",tmp_ind,"\n")
      set.ka=rbind(set.ka,c(ii,jj,tmp_ind,
                            c_BF_Sensitivity(tmp.bns[[jj]],tmp_ind),
                            c_BF_Effective(tmp.bns[[jj]],tmp_ind) ));
    }
  }
}

nets=BoolGRN_ThresholdModel;
set.tr=NULL;
for(ii in c(1:length(nets))){
  tmp.ind=nets[[ii]][[2]];
  tmp.bns=nets[[ii]][[4]];
  for(jj in c(1:length(tmp.ind))){
    tmp_ind=length(tmp.ind[[jj]]);
    if(11>tmp_ind&&tmp_ind>2){# Only check larger than 2!
      cat(ii,"~",jj,"~",tmp_ind,"\n")
      set.tr=rbind(set.tr,c(ii,jj,tmp_ind,
                            c_BF_Sensitivity(tmp.bns[[jj]],tmp_ind),
                            c_BF_Effective(tmp.bns[[jj]],tmp_ind) ));
    }
  }
}


# netinfo=rbind.data.frame(
#   cbind.data.frame(set.cc,t="Cc",c="#2156A6AA"),
#   cbind.data.frame(set.ka,t="Ka",c="#F4CD0BAA"),#"#EE312EAA"
#   cbind.data.frame(set.tr,t="Tr",c="#EE312EAA"))#"#ECB112AA""#39AF8499"
netinfo=rbind.data.frame(
  cbind.data.frame(set.cc,t="Cc",c="#f5b800AA"),
  cbind.data.frame(set.ka,t="Ka",c="#ff7bacAA"),
  cbind.data.frame(set.tr,t="Tr",c="#3db8b8AA"))
#netinfo=netinfo[netinfo[,1]>=5,]
colnames(netinfo)=c("net","node","d","s","e","t","c");


library(ggplot2)
ggplot(netinfo, aes(x = e/d, y = s/d, col=c, shape = t, size = d, fill = c)) +
  geom_point(stroke = 0.25) +   scale_shape_manual(values = c("Cc" = 24, "Ka" = 22, "Tr" = 21)) + 
  scale_size(range = c(2, 5)) +  
  scale_color_manual(values =rep("#000000FF",3))+
  scale_fill_identity(guide = "legend",labels=LETTERS[1:9])+
  coord_fixed(ratio = 0.8)+
  theme_minimal()
# 7*3.5

table(netinfo[,3])
