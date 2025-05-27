library(viridis)
library(ggtern)
# Ternary plot for Activity & Sensitivity (K=2) 
TernaryContour3<-function(SimDatas,whichIDs=NULL,Legends=NULL,inc=0.10,bins=0.1/3,VValues=c(0,1,2),CColors=c("red","white","blue"),
    text=FALSE,Utilize=c("aa","bb","cc")){
    count=1;Decimals=2;utilized=c("IDLabel",Utilize);
    # Each type of points
    points=data.frame();
    for(z in seq(0,1.0,inc)){
        x=1-z;y=0;
        while(x>0){
            points=rbind(points,c(count,x,y,z));
            x=round(x-inc,digits=Decimals);y=round(y+inc,digits=Decimals);
            count=count+1;}
        points=rbind(points,c(count,x,y,z));
        count=count+1;}
    colnames(points)=c("IDPoint","T","L","R");
    # Self-defined all triangle triwise node's coordinate info.
    polygons=data.frame();count=1;
    # Normal triangles T-L-R
    for(p in points$IDPoint){
        if(p%in%points$IDPoint[points$T!=0]){# Pickout the non-zeros of T-cols
            pT=points$T[points$IDPoint==p];
            pL=points$L[points$IDPoint==p];
            pR=points$R[points$IDPoint==p];
            polygons=rbind(polygons,c(count,p),
                c(count,points$IDPoint[abs(points$L-pL)<inc/2 & abs(points$R-pR-inc)<inc/2]),
                c(count,points$IDPoint[abs(points$L-pL-inc)<inc/2 & abs(points$R-pR)<inc/2]));
                count=count+1;}}
    # Upside down triangles L-R-D(T)
    for(p in points$IDPoint){
        if(!is.element(p,points$IDPoint[points$T==0])){
            if(!is.element(p,points$IDPoint[points$L==0])){
                pT=points$T[points$IDPoint==p];
                pL=points$L[points$IDPoint==p];
                pR=points$R[points$IDPoint==p];
                polygons=rbind(polygons,c(count,p),
                    c(count,points$IDPoint[abs(points$T-pT)<inc/2 & abs(points$R-pR-inc)<inc/2]),
                    c(count,points$IDPoint[abs(points$L-pL)<inc/2 & abs(points$R-pR-inc)<inc/2]));
                count=count+1;}}}
    # Important for correct ordering.
    polygons$PointOrder=1:nrow(polygons);# Samll triangles (nrow/3): [1+(2n-1)]*n/2 -> n=1/inc
    colnames(polygons)=c("IDLabel","IDPoint","PointOrder");
    df_tr=merge(polygons,points);# Integrate all triangle informations
    # The label's tri-coordinates:
    Labs=plyr::ddply(df_tr,"IDLabel",function(x){c(c(mean(x$T),mean(x$L),mean(x$R)))})
    colnames(Labs)=c("IDLabel","T","L","R");
    # Combine simulation datas.
    loaddata=data.frame();logis=NA;
    discr=round(Labs/bins);discr[,1]=Labs[,1];
    for(i in c(1:nrow(discr))){
        logis=(SimDatas[,1]==(discr$T)[i])&(SimDatas[,2]==(discr$L)[i])&(SimDatas[,3]==(discr$R)[i])
        loaddata=rbind(loaddata,c(discr[i,],SimDatas[logis,]));}#discretization
    FinalDatas=plyr::join(df_tr,loaddata[,utilized],by="IDLabel")
    #Labs=plyr::ddply(df,.(IDLabel,N),function(x){c(c(mean(x$T),mean(x$L),mean(x$R)))})
    #colnames(Labs)=c("Label","N","T","L","R")
    heat=ggtern();
    commds=paste0("heat=heat+geom_polygon(data=FinalDatas,aes(L,T,R,fill=",whichIDs,
            ",group=IDLabel),color='#FFFFFF',alpha=1,size=0,linetype=0)")
    eval(parse(text=commds));
    #heat=heat+geom_path(data=Paths,aes(x,y,z),color="black",size=1.20,linetype="dashed")
    Breaks=c(1/6,1/3,0.5,2/3,5/6);Labels=c("1/6","1/3","1/2","2/3","5/6");
    heat=heat+scale_T_continuous(limits=c(0,1.0),breaks=Breaks,labels=Labels)+
        scale_L_continuous(limits=c(0,1.0),breaks=Breaks,labels=Labels)+
        scale_R_continuous(limits=c(0,1.0),breaks=Breaks,labels=Labels);
    heat<-heat+Tlab(latex2exp::TeX("\\mathbf{P_0}"))+Llab(latex2exp::TeX("\\mathbf{P_1}"))+Rlab(latex2exp::TeX("\\mathbf{P_2}"));
    #labs(xarrow="P_0",yarrow="P_1",zarrow="P_2") + 
    if(length(CColors)==2){
        heat<-heat+scale_fill_gradientn(name=Legends,colours=ColorBars);}
    else if(length(CColors)==3){
        heat<-heat+scale_fill_gradient2(
            low=CColors[1],high=CColors[3],mid=CColors[2],
            midpoint = VValues[2]);}
    else if(length(CColors)==1&&is.na(CColors)){
      heat<-heat+scale_fill_viridis_c(
        #begin=0.2,end=0.85,option="G")
      begin=0.05,end=0.95,option="rocket")
    }
    else {
        cat("Error!\n");}
    heat<-heat+theme(
        axis.text=element_text(size=12,face="bold",color="#000000"),
        axis.line=element_line(linewidth =12,linetype=1,color="#000000"),
        axis.ticks=element_line(size=1.5,linetype=1,color="#000000"),
        axis.title=element_text(size=18,face="bold",color="#000000"),
        legend.title=element_text(size=12,face="bold",color="#000000"),
        legend.text=element_text(size=9,face="bold",color="#000000"),
        legend.key.height=unit(1,"cm")
    )
    return(heat);
}

shuju=read.table("./simulation_analysis/engaged_k2.txt",sep=",");
shuju=as.matrix(shuju);
shuju[,6]=log10(1+shuju[,6]);
colnames(shuju)=c("P","Q","M","s","u","re");
shuju=as.data.frame(shuju);
xxx=seq(0,0.666,0.001);tmp=sqrt(2*xxx-3*xxx*xxx);
tmp_k2=cbind(c(xxx,xxx),0.5*c(tmp-xxx+1,-tmp-xxx+1));
tmp_k2=cbind(tmp_k2,1-rowSums(tmp_k2));
colnames(tmp_k2)=c('x','y','z');tmp=nrow(tmp_k2);
tmp_k2=rbind.data.frame(tmp_k2[c(1:(tmp/2)),],tmp_k2[c(tmp:(tmp/2+1)),]);
ppt_k2=TernaryContour3(shuju,"re","Per.Ratio",inc=0.01,bins=0.01/3.0,
                               #CColors = c("#2156A6","white","#EE312E"),
                       CColors = NA,#c("#2156A6","white","#EE312E"),
                               VValues = c(min(shuju[,6]),0.5*(min(shuju[,6])+max(shuju[,6])),max(shuju[,6])),
                               Utilize=colnames(shuju)[3:6]);
ppt_k2_c=ppt_k2+geom_path(data=tmp_k2,aes(x,y,z),color="#3E9EFF",linewidth=1.00,linetype="solid") # inch 6*7
# #FEFE00

# Ternary plot for Activity & Sensitivity (K=3)
shuju=read.table("./simulation_analysis/engaged_k3.txt",sep=",");
shuju=as.matrix(shuju);
shuju[,6]=log10(1+shuju[,6]);
colnames(shuju)=c("P","Q","M","s","u","re");
shuju=as.data.frame(shuju);
xxx=seq(0,1,0.001);tmp=sqrt(3)*sqrt(1+6*xxx-9*xxx*xxx);
tmp_k3=cbind(c(xxx,xxx),c(tmp-3*xxx+3,-tmp-3*xxx+3)/6.00);
tmp=cbind(tmp_k3,1-rowSums(tmp_k3));
tmp=tmp[!(is.nan(tmp[,2])|tmp[,2]<0|tmp[,3]<0),];
tmp_k3=rbind(tmp,tmp[,c(2,3,1)],tmp[,c(3,1,2)]);
colnames(tmp_k3)=c('x','y','z');
tmp_k3=as.data.frame(tmp_k3);
tmp_k3=tmp_k3[order(atan((tmp_k3[,2]-0.5)/(tmp_k3[,1]-0.5))),]
ppt_k3=TernaryContour3(shuju,"re","Per.Ratio",inc=0.01,bins=0.01/3.0,
                              # CColors = c("#2156A6","white","#EE312E"),#
                          CColors = NA,#c("#
                               VValues =c(min(shuju[,6]),0.5*(min(shuju[,6])+max(shuju[,6])),max(shuju[,6])),
                               Utilize=colnames(shuju)[3:6]);
ids=seq(0,nrow(tmp_k3),nrow(tmp_k3)/3);
ppt_k3;
ppt_k3_c=ppt_k3+# #FEDA2E # EE4444
  geom_path(data=tmp_k3[(ids[1]+1):ids[2],],aes(x,y,z),color="#3E9EFF",size=1.25,linetype="solid")+
  geom_path(data=tmp_k3[(ids[2]+1):ids[3],],aes(x,y,z),color="#3E9EFF",size=1.25,linetype="solid")+
  geom_path(data=tmp_k3[(ids[3]+1):ids[4],],aes(x,y,z),color="#3E9EFF",size=1.25,linetype="solid");

ggsave("fig2e2.pdf",ppt_k2_c,width = 5.2,height = 3.5,units = "in");
ggsave("fig2e3.pdf",ppt_k3_c,width = 5.2,height = 3.5,units = "in");

# 6*7 inch HL-3.5*5.2