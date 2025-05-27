#' @title Show the percolation
#' @description This function displays the state distribution at each point under 
#' square/triangle/hexagon lattice. Nodes are generally classified into three 
#' categories: stable-1, stable-0, and oscillatory. This function depends on the 
#' frame of \code{ggplot} class.
#' @param Values integer vector, a valid file/address of target GenNet.
#' @param Length an integer, Length of lattice matrix.
#' @param Height an integer, Height of lattice matrix.
#' @param Type an integer, 4, Square; 3, triangle; 6, hexagon.
#' @param ColorFills character vector [3], [1] unstable, [2]0-stable, [3]1-stable
#' @param LineColor a character, color of polygonal boundary. (\code{NA}, no diaply)
#' @param Linewidth numeric value, line width of polygonal boundary.
#' @return a \code{ggplot} class, show the percolation.
#' @export 
#' 
FigLatticePlot<-function(Values=NULL, Length=NULL, Height=NULL, Type=c(4L,3L,6L),
    ColorFills=c("#efefef","#EE312E","#2156A6"),
    LineColor="#111111", Linewidth=0.7){
    if(!as.integer(Type)[1]%in%c(3L,4L,6L)){
        stop("Invalid geometrical parameter.\n");}
    if(4L==as.integer(Type[1])){# Square (4)
        tmpDig=ceiling(log10(Length*Height))+1;# Add one to avoid some digits error.
        values=data.frame(value=Values,ids=factor(
        paste("N",formatC(c(1:(Length*Height)),width=tmpDig,flag="0"),sep="")));
        nnn=Length;mmm=Height;
        Xs=seq(0,nnn-1,1);
        Xs1=as.numeric(rbind(Xs,Xs+1));
        Xs2=as.numeric(rbind(Xs+1,Xs));
        Xs=c(Xs1,Xs2);
        Xs=rep(Xs,time=mmm);
        Ys=rep(0,nnn*2);Ys=c(Ys,Ys+1);
        Ys=rep(Ys,time=mmm)+rep(c(0:(mmm-1)),each=4*nnn);
        id0=rep(rep(c(1:nnn),each=2),time=2);
        ids=id0+nnn*rep(c(0:(mmm-1)),each=length(id0));
        ids=paste("N",formatC(as.integer(ids),width=tmpDig,flag="0"),sep="");
        id_labels=data.frame(ids=factor(ids),x=Xs,y=Ys);
        id_labels=id_labels[order(id_labels[,1]),];
        datapoly=merge(values,id_labels,by=c("ids"));
    } else if(6L==as.integer(Type[1])){# Hexagon (6)
        tmpDig=ceiling(log10(Length*Height));
        values=data.frame(value=Values,ids=factor(
        paste("N",formatC(c(1:(Length*Height)),width=tmpDig,flag="0"),sep="")));
        nnn=Length;mmm=Height;
        idx1=seq(1,2*nnn-1,2);idx2=idx1+1;
        idx=rep(c(idx1,idx2),time=mmm/2);
        idx=idx*sqrt(3)*0.5;
        idy=rep(seq(1,0.5+1.5*mmm,1.5),each=nnn);
        x_mod=c(0,1,1,0,-1,-1)*sqrt(3)*0.5;x_mod=rep(x_mod,time=length(idx));
        y_mod=c(-1,-0.5,0.5,1,0.5,-0.5);y_mod=rep(y_mod,time=length(idy));
        Xs=rep(idx,each=6)+x_mod;
        Ys=rep(idy,each=6)+y_mod;
        ids=rep(c(1:(nnn*mmm)),each=6);
        ids=paste("N",formatC(ids,width=tmpDig,flag="0"),sep="");
        id_labels=data.frame(ids=factor(ids),x=Xs,y=Ys);
        datapoly=merge(values,id_labels,by=c("ids"));
    } else {# Triangle (3)
        tmpDig=ceiling(log10(Length*Height));
        values=data.frame(value=Values,ids=factor(
        paste("N",formatC(c(1:(Length*Height)),width=tmpDig,flag="0"),sep="")));
        nnn=Length;mmm=Height;
        Xs1=c(rep(seq(1,nnn-1,2),each=3),nnn+1);Xs1=Xs1[-1];
        Xs2=c(0,rep(seq(2,nnn,2),each=3));Xs2=Xs2[-length(Xs2)];
        Xs=c(Xs1,Xs2,Xs2,Xs1);
        Xs=rep(Xs,time=mmm/2);
        Ys=c(1:mmm);Ys=as.numeric(rbind(Ys-1,Ys));
        Ys=rep(sqrt(3)*Ys,each=(nnn/2*3));
        id1=seq(2,nnn,2);id1=as.numeric(rbind(id1-1,id1,id1));
        id2=seq(2,nnn,2);id2=as.numeric(rbind(id2-1,id2-1,id2));
        id0=c(id1,id2,id2+nnn,id1+nnn);
        ids=rep(id0,time=mmm/2)+2*Length*rep(c(0:(mmm/2-1)),each=length(id0));
        ids=paste("N",formatC(ids,width=tmpDig,flag="0"),sep="");
        id_labels=data.frame(ids=factor(ids),x=Xs,y=Ys);
        datapoly=merge(values,id_labels,by=c("ids"));
    }
    tmpp=ggplot(data=datapoly,aes(x=x,y=y)) +
    #tmpp=ggplot(data=datapoly, aes(x=.data[[x]],y=.data[[y]])) +
        geom_polygon(aes(fill=value, group=ids), colour=LineColor, linewidth=Linewidth)+
        scale_fill_manual(values=ColorFills);
    if(4L==as.integer(Type[1])){
        tmpp=tmpp+scale_x_continuous(position="bottom",expand=c(0,0), limits=c(0,Length))+
            scale_y_continuous(position="left",expand=c(0,0), limits=c(0,Height));
    } else if(6L==as.integer(Type[1])){
        tmpp=tmpp+scale_x_continuous(position="bottom",expand=c(0,0),limits=c(0,(Length+0.5)*sqrt(3)+0.1))+
            scale_y_continuous(position="left",expand=c(0,0),limits=c(0,Height*1.5+0.5));
    } else {
        tmpp=tmpp+scale_x_continuous(position="bottom",expand=c(0,0), limits=c(0,Length+1))+
            scale_y_continuous(position="left",expand=c(0,0), limits=c(0,Height*sqrt(3)));
    }
    tmpp=tmpp+guides(fill="none")+
        theme(
            axis.line=element_blank(),
            axis.ticks.x.top=element_blank(),
            axis.ticks.x.bottom=element_blank(),
            axis.ticks.y.left=element_blank(),
            axis.ticks.y.right=element_blank(),
            panel.grid=element_blank(), # Delete grid
            panel.background=element_blank(), # Delete bakcground
            legend.key=element_rect(fill="#FFFFFF"),
            legend.title=element_text(face="bold",size=14),
            legend.text=element_text(face="bold",size=12),
            plot.title=element_text(hjust=0.5),
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            axis.title.x.bottom=element_blank(),
            axis.title.y.left=element_blank() );
    return(tmpp);
}