#' @title Recover a multi-valued function from multi-valued model
#' @description Unlike the standard multi-valued logical paradigm, 
#' we adopt the definition of so-called qualitative networks; please refer to 
#' [\href{https://doi.org/10.1186/1752-0509-1-4}{Schaub2007}]
#' @param aMulVNet A four-element formatted list that records information of a multi-valued
#' @param GeneName A string that denotes the gene name.
#' @param MaxValue An integer that represents the max value of the gene (Default: 2).
#' @param MapTab A logical value. Is the mapping table also returned? 
#' (Default: \code{FALSE}) 
#' @return An IntegerVector or a mapping table.
#' @export 
#' 
MulVFun_Recovery<-function(aMulVNet, GeneName, MaxValue=2L, MapTab=FALSE){
  final=NULL;
  if(!(GeneName%in%aMulVNet[[1]])){
    stop("GeneName not in the Network. Please check it.\n");
  } else {
    if(is.na(aMulVNet$InEdge[[GeneName]][1])){# No input: 
      final=c(0:MaxValue);
      if(MapTab){
        final=cbind(c(0:MaxValue),final);
        colnames(final)=c(GeneName,paste0("f_",GeneName));
      }
    } else {# Have input:
      innode=aMulVNet$AllMember[aMulVNet$InEdge[[GeneName]]+1];
      if((GeneName%in%innode)){# Exist the self-loop
        idx=which(GeneName==innode);
        xxx=aMulVNet$BoolFun[[GeneName]];
        ttt=r_GenerateTruthTable_M(length(xxx), (MaxValue+1));
        index=ttt%*%c(xxx);
        is.all.n=all(xxx<0);
        if(is.all.n){# all negative
          index=NumLev+index;
          final=ttt[,idx]+sign(index-ttt[,idx]);
        } else {# has least one positive
          index[index<0]=0;
          final=ttt[,idx]+sign(index-ttt[,idx]);
        }
        final[final>MaxValue]=MaxValue;
        final[final<0]=0;
        if(MapTab){
          final=cbind(ttt,final);
          colnames(final)=c(innode,paste0("f_",GeneName));
        }
      } else {# GeneName not exist the self-loop
        xxx=aMulVNet$BoolFun[[GeneName]];
        ttt=r_GenerateTruthTable_M(length(xxx)+1, (MaxValue+1));# Highest is itself
        index=ttt%*%c(0,xxx);
        is.all.n=all(xxx<0);
        if(is.all.n){# all negative
          index=MaxValue+index;
          final=ttt[,1]+sign(index-ttt[,1]);
        } else {# has least one positive
          index[index<0]=0;
          final=ttt[,1]+sign(index-ttt[,1]);
        }
        final[final>MaxValue]=MaxValue;
        final[final<0]=0;
        if(MapTab){
          final=cbind(ttt,final);
          colnames(final)=c(GeneName,innode,paste0("f_",GeneName));
        }
      }
    }
  }
  return(final);
}