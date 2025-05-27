#' @title Convert Boolean expressions as mapping tables
#' @description The purpose is to reveal the relations between the mapping 
#' results of Boolean functions and their inputs. Its implementation relies on 
#' regular expressions.
#' @param BoolExpChr a valid Boolean expression, function without judgment capacity.
#' @param TableForm bool, should return table form? (\code{TRUE}, matrix; \code{FALSE}, a bool-vector)
#' @details Note that the function does not assess the validity of the expression.
#' @return a bool-vector or only-show the table.
#' @export
#' @examples
#' # Convert \code{X = ((NOT A AND B)) OR C}.
#' # BoolFun_Expression2MapTable("X = ((NOT A AND B)) OR C")
#' # Return bool-vec c(FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE).
#' boolexp="X = ((NOT A AND B)) OR C"
#' BoolFun_Expression2MapTable(boolexp)
#' 
#' # Show-only \code{Y = (NOT A_a OR NOT B) AND (cc OR dd)}.
#' # BoolFun_Expression2MapTable("Y = (NOT A_a OR NOT B) AND (cc OR dd)", TRUE)
#' # Return a truth table: (0,0,0,0 | 0), (0,0,0,1 | 1), ..., (1,1,1,1 | 0).
#' # Last column is c(0,1,1,1, 0,1,1,1, 0,1,1,1, 0,0,0,0)
#' boolexp="Y = (NOT A_a OR NOT B) AND (cc OR dd)"
#' BoolFun_Expression2MapTable(boolexp,TRUE)
BoolFun_Expression2MapTable<-function(BoolExpChr, TableForm=FALSE){
    exprs=BoolExpChr;
    exprs=gsub("\t", "  ", exprs, perl = TRUE);
    exprs=gsub("/", "_", exprs, perl = TRUE);
    exprs=gsub("\\.", "_", exprs, perl = TRUE);
    exprs=gsub("-", "_", exprs, perl = TRUE);
    exprs=gsub("@", "_", exprs, perl = TRUE);

    tmps=unlist(strsplit(exprs, "="))
    exprs=tmps[2];Outer=tmps[1];Outer=gsub(" ","",Outer);
    exprs=unlist(strsplit(exprs, " "));
    exprs=gsub("(?<![\\w_])(OR|or)(?![\\w_])", "|", exprs, perl = TRUE)
    exprs=gsub("(?<![\\w_])(AND|and)(?![\\w_])", "&", exprs, perl = TRUE)
    exprs=gsub("(?<![\\w_])(NOT|not)(?![\\w_])", "!", exprs, perl = TRUE)
    exprs=paste(exprs,collapse = " ");
    # Analysis dict
    dicts <- strsplit(exprs, "[^[:alnum:]_]");
    dicts <- lapply(dicts, function(x) x[x != ""]);
    dicts <- lapply(dicts, function(x) x[x != " "]);
    dicts=unlist(dicts);
    dicts=unique(dicts);
    n.word=length(dicts);
    n.leng=bitwShiftL(1,n.word);
    ttt=matrix(NA,n.leng,n.word);
    res=rep(NA,n.leng);
    colnames(ttt)=dicts;
    rownames(ttt)=paste0("[",c(1:n.leng)-1,"]");
    for(ii in c(0:(n.leng-1))){
      ttt[ii+1,]=as.numeric(intToBits(ii)[n.word:1]);
      tmpx=paste0(dicts,"=",ttt[ii+1,],";");
      tmpx=paste(tmpx,collapse = " ");
      res[ii+1]=eval(parse(text = paste0(tmpx,exprs)));
    }
    ttt=cbind.data.frame(ttt,Output=as.integer(res));
    colnames(ttt)[n.word+1]=Outer;
  if(TableForm){
    return(ttt);}
  else {
    return(as.logical(ttt[,n.word+1]));}
}

