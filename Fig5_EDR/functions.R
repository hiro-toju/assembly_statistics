######### taxaratio(211224追加)######
# seqtabとtaxonomylist(どっちもmatrix型)を突っ込んだら
# class~genusまでにまとめて返してくれる
# stdは抜いておくこと！！！！！！！！！
################################

######### sum_other(211224追加)######
# seqtab様のmatrix(genusとかasvとかなんでもよい)
# (taxaratioのoutputの一部でいい)、と切りあげたい要素の個数を突っ込んだら
# それ以下のものを"others"にして返す
#######################################

######### rare_merge(220107追加)######
# seqtab様のmatrix(genusとかasvとかなんでもよい)
# いくつでrarefuctionをかけるか、環境条件(実験条件)を引数に
#######################################

######### palettes(220124, 藤田scriptsからコピー)######
#######################################

######### remove.outliers(220414追加)######
# スミルノフ・グラブズ検定(正規分布を仮定)を繰り返し利用した手法
# listで出力，はじいた数値[[2]]とはじかれなかった数値[[1]]
# https://www.trifields.jp/how-to-remove-outliers-using-smirnov%E2%80%90grubbs-test-in-r-2114 からコピー
#######################################

######### forjac(220908追加)##########
# 普通のmatrixを0/1のmatrixに変換する関数
# thresholdsを指定出来て、デフォルト値は1
# dataframeで出力したいならdataframe=Tにする、デフォルトはF
#######################################

taxaratio <- function(seqtab,taxa){
  f <- function(seqtab,taxa,taxid){
    taxa_n <- cbind(ASV=rownames(taxa),ta=taxa[,(colnames(taxa)==taxid)])
    ls <- unique(taxa_n[,2]); ans <- matrix(0,nrow=nrow(seqtab),ncol=length(ls))
    for(i in 1:length(ls)){
      key <- taxa_n[(taxa_n[,2] == ls[i]),1]
      if(length(key)>1){
        ans[,i] <- rowSums(seqtab[,(colnames(seqtab) %in% key)])
      }else{
        ans[,i] <- seqtab[,(colnames(seqtab) %in% key)]
      }
    }
    rownames(ans) <- rownames(seqtab); colnames(ans) <- ls
    return(ans)
  }
  cla <- f(seqtab=seqtab,taxa=taxa,taxid="Class")
  ord <- f(seqtab=seqtab,taxa=taxa,taxid="Order")
  fam <- f(seqtab=seqtab,taxa=taxa,taxid="Family")
  gen <- f(seqtab=seqtab,taxa=taxa,taxid="Genus")
  
  return(list(cla,ord,fam,gen))
}

# numは上位何種を切り出すか
sum_other <- function(mat,num){
  a <- colSums(mat)
  b <- sort(a,T)
  c <- (a>b[num])|(a==b[num])
  out <- matrix(0, nrow=nrow(mat), ncol=(sum(c)+1))
  count <- 1
  outname <- c(rep(NA,sum(c)),"Others")
  for(i in 1:ncol(mat)){
    if(c[i]==T){
      out[,count] <- mat[,i]
      outname[count] <- colnames(mat)[i]
      count <- count+1
    }
    else{
      out[,num+1] <- out[,num+1]+mat[,i]
    }
  }
  rownames(out) <- rownames(mat)
  colnames(out) <- outname
  return(out)
}

rare_merge <- function(mat,csv,rarefuction){
  a <- rrarefy(mat, rarefuction)
  b <- cbind(rownames(a),as.data.frame(a))
  colnames(b) <- c("sample_name",colnames(b))
  c <- merge(b,csv,by="sample_name")
  
  return(c)
}

###############
lib <- 'RColorBrewer';library(package = lib, character.only=TRUE);packageVersion(lib) 
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col.vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

palettes <- function(x){ col.vector[1:length(unique(x))]}


palettes_2 <- function(x){
  unit <- length(x)
  pal <- c("#ff4b00","#804000","#03af7a","#005aff","#4dc4ff",
           "#ff8082","#f6aa00","#990099","#804000")
  ans <- pal[1:unit]; return(ans)
}

###############
remove.outliers <- function(x, conf.level = 0.95)
{
  x <- x[!is.na(x)]
  del.val <- NULL
  
  while (TRUE) {
    n <- length(x)
    if (n < 3) {
      break
    }
    
    r <- range(x)
    t <- abs(r - mean(x)) / sd(x)
    q2 <- (n - 2) / ((n - 1) ^ 2 / t ^ 2 / n - 1)
    q2[q2 < 0] <- 0
    q <- sqrt(q2)
    p <- n * pt(q, n - 2, lower.tail = FALSE)
    
    if (t[1] < t[2]) {
      if (p[2] < 1 - conf.level) {
        del.val <- c(del.val, r[2])
        x <- x[x != r[2]]
        next
      }
    } else {
      if (p[1] < 1 - conf.level) {
        del.val <- c(del.val, r[1])
        x <- x[x != r[1]]
        next
      }
    }
    break
  }
  return(list(x = x, del.val = del.val))
}

forjac <- function(mat,thr=1,dataframe=F){
  ans <- mat
  for(i in 1:nrow(ans)){
    for(j in 1:ncol(ans)){
      if(ans[i,j]>(thr-1)){
        ans[i,j] <- 1
      }else{
        ans[i,j] <- 0
      }
    }
  }
  if(dataframe==T){
    ans <- as.data.frame(ans)
  }
  return(ans)
}
