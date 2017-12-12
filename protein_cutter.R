# By default, Shiny limits file uploads to 5MB per file. You can modify this limit by using the shiny.
# maxRequestSize option. For example, adding options(shiny.maxRequestSize=30*1024^2) to the top of server.R
# would increase the limit to 30MB.
library(seqinr)
library(Biostrings)
library(foreach)
library(msa)
library(reshape2)
library(bios2mds)
library(stringr)
#import sequences of intererst
#pre_aln <- readAAStringSet("C:/Users/mjohn/Documents/R/KRP_unaligned_B_C.fasta")
pre_aln <- readAAMultipleAlignment("C:/Users/mjohn/Documents/R/KRP_hmmer.fa")
aln <- msa(pre_aln)
#aln <- pre_aln
print(aln, show="all")
writeXStringSet(unmasked(aln), "C:/Users/mjohn/Documents/R/alignmentfromR.fa")
data("BLOSUM62")
cc<-which(msaConservationScore(aln, BLOSUM62) > 100)
cm<-consensusMatrix(pre_aln, as.prob=TRUE)
cm<-cm[-1, ]
cc<-names(cm[apply(cm,2,function(x) which.max(x)),1])
aln_m<-as.matrix(aln)
aln_cc<-aln_m[,cc]
colnames(aln_cc) <- cc
aln_uniq<-apply(aln_cc, 2, function(col)prop.table(table(col)))
#create scoring master and reshape to correct format
master<-melt(aln_uniq)
colnames(master) <- c(col="aa",value="conc",L1="pos")
hold<-as.matrix(msaConservationScore(aln, BLOSUM62)) #append quality score
hold<-hold[cc]
#list(qual=hold)
#table(master$pos)
tep <- foreach (i = 1:length(unique(master$pos))) %do%
  {newcol<-rep(hold[i],table(master$pos)[i])}
tep2<-melt(tep)
master$qual<-tep2$value
master$aa <- as.character(master$aa)
master$conc<- as.numeric(master$conc)
master$pos<- as.numeric(master$pos)
master$qual<- as.numeric(master$qual)
#reference genome import
s <- readAAStringSet("C:/Users/mjohn/Documents/R/KRP_hmmer.fa")
seqID <- names(s) #gets names
seq <- paste0(s) #gets sequences
seq <- gsub("-*", "", seq)
fa <- cbind(seqID,seq)
#amino acid chemistry
neg_charge <- c("DE")
hydrophobic <- c("AILMFWVP")
pos_charge <- c("KRH")
polar <- c("NQSTGCY")
neg_charge_logical <- grepl(paste0("[",neg_charge,"]"), master$aa)
hydrophobic_logical <- grepl(paste0("[",hydrophobic,"]"), master$aa)
pos_charge_logical <- grepl(paste0("[",pos_charge,"]"), master$aa)
polar_logical <- grepl(paste0("[",polar,"]"), master$aa)
master$chem <- ifelse(neg_charge_logical,"N",(ifelse(hydrophobic_logical,"H",(ifelse(pos_charge_logical,"Pos",(ifelse(polar_logical,"Pol","Gap")))))))
#score function
#score <- function(sequences, x){#each sequence
  #mn<-4 #softcodeme
  #pattern <- paste0("^.{",x$pos,"}",x$aa)
 # foreach (i = 1:nrow(master)) %do% #each master row
#  {foreach (ii = 1:mn) %do% {#each motif cluster
    #if ((x[i,6]==ii)){
      #ii = motif number (4)
      #i = each row (150)
      #return(x$pos)
  #  }
   #}
  #}
    #if(grepl(pattern[i], seq_list)){
  #    obs <- x$conc[i]*x$qual[i]
   #   expected <- x$chem[i]*x$qual[i]
  #    v <- (obs-expected)/obs
  #  }
   # else{v <- 0}
  #  return(v)}
  #}
#lapply(seq_list[[1]], score, x=master)
#cluster motifs
library(ggplot2)
library(dbscan)
p <- ggplot(master, aes(x=pos,y=conc,color=chem)) +
  geom_point()+
  facet_wrap(~chem)

scan <- data.frame(x=master$pos,y= master$conc)
cluster.dbscan <- dbscan(scan, eps = 7, minPts = 5, borderPoints = T)
plot(y ~ x, data = scan, col = cluster.dbscan$cluster + 1L, pch = 20)
master$cluster<-cluster.dbscan$cluster

db_c<-data.frame(x=1:928,y=apply(cm, 2, function(x) {ifelse(max(x>0.1), max(x), NA)}))[complete.cases(data.frame(x=1:928,y=apply(cm, 2, function(x) {ifelse(max(x>0.1), max(x), NA)}))), ]
cluster.dbscan <- dbscan(db_c, eps = 7, minPts = 5, borderPoints = T)
plot(y ~ x, data = db_c, col = cluster.dbscan$cluster + 1L, pch = 20)


#generates regular expressions
pattern2 <- list()
for (i in unique(master$cluster)){
  current_cluster <- master[master$cluster == i, ]
  pos_list <- unique(current_cluster$pos)
  n_col <- length(pos_list)
  name <- paste('cluster:',i,sep='')
  for (ii in 1:n_col){
    expression <- str_c(unlist(current_cluster[current_cluster$pos==pos_list[ii], ]$aa), collapse='')
    dist <- diff(pos_list,1)[ii]
    expression_with_distance <- paste("[",expression,"]", ".{", 0, ",", dist, "}", sep='')
    pattern2[[name]][ii] <- expression_with_distance
  }
}
pattern3<-lapply(pattern2, function(x) str_c(x,collapse=''))
pattern4<-lapply(pattern3, function(x) gsub("\\.\\{0,NA\\}", "", x))
pattern5<-lapply(pattern4, function(x) gsub("-", "", x))
#info of where matches are
info <-list()
for (i in 1:5){#soft code
name <- paste('cluster_info:',i,sep='')
extract_motif <- function(fasta, grep_lst){#ITERATING THROUGH GENOME
 regexpr(grep_lst[i], fasta)
 }
info[[name]]<-lapply(fa[,2], extract_motif, grep_lst = pattern5)
}

block<-lapply(info, function(x) regmatches(fa[,2],x))

#extract <- function(fasta, anchor_info){
#  locations <- which(anchor_info != -1)
#  begin <- as.numeric(anchor_info[locations])
#  location_of_size <- anchor_info + attr(anchor_info, "match.length")
#  size <- location_of_size[locations]-1
#  hits.seq <- substr(fasta[locations], begin, size)
#  return(hits.seq)
#}
#seq_list <- lapply(anchor_info, extract, fasta = fa_v2$seq)

