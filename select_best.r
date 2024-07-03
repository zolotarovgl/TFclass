# Select the best classifications per query 
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
options(max.print = 30)

#infile_hg = 'results/Octsin.HG'
#infile_og = 'results/Octsin.OG'
#outfile = 'classification.tsv'


args <- commandArgs(trailingOnly = TRUE)
#print(args)

if(length(args)<3){
    message('script.R infile_hg infile_og outfile')
    quit()
}else{
    infile_hg = args[1]
    infile_og = args[2]
    outfile = args[3]
}

thr = 1 # log10 E-value difference between the best hit and the second best hit to report. 

select_best = function(infile,thr = 1){
    # specify the columns!!!
    inf = read.table(infile)
    message(sprintf('%s queries; %s targets',length(unique(inf$V3)),length(unique(inf$V1))))
    # select top n hits per query
    f = inf%>%group_by(V3)%>%arrange(V3,V5)%>%top_n(2,-V5)
    f = lapply(split(f,f$V3),FUN = function(x) sort(setNames(x$V5,x$V1)))
    f = f[sapply(f,length)>1]
    log10_enrich = sapply(f,FUN = function(x) log10(x[2])-log10(x[1]))
    names(log10_enrich) = names(f)
    message(sprintf('%s/%s (%s%%) predictions have a top E-value hit with log10 difference above %s',sum(log10_enrich>=thr),length(log10_enrich),round(sum(log10_enrich>=thr)/length(log10_enrich),1)*100,thr))
    o = do.call(rbind,lapply(seq_along(f),FUN = function(i) c(names(f)[i],names(f[[i]])[1],as.numeric(log10_enrich[i]))))
    o = as.data.frame(o)
    o$V3 = as.numeric(o$V3)
    colnames(o) = c('ID','Predicted','Log10_first_second')
    # so, you of course can output the best classifications like that. 
    return(o)
}


hg = select_best(infile_hg)
og = select_best(infile_og)
u = left_join(hg,og,by ='ID')
message(sprintf('Missing OG predictions for %s proteins!',sum(is.na(u$Predicted.y))))
# on average, we are more sure about classifying into HGs - as expected. 
# how about the representative orthologs from other species? 
# we need to build a simple classifier from this. 
#sort(table(sapply(str_split(u$Predicted.x,'\\.'),FUN = function(x) x[1])))
colnames(u) = c('ID','Predicted.HG','Log10.HG','Predicted.OG','Log10.OG')
u$Log10.HG[is.na(u$Log10.HG)] = 0
u$Log10.OG[is.na(u$Log10.OG)] = 0

message(sprintf('%s HG assignments exceeding threshold',sum(u$Log10.HG>=thr)))
message(sprintf('%s OG assignments exceeding threshold',sum(u$Log10.OG>=thr)))
if(any(u$Log10.HG<thr)){
	u[u$Log10.HG<thr,]$Predicted.HG = ''
}
if(any(u$Log10.OG<thr)){
	u[u$Log10.OG<thr,]$Predicted.OG = ''
}
u$Log10.HG = round(u$Log10.HG,2)
u$Log10.OG = round(u$Log10.OG,2)
write.table(u,outfile,sep = '\t',quote = F,row.names = F)
