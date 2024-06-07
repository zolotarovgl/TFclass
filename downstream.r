# pick the best hit per query 
library(stringr)
library(dplyr)
library(data.table)
list.files()
infile = 'results/Octsin.HG.best'
options(max.print =30)
inf = read.table(infile)
rownames(inf) = inf[,3]
length(unique(inf[,1]))
length(unique(inf[,3]))
# predicted more Homeodomains than there are!
#################### HG level ###########################
# 104 queries against 37 models. 
# compare to the real classifications 
q2t.best = setNames(inf[,1],inf[,3])
hg_q2t.best = q2t.best
length(q2t.best)
truth = as.data.frame(fread('data/classifications/tfs.Octsin_genes.curated.csv',fill = T)[,1:2])
truth = setNames(truth[,2],truth[,1])
truth = truth[names(truth) %in% names(q2t.best)]
truth_hg = setNames(sapply(str_split(truth,'\\.'),FUN = function(x) paste0(x[2],'.',x[3])),names(truth))
truth_og = setNames(sapply(str_split(gsub('tfs.','',truth),':'),FUN = function(x) x[1]),names(truth))

# FALSE positives
FP = names(q2t.best[!names(q2t.best) %in% names(truth_hg)])

summary(inf[FP,]$V5)
summary(inf[!rownames(inf) %in% FP,]$V5)
# ok, they clearly differ from their counterparts. 

writeLines(FP,'overpred')
library(ComplexHeatmap)
library(reshape2)
library(tidyr)
library(fossil)

df = data.frame(pep = names(truth_hg),truth = truth_hg,predicted = q2t.best[names(truth_hg)])
df$truth = factor(df$truth,levels = unique(df$truth))
df$predicted = factor(df$predicted,levels = unique(df$truth))
adj.rand.index(as.integer(df$truth),as.integer(df$predicted))
# not perfect but not horrible either. 
# hmmm
# we can also get informed by the stability of the 

# Cross-class classification
################### OG-level classifications ##################################
# Are the best OG classifications concordant with the best HG classifications? 

infile = 'results/Octsin.OG.best'
options(max.print =30)
inf = read.table(infile)
rownames(inf) = inf[,3]
length(unique(inf[,1]))
length(unique(inf[,3]))
# predicted more Homeodomains than there are!
#################### HG level ###########################
# 104 queries against 37 models. 
# compare to the real classifications 
q2t.best = setNames(inf[,1],inf[,3])
og_q2t.best = q2t.best


df = data.frame(pep = names(truth_og),truth = truth_og,predicted = q2t.best[names(truth_og)])
df$truth = factor(df$truth,levels = unique(df$truth))
df$predicted = factor(df$predicted,levels = unique(df$truth))
adj.rand.index(as.integer(df$truth),as.integer(df$predicted))
# 97% adjusted rand index - woooow 

# are they consistent? 

########################################################
# OG-HG consistency 
int = unique(names(og_q2t.best),names(hg_q2t.best))
d = data.frame(gene = int, hg_pred = hg_q2t.best[int],og_pred = og_q2t.best[int])
d$og_pred_hg = sapply(str_split(d$og_pred,'\\.'),FUN = function(x) paste0(x[1:2],collapse = '.'))
d$hg_pred = factor(d$hg_pred)
d$og_pred_hg = factor(d$hg_pred,levels = levels(d$hg_pred))
adj.rand.index(as.integer(d$hg_pred),as.integer(d$og_pred_hg))
# 1?? 
# table(d$hg_pred==d$og_pred_hg)
# in all cases, the picked hg is fine!
# now, what about OG-level classifications? 
# we can control precision by looking at the delta between the best hits 



infile = 'results/Octsin.OG'
options(max.print =30)
inf = read.table(infile)
inf = inf[inf$V3 %in% names(truth),]
# there are some missing proteins? 
# why? 
# aha! C2H2s? 

# select top n hits per query
f = inf%>%group_by(V3)%>%arrange(V3,V5)%>%top_n(2,-V5)
# what do the fields here mean??? 
# you need to format your inputs / outputs names
f = lapply(split(f,f$V3),FUN = function(x) sort(setNames(x$V5,x$V1)))
f = f[sapply(f,length)>1]
# so, for the cases, where there are more than one good hits 
log10_enrich = sapply(f,FUN = function(x) log10(x[2])-log10(x[1]))
names(log10_enrich) = names(f)
plot(density(log10_enrich))
# so we can cut the E-value here. 
# I mean, based on this threshold, we can decide how precise we want to be in our classifications. 
# Check proteins with only one hit:
# true positive rate 
# false positive rate 
# we can see the ROC
# this is a multi-class classification problem actually. 
#message(sprintf('%s queries with single OG',length(single_best)))
# how well does it agree with the ground truth? 
# the One-vs-Rest scheme compares each class against all the others (assumed as one);
thr = 1
table(log10_enrich>=thr)
message(sprintf('%s/%s (%s%%) predictions have a top E-value hit with log10 difference above %s',sum(log10_enrich>=thr),length(log10_enrich),round(sum(log10_enrich>=thr)/length(log10_enrich),1)*100,thr))
# so, you of course can output the best classifications like that. 
