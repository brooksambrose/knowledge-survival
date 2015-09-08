#Port fuzzy set functions to here
if(F){if(F){
cr2zcr.f<-function( #utility for resolving identity uncertainty for WOK CR field
	cr=stop("cr = supply a character vector of unique and uncertain CR ids") #unprocessed bimodal edgelist in UT CR order
	,just.normalize=T # TRUE will remove DOI and capitalize. FALSE will perform fuzzy set replacement.
)
{
require(stringdist)
require(randomForest)
load("/Users/bambrose/Dropbox/2014-2015/Sprints/1/BOURDIEU, 1985, THEOR SOC/out/wok2dbl.RData")
if(!is.data.table(wok2dbl)) wok2dbl<-data.table(wok2dbl)
setkey(wok2dbl,field,id)

### impose formatting and nomenclature ###
rawbel<-wok2dbl["CR",c(1,3),with=F] #delete
#save.rawbel<-copy(rawbel)
setnames(rawbel,1:2,c("ut","cr"))
rawbel[,`:=`(
	cr=as.character(cr)
	,ut=as.character(ut)
)]
(rawbel[,cr:=sub(", DOI .+","",cr)]) #remove DOI
(rawbel[,cr:=gsub("(\\w)","\\U\\1",rawbel$cr,perl=T)]) #capitalize
(rawbel[,cr:=sub("\\[ANONYMOUS\\], ","",cr)]) #remove anonymous
setkey(rawbel,ut,cr) #sort

#if(!just.normalize){ #fuzzy set replacement


	### build database to define features of pairwise string comparison model ###

	## first define new database of citations, including original degree info
	setkey(rawbel,cr)
	cr<-rawbel[,.N,by=cr]

	## greedy comparitor and pick threshhold

	#### ok  match concept is good, but using amatch throws away information that we'll need later

bmnl.f<-function(jcr,jw.thresh=.1,jw.penalty=.1){ # records best match nodelist (bmnl); nodelist of best match(es) only
		jcr<-as.character(jcr)
		require(stringdist)
		require(data.table)
		x<-list()
		cat("\n")
		l<-length(jcr)
		for(i in 1:l){
			cat("\r",round(i/l,3),"\t")
			x[[length(x)+1]]<-stringdist(jcr[i],jcr[-i],method="jw",p=jw.penalty) # compute jw distance
			if(!any(x[[length(x)]]<=jw.thresh)){x[[length(x)]]<-NULL;next} # if none pass threshold, clear and move on to the next
			w<-which(x[[length(x)]]==min(x[[length(x)]]))
			x[[length(x)]]<-data.table("Target"=jcr[-i][w],"jw"=x[[length(x)]][w]) # data table of node and distance as edge weight
			rm(w)
			names(x)[length(x)]<-jcr[i]
		}
		x
	}

	library(microbenchmark)
	(mb<-microbenchmark(bmnl_1<-bmnl.f(cr$cr),bmnl_2<-bmnl.f(cr$cr,jw.thresh=.2),times=1))

	bm.el1<-data.table("Source"=factor(rep(names(bmnl_1),sapply(bmnl_1,function(x) dim(x)[1]))),rbindlist(bmnl_1))
	setnames(bm.el1,2:3,c("Target","jw"))
	setkey(bm.el1,Source,Target)
	(bm.el1[,cjw:=1-jw])

	bm.el2<-data.table("Source"=factor(rep(names(bmnl_2),sapply(bmnl_2,function(x) dim(x)[1]))),rbindlist(bmnl_2))
	setnames(bm.el2,2:3,c("Target","jw"))
	setkey(bm.el2,Source,Target)
	(bm.el2[,cjw:=1-jw])

	setwd("/Users/bambrose/Dropbox/2014-2015/Sprints/1/BOURDIEU, 1985, THEOR SOC/")
	write.table(bm.el1,file="out/bm.el1.tab",quote=F,sep="\t",row.names=F,col.names=T)
	write.table(bm.el2,file="out/bm.el2.tab",quote=F,sep="\t",row.names=F,col.names=T)


	library(igraph)
	ig<-graph.adjacency(b.net[,])
	cig<-optimal.community(ig)
	pdf("/Users/bambrose/Dropbox/2014-2015/Sprints/1/BOURDIEU, 1985, THEOR SOC/out/ml_igraph_opt.pdf")
	plot(cig,ig)
	dev.off()



	## then extract basic data on each citation
	rawbel[,`:=`(
		y=as.integer(sub("(^|(.*, ))((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9]))($|, .*)","\\3",cr)) # year number, or na
		,v=as.integer(sub(".*, V([0-9]+).*","\\1",cr)) # volume number, or na
		,p=as.integer(sub(".*, V([0-9]+).*","\\1",cr)) # page number, or na
		,s=as.numeric(strptime(
	paste(sub("(^|(.*, ))((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9]))($|, .*)","\\3",cr),sub(".*(((0[1-9])|(1[0-2]))((0[1-9])|([1-2][0-9])|(3[0-1]))).*","\\1",cr),sep="")
	,"%Y%m%d"))/60/60/24 # periodical date in days
		,l=nchar(cr) # length of string
		,ca=grepl("^\\*",cr)
		,b09=grepl("^((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9]))",cr)
	)]
}

return(rawbel)
}

if(F){
t2<-proc.time()
trg<-data.table(dl=sapply(hd,function(x) any(grepl(" ?((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9])),",labels(x))))) # has date (all true)
trg[,vl:=sapply(hd,function(x) any(grepl(", V[0-9]",labels(x))))] #has volume
trg[,pl:=sapply(hd,function(x) any(grepl(", P[0-9]+[A-Z]?$",labels(x))))] #has page (always last)
trg[,sl:=sapply(hd,function(x) any(grepl(".*(((0[1-9])|(1[0-2]))((0[1-9])|([1-2][0-9])|(3[0-1]))).*",labels(x))))] #daily serial
trg[,f3l:=sapply(hd,function(x) length(unique(substr(labels(x),1,3)))==1)]

trg[,"d":=list(lapply(hd,function(x) {x<-na.omit(as.integer(sub(".*((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9])),.+","\\1",labels(x))));attributes(x)<-NULL;x}))]
trg[,"v":=list(lapply(hd,function(x) {x<-na.omit(as.integer(sub(".+, V([0-9]+).*","\\1",labels(x))));attributes(x)<-NULL;x}))]
trg[,"p":=list(lapply(hd,function(x) {x<-na.omit(as.integer(sub(".+, P([0-9]+).*","\\1",labels(x))));attributes(x)<-NULL;x}))]
trg[,"s":=list(lapply(hd,function(x) {x<-na.omit(as.numeric(strptime(
	paste(sub(".*((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9])).*","\\1",labels(x)),sub(".*(((0[1-9])|(1[0-2]))((0[1-9])|([1-2][0-9])|(3[0-1]))).*","\\1",labels(x)),sep="")
	,"%Y%m%d")));attributes(x)<-NULL;x<-x/60/60/24;x
}))] ## daily serial in days
trg[,"h":=list(lapply(hd,get_branches_heights,sort=F))]
trg[,"nc":=list(lapply(hd,function(x) nchar(labels(x))))]

trg[,mxd:=sapply(d,max)]
trg[,mnd:=sapply(d,min)]
trg[,mxv:=sapply(v,max)]
trg[,mnv:=sapply(v,min)]
trg[,mxp:=sapply(p,max)]
trg[,mnp:=sapply(p,min)]
trg[,mxh:=sapply(h,max)]
trg[,mnh:=sapply(h,min)]
trg[,mxnc:=sapply(nc,max)]
trg[,mnnc:=sapply(nc,min)]
trg[,l:=sapply(hd,function(x) length(labels(x)))]

trg[is.infinite(mxv),mxv:=NA]
trg[is.infinite(mnv),mnv:=NA]
trg[is.infinite(mxp),mxp:=NA]
trg[is.infinite(mnp),mnp:=NA]


trg[,`:=`(
	dsd=sapply(d,sd)
	,vsd=sapply(v,sd)
	,psd=sapply(p,sd)
	,ssd=sapply(s,sd)
	,hsd=sapply(h,sd)
	,ncsd=sapply(nc,sd)
)]

trg[,baz:=sapply(hd,function(x) any(grepl("^[A-Z]",labels(x))))] #begin letter
trg[,b09:=sapply(hd,function(x) any(grepl("^((((17)|(18)|(19))[0-9]{2})|(((200)|(201))[0-9]))",labels(x))))] #begin date
trg[,ca:=sapply(hd,function(x) all(grepl("^\\*",labels(x))))] #corp author source
trg[,aca:=sapply(hd,function(x) any(grepl("^\\*",labels(x))))] #mixed corp author source
trg[,aca:=aca&!ca]

samp.c<-list()
samp.c[which((is.na(trg$vsd))&is.na(trg$psd)&!trg$sl)]<-"nvnp"
samp.c[which((!is.na(trg$vsd))&is.na(trg$psd)&!trg$sl)]<-"yvnp"
samp.c[which((is.na(trg$vsd))&!is.na(trg$psd)&!trg$sl)]<-"nvyp"
samp.c[which((!is.na(trg$vsd))&!is.na(trg$psd)&!trg$sl)]<-"yvyp"
samp.c[which(trg$sl)]<-"serial"
samp.c<-factor(unlist(samp.c))
trg[,samp.c:=samp.c]

trg[,dsd0:=dsd==0]
trg[,vsd0:=vsd==0]
trg[,psd0:=psd==0]
trg[,ssd0:=ssd==0]
trg[,l2:=l==2]
trg[,mxv3:=mxv<=3]
trg[,osamp:=sapply(hd,function(x) sort(unique(unlist(attributes(x)$osamp))))]
trg[,og:=sapply(hd,function(x) attributes(x)$og)]

t3<-proc.time()
round((t3-t2)/60,2)
round((t3-t0)/60,2)

trg[,osampl1:=unlist(sapply(osamp,min))]
trg[sapply(osamp,function(x) length(x)>1),osampl1:=list(NA)]
trg[,osampl1:=unlist(osampl1)]
setkey(trg,osampl1)
trg[,nsamp:=1:nrow(trg)]
}

if(F) {
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/sample.RData")
	attributes(trg)$samp<-samp
	save(trg,file="/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/triage.RData")
}

##########

#calculate quantiles
if(F) {
cen<-lapply(trg[,list(dsd,vsd,psd,ssd,hsd,ncsd,mxd,mnd,mxv,mnv,mxp,mnp,mxh,mnh,mxnc,mnnc,l)],function(x) quantile(x[is.finite(x)],p=seq(0.01,1,.01),na.rm=T))

names(cen)<-c("d","v","p","s","h","nc","mxd","mnd","mxv","mnv","mxp","mnp","mxh","mnh","mxnc","mnnc","l")
trg[,`:=`(
	dsp=sapply(dsd,function(x) max(which(cen$d<=x))) # lower is better
	,vsp=sapply(vsd,function(x) max(which(cen$v<=x))) # lower is better
	,psp=sapply(psd,function(x) max(which(cen$p<=x))) # lower is better
	,ssp=sapply(ssd,function(x) max(which(cen$s<=x))) # lower is better
	,hsp=sapply(hsd,function(x) max(which(cen$h<=x))) # lower might not be better, but absolute is more useful
	,ncsp=sapply(ncsd,function(x) max(which(cen$nc<=x))) # lower might not be better, but absolute

	,mxvsp=sapply(mxv,function(x) max(which(cen$mxv<=x))) # if higher than three, assume not a book, percentile not very useful
	,mnvsp=sapply(mnv,function(x) max(which(cen$mnv<=x))) # if higher than three, assume not a book, percentile not very useful
	,mxpsp=sapply(mxp,function(x) max(which(cen$mxp<=x))) # won't use, if one assume intro
	,mnpsp=sapply(mnp,function(x) max(which(cen$mnp<=x))) # won't use
	,mxhsp=sapply(mxh,function(x) max(which(cen$mxh<=x))) # lower is better
	,mnhsp=sapply(mnh,function(x) max(which(cen$mnh<=x))) # lower is better
	,mxncsp=sapply(mxnc,function(x) 100-max(which(cen$mxnc<=x))) # higher is better, so subtracted percentile from 100, so lower is better. Min is more useful
	,mnncsp=sapply(mnnc,function(x) 100-max(which(cen$mnnc<=x))) # higher is better, so subtracted percentile from 100, so lower is better. Use instead of max

	,lsp=sapply(l,function(x) 100-max(which(cen$l<=x))) # higher is better, so subtracted percentile from 100, so lower is better.
)]

trg[dsd==0,dsp:=0]
trg[vsd==0,vsp:=0]
trg[psd==0,psp:=0]
trg[ssd==0,ssp:=0]

trg[is.na(dsd),dsd:=Inf]
trg[is.na(vsd),vsd:=Inf]
trg[is.na(psd),psd:=Inf]
trg[is.na(ssd),ssd:=Inf]
trg[is.na(hsd),hsd:=Inf]
trg[is.na(mxv),mxv:=Inf]
trg[is.na(mnv),mnv:=Inf]
trg[is.na(mxp),mxp:=Inf]
trg[is.na(mnp),mnp:=Inf]

trg[dsp==-Inf,dsp:=Inf]
trg[vsp==-Inf,vsp:=Inf]
trg[psp==-Inf,psp:=Inf]
trg[ssp==-Inf,ssp:=Inf]
trg[hsp==-Inf,hsp:=Inf]
trg[ncsp==-Inf,ncsp:=Inf]
trg[mxvsp==-Inf,mxvsp:=Inf]
trg[mnvsp==-Inf,mnvsp:=Inf]
trg[mxpsp==-Inf,mxpsp:=Inf]
trg[mnpsp==-Inf,mnpsp:=Inf]
trg[mxhsp==-Inf,mxhsp:=Inf]
trg[mnhsp==-Inf,mnhsp:=Inf]
trg[mxncsp==-Inf,mxncsp:=Inf]
trg[mnncsp==-Inf,mnncsp:=Inf]

trg[,ap:=mapply(
	function(dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp)
		mean(c(dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp)[
			is.finite(c(dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp))
		])
	,dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp
)]

trg[,mp:=mapply(
	function(dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp)
		max(c(dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp)[
			is.finite(c(dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp))
		])
	,dsp,vsp,psp,ssp,hsp,ncsp,mxvsp,mnvsp,mxpsp,mnpsp,mxhsp,mnhsp,mxncsp,mnncsp,lsp
)]

t1<-proc.time()
t1-t0
lapply(hd[head(order(trg$ap),10)],labels)
trg[head(order(trg$ap),10)]
lapply(hd[tail(order(trg$ap),10)],labels)
trg[tail(order(trg$ap),10)]
#attributes(trg)$samp<-sample(1:nrow(trg),1000)
for(i in c("dsd","vsd","psd")) trg[,paste(i,"z",sep=""):=trg[[i]]==0,with=F]
data.frame(cen)
lapply(trg[,list(dsd,vsd,psd,ssd,hsd)], hist)
#surefire?
lapply(hd[sample(which(trg$dsd==0&trg$vsd==0&trg$psd==0&trg$ncsd==0&trg$l>2),10)],function(x) cat(c(labels(x),"\n"),sep="\n"))
c<-75
t<-function(c) (trg$dsd<=cen$d[c]|is.infinite(trg$dsd))&(trg$vsd<=cen$v[c]|is.infinite(trg$vsd))&(trg$psd<=cen$p[c]|is.infinite(trg$psd))&(trg$hsd<=cen$h[c])&(trg$ncsd<=cen$nc[c])
lapply(hd[sample(which(t(c)),10)],function(x) cat(c(labels(x),"\n"),sep="\n"))
where2drawtheline<-data.frame(m=round(sapply(1:100,function(c) mean(t(c)))*100,1),n=sapply(1:100,function(c) sum(t(c))))
plot(where2drawtheline$n,type="l")
t75<-t(c)
data.frame(lapply(cen,round,3))

}

if(F){
#dl has date
#vl has volume
#pl has page
#sl daily serial
#d date list
#v volume list
#p page list
#s date in seconds
#h dedrogram heights list
#mxv max volume
#mnv min volume
#mxh max height
#mnh min height
#l length, number of leaves/cases/citations
#dsd date standard deviation
#vsd volume standard deviation
#psd page standard deviation
#ssd seconds standard deviation
#hsd height standard deviation
#dsp percentile of date standard deviation
#vsp percentile of volume standard deviation
#psp percentile of page standard deviation
#ssp percentile of seconds standard deviation
#hsp percentile of heights standard deviation
#mxvsp percentile of max volume
#mnvsp percentile of min volume
#mxhsp percentile of max height
#mnhsp percentile of min height
#lsp percentile of length/leaves/citations
#ap average of percentiles, omitting Inf
#baz log begins with a letter
#b09 log beings with date
#ca corp author source
#aca mixed corporate author source


##


#######

load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/FirstTraining_nvnp+/mastersets.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/triage.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hd.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/sample.RData")
library(data.table)

hd<-hd[attributes(trg)$samp]
train<-trg[attributes(trg)$samp]
train[,samp:=samp]
train[,good:=sapply(sets,length)==1&sapply(sets,function(x) length(unlist(x)))==sapply(hd,function(x) attributes(x)$members)]
train<-train[,list(good,vl,pl,sl,dsd0,vsd0,psd0,dsd,vsd,psd,l,l2,f3l,mxv3,ncsd,baz,mxh,mnnc,samp)]

#### really just do 4 different models to avoid missingness

(less<-table(data.frame(vna=is.na(train$vsd)&!train$sl,pna=is.na(train$psd)&!train$sl)))
round(prop.table(less)*100,2)
}
### expanded sample
if(F){
nvnp<-list(trgl=(is.na(trg$vsd))&is.na(trg$psd)&!trg$sl,trainl=(is.na(train$vsd))&is.na(train$psd)&!train$sl)
yvnp<-list(trgl=(!is.na(trg$vsd))&is.na(trg$psd)&!trg$sl,trainl=(!is.na(train$vsd))&is.na(train$psd)&!train$sl)
nvyp<-list(trgl=(is.na(trg$vsd))&!is.na(trg$psd)&!trg$sl,trainl=(is.na(train$vsd))&!is.na(train$psd)&!train$sl)
yvyp<-list(trgl=(!is.na(trg$vsd))&!is.na(trg$psd)&!trg$sl,trainl=(!is.na(train$vsd))&!is.na(train$psd)&!train$sl)

nvnp$samp$old<-intersect(samp,which(nvnp$trgl))
nvnp$samp$new<-integer(0)

yvnp$samp$old<-intersect(samp,which(yvnp$trgl))
yvnp$samp$new<-sample(setdiff(which(yvnp$trgl),samp),500-length(yvnp$samp$old))

nvyp$samp$old<-intersect(samp,which(nvyp$trgl))
nvyp$samp$new<-sample(setdiff(which(nvyp$trgl),samp),500-length(nvyp$samp$old))

yvyp$samp$old<-intersect(samp,which(yvyp$trgl))
yvyp$samp$new<-sample(setdiff(which(yvyp$trgl),samp),500-length(yvyp$samp$old))

cleans<-list(nvnp=nvnp,yvnp=yvnp,nvyp=nvyp,yvyp=yvyp)
#whoops<-cleans
if(F) save(whoops,file="/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/whoops.RData")

}



if(F){

######## oh dear
bad<-read.table(file = "1900-2010FuzzySets/yvnp/whoops_bad.tab", sep = "\t", header = F, stringsAsFactors = FALSE)
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/yvnp/mastersets.RData")

source("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/nvyp/nvyp_bads.R")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/nvyp/mastersets.RData")
bad<-nvypb

search<-sapply(sets,function(x) unlist(x)[1])
if(all(bad[[1]]==which(sapply(sets,function(x) length(unlist(x)))<2))){
	for(i in 1:nrow(bad)) if(is.null(search[[bad[i,1]]])) search[[bad[i,1]]]<-as.character(bad[i,2])
}
search<-unlist(search)

#load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hd.RData")
hdi<-sapply(hd,function(x) labels(x)[1])
library(data.table)
hdi<-data.table(i=1:length(hdi),cr=unlist(hdi))
#whps<-hdi[search]

library(stringdist)
sdm<-stringdistmatrix(search,hdi$cr,useNames=T,ncores=4,method="jw",p=.1)
sdml<-apply(sdm,1,function(x) which(x==min(x)))
setkey(trg,osampl1)
trg[list(sapply(sdml,function(x) x[1])),.N,by=og]
nsdml<-list()
for(i in which(!sdml%in%trg$osampl1)){
	nsdml[[i]]<-min(unlist(trg$osamp[sapply(trg$osamp,function(x) sdml[[i]][1]%in%x)]))
}
if(all(which(!sdml%in%trg$osampl1)==which(sapply(nsdml,length)==1))) sdml[!sdml%in%trg$osampl1]<-unlist(nsdml)
trg[list(sapply(sdml,function(x) x[1])),.N,by=og]

### update last whoops
#load(grep("whoops.RData",list.files("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data",recursive=T,full.names=T),value=T))
whoops$nvnp$samp$new<-integer(0)
lapply(whoops,function(x) sapply(x$samp,length))
setkey(trg,nsamp)
#whoops$yvnp$samp$new<-trg[list(sdml)]$osampl1
whoops$nvyp$samp$new<-trg[list(sdml)]$osampl1
setkey(trg,osampl1)
trg[,hand:=F]
trg[list(unlist(lapply(whoops,function(x) x$samp))),hand:=T]

##### attach samples to trg
setl<-grep("mastersets",list.files("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets",recursive=T,full.names=T),value=T)
for(i in setl){
	load(i)
	whoops[[sub(".+([ny]v[ny]p).+","\\1",i)]]$sets<-sets
}





save(trg,file="triage.RData")




read.table()

	sets<-list()
	up<-list()
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/cleans.RData")
#	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/mastersets.RData")
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/progress.RData")
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hd.RData")
	library(dendextend)
	library(data.table)
	hd<-hd[yvyp$samp$new]

	nlhd<-sort(sapply(hd,nleaves),decreasing=T)
	tt<-cumsum(nlhd)
	tt<-sum(nlhd)-tt

	source("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/strdist.dend.picker.R")
	for(i in (length(sets)+1):length(hd)){
		up$beg[[i]]<-as.character(Sys.time())
		t0<-proc.time()
		#Sys.sleep(rpois(1, lambda = .5)+1)
		sets[[i]]<-strdist.dend.picker(hd[[i]],out=paste("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets",.Platform$file.sep,"hd",i,"-",sep=""))
		t1<-proc.time()
		up$et[[i]]<-t1-t0
		up$rps[[i]]<-nlhd[i]/up$et[[i]]["elapsed"] #rate (per second)

		mph<-mean(unlist(up$rps))*60*60
		rph<-up$rps[[i]]*60*60
		cat(
		ifelse(rph>=mph,"\n:) ","\n:( "),round(i/length(hd)*100,1)," Rate(Avg): ",round(rph,1),"(",round(mph,1),") /hr\tFinished in ",round(tt[i]/rph,1),"(",round(tt[i]/mph,1),")"," hrs or ",round(tt[i]/rph/4,1),"(",round(tt[i]/mph/4,1),")"," 4 hr days"
		,sep="")
		flush.console()
		save(sets,file="/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/mastersets.RData")
		save(up,file="/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/progress.RData")
	}

### analysis

load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/sample.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hd.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/triage.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/cleans.RData")
load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets/FirstTraining_nvnp+/mastersets.RData")

for(i in c("nvnp","yvnp","nvyp","yvyp")){
traino<-train[cleans[[i]]$trainl,list(good,dsd0,vsd0,psd0,dsd,vsd,psd,l,l2,f3l,mxv3,ncsd,baz,mxh,mnnc)]
load(paste(grep(i,list.dirs("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1900-2010FuzzySets",recursive=F),value=T),"mastersets.RData",sep=.Platform$file.sep))
hdn<-hd[cleans[[i]]$samp$new]
trainn<-trg[cleans[[i]]$samp$new]
trainn[,good:=sapply(sets,length)==1&sapply(sets,function(x) length(unlist(x)))==sapply(hdn,function(x) attributes(x)$members)]
trainn<-trainn[,list(good,dsd0,vsd0,psd0,dsd,vsd,psd,l,l2,f3l,mxv3,ncsd,baz,mxh,mnnc)]
trainn[,samp:=cleans[[i]]$samp$new]
cleans[[i]]$dat<-rbindlist(list(traino,trainn),use.names=T,fill=T)
narms<-names(cleans[[i]]$dat)[sapply(cleans[[i]]$dat,function(x) any(is.na(x)))]
for(j in narms) cleans[[i]]$dat[,j:=NULL,with=F]
}

library(MASS)

for(i in c("nvnp","yvnp","nvyp","yvyp")){
	cleans[[i]]$fit<-glm(good~.
		#as.formula(paste(sub("\\+","~",paste(grep("0",names(cleans[[i]]$dat),value=T,invert=T),collapse=" + ")),paste("(",paste(grep("0",names(cleans[[i]]$dat),value=T),collapse=" + "),")^2",sep=""),"0",sep=" + "))
		,family=binomial(link="cloglog"),data=cleans[[i]]$dat)
	cleans[[i]]$stp<-stepAIC(cleans[[i]]$fit,direction="both")
	#cleans[[i]]$stp<-cleans[[i]]$fit
	cleans[[i]]$hat<-predict(cleans[[i]]$stp,cleans[[i]]$dat,type="response")
	cleans[[i]]$hatc<-cut(cleans[[i]]$hat,breaks=seq(0,1,.1))
	cleans[[i]]$tab$c<-table(data.frame(cleans[[i]]$hatc,cleans[[i]]$dat$good))
	cleans[[i]]$tab$r<-round(prop.table(cleans[[i]]$tab$c,margin=1)*100,2)
	cleans[[i]]$tab$t<-round(prop.table(cleans[[i]]$tab$c)*100,2)
}

for(i in c("nvnp","yvnp","nvyp","yvyp")) {
	cat("\n\n#################################",i,round(mean(cleans[[i]]$trgl)*100,2),"\b%",sum(cleans[[i]]$trgl),"#################################\n\n")
	print(cbind(
	r=cleans[[i]]$tab$r
	,t=cleans[[i]]$tab$t
	,c=cleans[[i]]$tab$c))
	slt<-split(lapply(hd[c(cleans[[i]]$samp$old,as.integer(cleans[[i]]$samp$new))],labels),f=list(cleans[[i]]$dat$good,cleans[[i]]$hatc))
	for(j in names(slt)) {cat("\n");cat(j,unlist(sample(ifelse(!length(slt[[j]]),"none",slt[[j]]),1)),sep="\n")}
	print(summary(cleans[[i]]$fit))
	print(summary(cleans[[i]]$stp))
	print(cleans[[i]]$stp$anova)
}

save(cleans,file="clean_fits.RData")




for(i in c("nvnp","yvnp","nvyp","yvyp")){
	cleans[[i]]$tab$c<-table(data.frame(cut(cleans[[i]]$hat,10),cleans[[i]]$dat$good))
	cleans[[i]]$tab$r<-round(prop.table(cleans[[i]]$tab$c,margin=1)*100,2)
	cleans[[i]]$tab$t<-round(prop.table(cleans[[i]]$tab$c)*100,2)
}


library(mi)
mitrainog<-data.table(data.frame(lapply(trainog, function(x) as.numeric(replace(x, is.infinite(x),NA)))))
mitrainog<-mi(mitrainog,n.imp=10,n.iter=100)

library(SuperLearner)
slfit<-SuperLearner(Y=as.numeric(trainog$good),X=,SL.library=c("SL.glm"),family=binomial(link="cloglog"),verbose=T)
}

if(F){
terms<-list()
for(j in 1:10){
#impute missing data, badly
for(i in which(sapply(train,function(x) any(is.infinite(x))))) train[is.infinite(train[[i]]),i:=sample(train[[i]][is.finite(train[[i]])],sum(is.infinite(train[[i]])),replace=T),with=F]
#http://plantecology.syr.edu/fridley/bio793/glm.html
terms[[j]]<-cllfit.stp$formula
cat("\r",round(j/10,3),sep="")
}
terml<-lapply(terms,function(x) attributes(terms(x))$term.labels)

(form<-as.formula(paste("good",paste(unique(unlist(terml)),collapse=" + "),sep=" ~ ")))

train<-copy(trainog)
for(i in which(sapply(train,function(x) any(is.infinite(x))))) train[is.infinite(train[[i]]),i:=sample(train[[i]][is.finite(train[[i]])],sum(is.infinite(train[[i]])),replace=T),with=F]


cllfit.stp<-glm(form,family=binomial(link="cloglog"),data=train)
coef<-list(cllfit.stp$coefficients)
train<-copy(trainog)
t2<-proc.time()
for(j in 2:1000) {
	for(i in which(sapply(train,function(x) any(is.infinite(x))))) train[is.infinite(train[[i]]),i:=sample(train[[i]][is.finite(train[[i]])],sum(is.infinite(train[[i]])),replace=T),with=F]
	coef[[j]]<-glm(form,family=binomial(link="cloglog"),data=train)$coefficients
	train<-copy(trainog)
	cat("\r",round(j/1000,3),sep="")
}
t3<-proc.time()
(t3-t2)/60
coef<-data.frame(do.call(rbind,coef))
pkdens.coef<-apply(coef,2,function(x) {y<-density(x);y<-y$x[which.max(y$y)];y})
den<-apply(coef,2,density)
}
}
