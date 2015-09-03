wok2db.f<-function(
	dir=stop("Choose input directory containing WOK batches.")
	,out=stop("Specify output directory for your project.")
	,art_rev_only=T
	,sample.batches=F
	,sample.size=1000
	,save=T
	,verbose=T
	,check.for.saved.output=F
)
{
	if(check.for.saved.output) if(any(grepl('wok2db.RData',dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning('Loading and returning saved wok2db.RData.',call. = F)
		load(dir(pattern='wok2db.RData',recursive=T,full.names=T,ignore.case=F)[1])
		return(wok2db)
	}

	#NEWER WOK Database Import
	require(data.table)

	files<-list.files(dir,full.names=T,recursive=T)

	c<-0
	n<-length(files)

	cat("\n",n," batches detected.",sep="")

	if(sample.batches) {
		files<-sort(sample(x=files,size=sample.size))
		cat("\n",sample.size," or ",round(sample.size/n*100,3)," % of batches drawn at random.\n\n",sep="")
	}

	n<-length(files)

	wok2db<-list()
	for(i in files){
		c<-c+1
		if(verbose) {flush.console();cat("\r",round(c/n,3),i,sep=" ")}
		b<-readLines(i,warn=F)
		field<-sub("^(.{2}).+","\\1",b)
		cut<-field%in%c("FN","VR","","ER","EF")
		b<-b[!cut]
		field<-field[!cut]
		t<-field=="  "
		x<-which(t)
		y<-which(!t)
		if(any(t)) {
			for(j in 1:length(x)) x[j]<-y[sum(y<x[j])]
			field[t]<-field[x]
		}
		b<-sub("^.. ?(.*)","\\1",b)
		ind<-1:length(b)
		t<-field!="UT"
		x<-which(t)
		y<-which(!t)
		ind[t]<-ind[y[sapply(lapply(lapply(as.list(x),"<",y),which),min)]]
		d<-data.table(b.ind=b[ind],field,b)
		if(art_rev_only){
			dt<-which(d$field=="DT")
			dtt<-grepl("(Article)|(Review)",d$b[dt])
			if(!all(dtt)) {
			dt<-d$b.ind[dt[dtt]]
			setkey(d,b.ind)
			d<-d[dt]
			}
		}
		wok2db[[i]]<-copy(d)
	}
	wok2db<-rbindlist(wok2db)
	setkey(wok2db,b.ind,field)
	dup<-duplicated(wok2db)
	if(any(dup)){
		ud<-unique(wok2db$b.ind[dup])
		cat('\n\n',length(ud),' duplicate records detected and removed, e.g.:\n',sep='')
		if(length(ud)>5) {cat(sample(ud,5),sep='\n')} else {cat(ud,sep='\n')}
		wok2db<-unique(wok2db)
	}
	wok2db<-droplevels(wok2db)
	setnames(wok2db,c("b.ind","b"),c("id","val"))
	if(save) save(wok2db,file=paste(out,.Platform$file.sep,"wok2db.RData",sep=""))
	wok2db
}

jstor2db.f<-function(
	dir=stop("Choose input directory containing DFR.JSTOR.ORG batches of zip archives.")
	,out=stop("Specify output directory for your project.")
	,sample.batches=T
	,sample.size=1
	,in.parallel=F # import batches in parallel
	,drop.nchar1=T # drop ngrams of 1 character
	,drop.freq1=T # drop ngrams that appear only once
	,check.for.saved.output=F # will scan output directory for a 'jstor2db.RData' file, load it, and return it instead of running a new import
)
{

	if(check.for.saved.output) if(any(grepl('jstor2db.RData',dir(recursive=T,full.names=T,ignore.case=T)))) {
		warning('Loading and returning saved jstor2db.RData.',call.=F)
		load(dir(pattern='jstor2db.RData',full.names=T,recursive=T,ignore.case=F)[1])
		return(wok2db)
	}

	zips<-grep('\\.zip$',list.files(dir,full.names=T,recursive=T,include.dirs=T),value=T)
	n<-length(zips)
	if(sample.batches) {
		zips<-sort(sample(x=zips,size=sample.size))
		cat("\n",sample.size," or ",round(sample.size/n*100,3)," % of batches drawn at random.\n\n",sep="")
	}

	import.dfr.jstor.f<-function(zip){
		require(tm)
		require(SnowballC)
		require(data.table)

		temp <- tempdir()
		unzip(zip,exdir=temp)
		f<-list.files(temp,recursive=T,include.dirs=T,full.names=T)
		dat<-try(data.table(read.table(grep('citations.tsv',f,value=T),header=F,skip=1,sep='\t',quote='',comment.char = "")))
		if(class(dat)[1]=="try-error") return(dat)
		dat[,ncol(dat):=NULL,with=F]
		setnames(dat,as.character(read.table(grep('citations.tsv',f,value=T),header=F,nrows=1,sep='\t',quote='',as.is=T)))
		dat[,id:=sub('/','_',id)]
		setkey(dat,id)

		bow.f<-function(x){
			x<-data.table(read.csv(grep(x,f,value=T),as.is=T))
			x<-x[!removeWords(x$WORDCOUNTS,stopwords('english'))=='']
			x[,WORDCOUNTS:=tolower(WORDCOUNTS)]
			x[,WORDCOUNTS:=removePunctuation(WORDCOUNTS)]
			x[,WORDCOUNTS:=removeNumbers(WORDCOUNTS)]
			x[,WORDCOUNTS:=stemDocument(WORDCOUNTS,language='english')]
			x<-x[,list(WEIGHT=sum(WEIGHT)),by=WORDCOUNTS]
			x<-x[,c(rep(WORDCOUNTS,WEIGHT))]
			x<-table(x)
			x
		}

		dat[,bow:=lapply(id,function(x) try(bow.f(x)))]
		unlink(temp)
		dat
	}

	if(in.parallel){
		require(doParallel)
		cl <- makeCluster(detectCores() )
		registerDoParallel(cl, cores = detectCores() )
		jstor2db <- foreach(i = zips,.packages = c('data.table','tm','SnowballC'),.inorder=F) %dopar% try(import.dfr.jstor.f(zip=i))
		stopCluster(cl)
	} else {jstor2db<-list();for(i in zips) jstor2db[[i]]<-try(import.dfr.jstor.f(zip=i))}

jstor2db<-rbindlist(jstor2db[sapply(jstor2db,is.data.table)])

	# condense vocab
if(drop.nchar1) jstor2db[,bow:=list(lapply(bow,function(x) as.table(x[nchar(names(x))>1])))]
if(drop.freq1) jstor2db[,bow:=list(lapply(bow,function(x) as.table(x[x>1])))]
	vocab<-sort(unique(unlist(lapply(jstor2db$bow,function(x) names(x)))))
jstor2db[,bow:=list(lapply(bow,function(x) {x<-matrix(c(which(vocab%in%names(x)),as.integer(x)),nrow=2);rownames(x)<-c('vocab.index','freq');x}))]

	dfr.jour<-jstor2db[,.N,by=journaltitle]
	print(dfr.jour)
	attributes(jstor2db)$vocab<-vocab
	save(jstor2db,file=paste(out,'jstor2db.RData',sep=.Platform$file.sep))
	jstor2db
}


#Port fuzzy set functions to here
if(F){if(F){
cr2zcr.f<-function( #utility for resolving identity uncertainty for WOK CR field
	cr=stop("cr = supply a character vector of unique and uncertain CR ids") #unprocessed bimodal edgelist in UT CR order
	,just.normalize=T # TRUE will remove DOI and capitalize. FALSE will perform fuzzy set replacement.
)
{
require(stringdist)
require(randomForest)
load("/Users/bambrose/Dropbox/2014-2015/Sprints/1/BOURDIEU, 1985, THEOR SOC/out/wok2db.RData")
if(!is.data.table(wok2db)) wok2db<-data.table(wok2db)
setkey(wok2db,field,id)

### impose formatting and nomenclature ###
rawbel<-wok2db["CR",c(1,3),with=F] #delete
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

db2bel.f<-function(
	wok2db
	,out=stop("Specify output directory for your project.")
	,saved_recode=NULL
	,save_og4recode=F
	,trim_doi=T
	,capitalize=T
	,trim_anon=T
	,cut_samp_def=0 # if sample is based on a common citation set that should be removed for analysis
	,trim_pendants=T
)
{

### check function requirements ###
require(data.table)
out

### draw only citation edge information from wok2db ###
wok2db<-data.table(wok2db)
setkey(wok2db,field)
db2bel<-wok2db["CR",list(id,val)]
rm("wok2db")

### impose formatting and nomenclature ###
setnames(db2bel,c("ut","cr"))
if(trim_doi) db2bel[,cr:=sub(", DOI .+","",cr)] #remove DOI
if(capitalize) db2bel[,cr:=gsub("(\\w)","\\U\\1",cr,perl=T)] #capitalize
if(trim_anon&capitalize) db2bel[,cr:=sub("^\\[ANONYMOUS\\], ","",cr)] #remove ANONYMOUS author

### recoding ###
setkey(db2bel,cr)
if(save_og4recode){ # save original codes to disk
	cat('\nSaving normalized original CR codes to pass to fuzzy set routine.')
	original.cr<-unique(dob2bel$cr)
	save(original.cr,file=paste(out,.Platform$file.sep,'original-cr.RData',sep=""))
}
if(!is.null(saved_recode)){ # recode using fuzzy sets
	lfs<-length(saved_recode)
	cat('\nRecoding CR from',lfs,'sets.\n')
	db2bel[,zcr:=cr]
	for(i in 1:lfs) {
		cat('\r',round(i/lfs*100,4),'%\t\t',sep='')
		ix<-saved_recode[[i]]
		recode<-db2bel[ix,list(cr)][,.N,by=cr]
		recode<-recode$cr[recode$N==max(recode$N)]
		recode<-tolower(recode[which.min(nchar(recode))])
		db2bel[ix,zcr:=recode]
	}
	cat('\n')
}
### pendants ###
if(trim_pendants){
db2bel[,pend:=!duplicated(cr)|duplicated(cr,fromLast=T)]
if(!is.null(saved_recode)) db2bel[,zpend:=!duplicated(zcr)|duplicated(zcr,fromLast=T)]
}

### cut sample ###
crd<-db2bel[,list(cr)][,.N,by=cr]
setkey(crd,N)
if(!is.null(saved_recode)) {
	setkey(db2bel,ut,zcr)
	db2bel[,zdup:=duplicated(db2bel)]
	zcrd<-db2bel[!db2bel$zdup,list(zcr)][,.N,by=zcr]
	setkey(zcrd,N)
	setkey(db2bel,cr)
}

if(cut_samp_def>0){
	#cut highest degree citations, argument is the length of the list in descending order of frequency
	cat("\nEnter -indices separated by spaces- to reject high degree citations such as those defining the sample, or -enter- to reject none.\n")
	cut<-crd[nrow(crd):(nrow(crd)-cut_samp_def),list(cr=cr,index=1:cut_samp_def,N=N)]
	print(cut)
	u<-readLines(n=1)
	if(u!=""){
		u<-unlist(strsplit(u," "))
		u<-as.integer(u)
		if(any(is.na(u))) {stop("\nTry again.")} else {db2bel<-db2bel[!list(cut[u]$cr)]}
	}
}

### results ###
print(crd)
cat('\n')
if(!is.null(saved_recode)) print(zcrd)
cat('\n')

setkey(crd,N)
if(!is.null(saved_recode)) setkey(zcrd,N)

res<-data.frame(case=c('Total acts of reference'
	,'Unique citations'
	,'Citations referenced once'
	,'Total referenced twice or more'
	,'  Unique referenced twice or more'
	)
,cr=c(crd[,sum(N)]
	,nrow(crd)
	,nrow(crd[list(1)])
	,crd[!list(1),sum(N)]
	,nrow(crd[!list(1)])
 )
)
if(!is.null(saved_recode)){
res$zcr<-c(zcrd[,sum(N)]
	,nrow(zcrd)
	,nrow(zcrd[list(1)])
	,zcrd[!list(1),sum(N)]
	,nrow(zcrd[!list(1)])
)
res$change<-res$zcr-res$cr
res$perc.change<-round(res$zcr/res$cr*100,5)
res$alt.change<-res$perc.change-100
}
print(res)
mres<-data.frame(case=
	c('Pendants as % of total references'
	 ,'Pendants as % of unique citations')
	,cr=round(c(res$cr[3]/res$cr[1],res$cr[3]/res$cr[2])*100,4)

)
if(!is.null(saved_recode)) mres$zcr<-round(c(res$zcr[3]/res$zcr[1],res$zcr[3]/res$zcr[2])*100,4)
print(mres)
# db2bel[,`:=`(ut=factor(ut),cr=factor(cr))] # better to factor after selection is made
# if(!is.null(saved_recode)) db2bel[,zcr:=factor(zcr)]
save(db2bel,file=paste(out,.Platform$file.sep,"db2bel.RData",sep=""))
db2bel
}


bel2mel.f<-function(
	db2bel=NULL
	,subset=NULL # a vector of UT
	,type=c("utel","crel")
	,out=stop("Specify output directory for your project.")
	,trim_pendants=T
)
{
	require(data.table)
	if(ncol(db2bel)>2) stop('db2bel must be a bimodal edgelist as a two column data.table. Selection on pendants, etc. should be made prior to passing to bel2mel.f.')

	setnames(db2bel,c('ut','cr'))
	db2bel[,ut:=as.character(ut)]
	db2bel[,cr:=as.character(cr)]

	bel2mel<-list()

	if('utel'%in%type){
		setkey(db2bel,ut,cr)
		utd<-db2bel[,.N,by=ut]
		setkey(utd,N,ut)
		utiso<-utd[list(1),ut]
		bel2mel$utel<-db2bel[!list(utiso),list(cr=combn(cr,m=2,FUN=sort,simplify=F)),by=ut]
		bel2mel$utel[,`:=`(cr1=sapply(cr,function(x) x[1]),cr2=sapply(cr,function(x) x[2]))]
		bel2mel$utel[,cr:=NULL]
		bel2mel$utel<-bel2mel$utel[,list(ew=.N,ut=list(ut)),keyby=c('cr1','cr2')]
	}
	if('crel'%in%type){
		setkey(db2bel,cr,ut)
		crd<-db2bel[,.N,by=cr]
		setkey(crd,N,cr)
		criso<-crd[list(1),cr]
		bel2mel$crel<-db2bel[!list(criso),list(ut=combn(ut,m=2,FUN=sort,simplify=F)),by=cr]
		bel2mel$crel[,`:=`(ut1=sapply(ut,function(x) x[1]),ut2=sapply(ut,function(x) x[2]))]
		bel2mel$crel[,ut:=NULL]
		bel2mel$crel<-bel2mel$crel[,list(ew=.N,cr=list(cr)),keyby=c('ut1','ut2')]
	}
	save(bel2mel,file='bel2mel.RData')
	bel2mel
}

mel2cfinder.f<-function(
)
{
}

cfinder2all.f<-function(
	cf.in=stop('cf.out.dir = Path to CFinder output.')
	,which=c('communities_links'
					 ,'communities'
					 ,'communities_cliques'
					 ,'degree_distribution'
					 ,'graph_of_communities'
					 ,'membership_distribution'
					 ,'overlap_distribution'
					 ,'size_distribution')
)
{
	require(data.table)
	fileps<-list.files(cf.in,recursive=T,full.names=T)
	ret<-list()
	for(i in which){
		if(c('communities_links')%in%i){
			cat('compiling',i,'\n')

			makerocketgonow<-function(raw){
				raw<-readLines(raw)
				raw<-strsplit(raw[grep('^[0-9][^:]+$',raw)],' ')
				raw<-data.table(do.call(rbind,lapply(raw,function(x) sort(as.integer(x)))))
				setnames(raw,c('src','tgt'))
				raw<-unique(raw)
				raw
			}

			require(doParallel)
			cl <- makeCluster(detectCores() )
			registerDoParallel(cl, cores = detectCores() )

			ret[[i]] <- foreach(j = grep(paste('.+',.Platform$file.sep,i,'$',sep=''),fileps,value=T),.packages = c("data.table"),.inorder=F) %dopar% {
				makerocketgonow(raw=j)
			}

			stopCluster(cl)

			ret[[i]]<-rbindlist(ret[[i]])[,list(ew=.N),by=c('src','tgt')]
			setkey(ret[[i]],ew)
		}
		if(c('communities')%in%i){
			cat('compiling',i,'\n')

			makerocketgonow<-function(raw){
				stub<-sub('^.+k=([0-9]+).+$','k\\1',raw)
				raw<-readLines(raw)
				raw<-strsplit(grep('^[0-9]',raw,value=T),':? ')
				nams<-sapply(raw,function(x) x[1])
				raw<-lapply(raw,function(x) as.integer(x[-1]))
				names(raw)<-paste(stub,nams,sep='-')
				raw
			}

			require(doParallel)
			cl <- makeCluster(detectCores() )
			registerDoParallel(cl, cores = detectCores() )

			ret[[i]] <- foreach(j = grep(paste('.+',.Platform$file.sep,i,'$',sep=''),fileps,value=T),.packages = c("data.table"),.inorder=F) %dopar% {
				makerocketgonow(raw=j)
			}

			stopCluster(cl)

			names(ret[[i]])<-sapply(ret[[i]],function(x) sub('^([^-]+).+$','\\1',names(x[1])))
		}
	}
ret
}

perm.pois.f<-function(
	mel2net
	,nsim
)
{
	require(network)
	cat("\nPermuting random poisson edge distribution\nSimulating...")
	perm<-list()
	s<-network.size(mel2net)
	maxcombo<-s*(s-1)/2
	combos<-1:maxcombo
	choices<-sum(mel2net%e%"ew")
	t1<-proc.time()
	for(i in 1:nsim){
		cat("\r",i,sep="")
		rdist<-sample(combos,size=choices,replace=T)
		rdist<-table(rdist)
		z<-maxcombo-length(rdist)
		rdist<-table(rdist)
		n<-names(rdist)
		rdist<-c(z,rdist)
		names(rdist)<-c("0",n)
		rdist<-as.table(rdist)
		perm[[i]]<-rdist
	}
	cat("\nSeconds to simulate:")
	print(proc.time()-t1)

	edist<-table(mel2net%e%"ew")
	n<-length(unique(unlist(bel2mel[[i]][,1:2])))
	z<-(n*(n-1)/2)-sum(edist)
	n<-names(edist)
	edist<-c(z,edist)
	names(edist)<-c("0",n)
	edist<-as.table(edist)

	cn<-sort(unique(unlist(lapply(perm,names))))
	permdb<-data.frame(matrix(0,nrow=length(perm),ncol=length(cn)))
	permdb<-data.frame(permdb,matrix(0,nrow=length(perm),ncol=length(edist)-length(cn)))
	names(permdb)<-names(edist)
	for(i in 1:length(perm)) permdb[i,names(perm[[i]])]<-perm[[i]]

	dif<-matrix(edist,nrow=dim(permdb)[1],ncol=dim(permdb)[2],byrow=T)-permdb
	md<-apply(dif,2,mean)
	sdd<-apply(dif,2,sd)
	cid<-apply(dif,2,quantile,prob=c(.05,.95))
	tad<-apply(dif,2,table)
	num<-lapply(tad,names)
	num<-lapply(num,as.numeric)

	dens<-1-(permdb$`0`/maxcombo)
	mean((1-edist["0"]/maxcombo)-dens)*100
	sd((1-edist["0"]/maxcombo)-dens)*100
	mean((1-edist["0"]/maxcombo)-dens)/sd((1-edist["0"]/maxcombo)-dens)

	apply(t(apply(permdb,1,"*",as.numeric(colnames(permdb)))),sum)/maxcombo

	round(cbind(md,sdd,t=md/sdd,p=0,t(cid)),digits=3)

	round(cbind(md=md/maxcombo,sdd=sdd/maxcombo,t=md/sdd,p=0),digits=4)[1:7,]
	sum(edist[-(1:7)])/maxcombo

	plot(as.table((md/sdd)[1:7]),type="l",ylim=range((md/sdd)[1:7]),lwd=3)
	abline(h=0,lty="dotted",lwd=3)

	round(cbind(md,sdd,t=md/sdd,p=0),digits=1)[1:7,]
	sum(edist[-(1:7)])

	mp<-apply(dif/maxcombo,2,mean)
	sdp<-apply(dif/maxcombo,2,sd)
	cip<-apply(dif/maxcombo,2,quantile,prob=c(.05,.95))

	round(cbind(mp,sdp,p=0,t(cip)),digits=5)[1:7,]
	sum(edist[-(1:7)])

	mvd<-mean((sum(((as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2)*edist)/maxcombo)-(apply(apply(permdb,1,"*",(as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2),2,sum)/maxcombo))
	sddv<-sd((sum(((as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2)*edist)/maxcombo)-(apply(apply(permdb,1,"*",(as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2),2,sum)/maxcombo))
	cidv<-quantile((sum(((as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2)*edist)/maxcombo)-(apply(apply(permdb,1,"*",(as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2),2,sum)/maxcombo),prob=c(.05,.95))
	round(c(mvd=mvd,sddv=sddv,tvd=mvd/sddv,p=0,cidv=cidv),digits=4)

	m<-apply(permdb,2,mean)
	sd<-apply(permdb,2,sd)
	ci<-apply(permdb,2,quantile,prob=c(.05,.95))
	cbind(m,sd,ci)

	zscores<-(edist[1:7]-apply(permdb,2,mean))/apply(permdb,2,sd)

	apply(matrix(edist,nrow=dim(permdb)[1],ncol=dim(permdb)[2],byrow=T)-permdb,2,mean)
	apply(matrix(edist,nrow=dim(permdb)[1],ncol=dim(permdb)[2],byrow=T)-permdb,2,sd)
}

mel2net.f<-function(
	bel2mel
	,count=T
	,rcount=T
	,out=NULL
)
{
	cat("mel2net.f aka Plagiat!","Written by Brooks Ambrose\n",sep="\n")
	if(!count&rcount) warning("\nrcount iff count=T",call.=F)
	library(network)
	mel2net<-list()
	for(i in names(bel2mel)){
		if(is.na(bel2mel[[i]])) {mel2net[[i]]<-NA;next}
		bel2mel[[i]]<-bel2mel[[i]][order(bel2mel[[i]][,1],bel2mel[[i]][,2]),]
		mel2net[[i]]<-network(bel2mel[[i]][,1:2],matrix.type="edgelist",directed=F)
		mel2net[[i]]%e%"ew"<-bel2mel[[i]]$ew
		if(count) if("count"%in%names(attributes(bel2mel[[i]]))) {
			attributes(mel2net[[i]])$count<-attributes(bel2mel[[i]])$count
			if(rcount){

			}
		}

	}
	names(mel2net)<-sub("el$","",names(mel2net))
	if(!is.null(out)) save(mel2net,file=paste(out,"mel2net.RData",sep=.Platform$file.sep))
	mel2net
}

as.edgelist<-function(
	net
)
{
	require(network)
	el<-matrix(unlist(do.call(rbind,net$mel)[,2:1]),ncol=2)
	el<-cbind(
		s=network.vertex.names(net)[el[,1]]
		,r=network.vertex.names(net)[el[,2]]
	)
	el
}

plot.mode.projection<-function(
	db2bel
	,m1.stub="^m1"
	,out="/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/" #may add a stub at the end of this line
	,vcx=1
	,vlw=1
	,vlt="solid"
	,elw=1
	,elt="solid"
	,loopcx=1
	,m1vsides=3
	,m2vsides=4
	,m1col="black"
	,m2col="white"
	,ecol="gray"
	,trim.pendants=c("don't","m1","m2")
	,bmain="Bimodal"
	,m1main="Mode 1"
	,m2main="Mode 2"
	,mar=1
	,omi=.25
	,cex=1
	,pnt.v=list() #named list of character vectors where name is color and contents are vertices to paint
	,pnt.e=list() #named list of directed character edgelists where name is color and contents are edges to paint
	,pnt.vlty=list() #named list of character vectors where name is lty and contents are vertices to paint
	,pnt.lty=list() #named list of directed character edgelists where name is lty and contents are edges to paint
	,pnt.v2e=T
	,layout=list()
)
{
	require(network)

	el<-data.frame(s=as.character(el[,1]),r=as.character(el[,2]),stringsAsFactors=F)
	el<-el[order(el$s,el$r),]
	suppressWarnings(if(trim.pendants=="m1") {
		t<-table(el[[1]])
		elo<-el
		el<-el[el[[1]]%in%names(t)[t>1],]
	})
	suppressWarnings(if(trim.pendants=="m2") {
		t<-table(el[[2]])
		elo<-el
		el<-el[el[[2]]%in%names(t)[t>1],]
	})
	trim<-grepl("^m",trim.pendants[1])
	if(trim){
		om1n<-sort(unique(elo[[1]]))
		om2n<-sort(unique(elo[[2]]))
		om1l<-length(om1n)
		om2l<-length(om2n)
		(pbmnet<-network(elo,bipartite=om1l,directed=T,matrix.type="edgelist"))
	}

	### paint logic
	if(length(pnt.v)&pnt.v2e) for(i in 1:length(pnt.v)) {
		pnt.e[[length(pnt.e)+1]]<-t(apply(combn(
		setdiff(unlist(el[apply(apply(el,2,is.element,pnt.v[[i]]),1,any),]),pnt.v[[i]])
		,2),2,sort))
		names(pnt.e)[length(pnt.e)]<-names(pnt.v[i])
	}

	if(length(pnt.e)) for(i in 1:length(pnt.e)) colnames(pnt.e[[i]])<-c("s","r")

	if(length(pnt.vlty)&pnt.v2e) for(i in 1:length(pnt.vlty)) {
		pnt.lty[[length(pnt.lty)+1]]<-t(apply(combn(
		setdiff(unlist(el[apply(apply(el,2,is.element,pnt.vlty[[i]]),1,any),]),pnt.vlty[[i]])
		,2),2,sort))
		names(pnt.lty)[length(pnt.lty)]<-names(pnt.vlty[i])
	}
	if(length(pnt.lty)) for(i in 1:length(pnt.lty)) colnames(pnt.lty[[i]])<-c("s","r")

		print(pnt.v)
		print(pnt.e)
		print(pnt.vlty)
		print(pnt.lty)

	m1n<-sort(unique(el[[1]]))
	m2n<-sort(unique(el[[2]]))
	m1l<-length(m1n)
	m2l<-length(m2n)
	(bmnet<-network(el,bipartite=m1l,directed=T,matrix.type="edgelist"))

	m<-bmnet[,]
	w<-grepl(m1.stub,colnames(m))
	m1<-m[w,!w]%*%t(m[w,!w])
	m2<-t(m[w,!w])%*%m[w,!w]
	m1[lower.tri(m1)]<-0
	m2[lower.tri(m2)]<-0

	m1<-data.frame(s=rownames(m1)[which(!!m1,arr.ind=T)[,1]],r=rownames(m1)[which(!!m1,arr.ind=T)[,2]],ew=m1[which(!!m1)],stringsAsFactors=F)
	m2<-data.frame(s=rownames(m2)[which(!!m2,arr.ind=T)[,1]],r=rownames(m2)[which(!!m2,arr.ind=T)[,2]],ew=m2[which(!!m2)],stringsAsFactors=F)

	m1<-m1[order(m1$s,m1$r),]
	m2<-m2[order(m2$s,m2$r),]

	(m1net<-network(m1[,1:2],loops=T,directed=F,matrix.type="edgelist"))
	(m2net<-network(m2[,1:2],loops=T,directed=F,matrix.type="edgelist"))

	network.vertex.names(m1net)<-m1n
	network.vertex.names(m2net)<-m2n

	set.edge.attribute(m1net, "ew", m1$ew)
	set.edge.attribute(m2net, "ew", m2$ew)

	#mats<-matrix(nrow=m1l,ncol=m2l,dimnames=list(m1n,m2n))
	#system.time(for(i in 1:nrow(el)) mats[el[[1]],el[[2]]]<-1)
	#matr<-matrix(sample(c(0,0,0,0,1),replace=T,size=24),ncol=6,nrow=4)
	#mat<-matr
	#mat<-mats

	pdf(paste(out,"mode_projection.pdf",sep=""),h=(8.5-2)/3*ifelse(trim,2,1),w=(8.5-2)/3*ifelse(trim,2,3))
	if(trim)  par(mfrow=c(2,2),mar=rep(mar,4),omi=rep(omi,4)) else par(mfrow=c(1,3),mar=rep(mar,4),omi=rep(omi,4))
	if(trim){
		vcol<-c(rep(m1col,om1l),rep(m2col,om2l))
		el<-data.frame(as.edgelist(pbmnet))
		if(!is.directed(pbmnet)) el<-data.frame(t(apply(el,1,sort)),stringsAsFactors=F)
		ecolp<-rep(ecol,nrow(el))
		eltp<-rep(elt,nrow(el))
		vltp<-rep(vlt,network.size(pbmnet))

		if(length(pnt.v)) for(i in 1:length(pnt.v)) vcol[network.vertex.names(pbmnet)%in%pnt.v[[i]]]<-names(pnt.v[i])
		if(length(pnt.e)) for(i in 1:length(pnt.e)) ecolp[which(duplicated(rbind(do.call(cbind,el),pnt.e[[i]]),fromLast=T))]<-names(pnt.e[i])
		if(length(pnt.vlty)) for(i in 1:length(pnt.vlty)) vltp[network.vertex.names(pbmnet)%in%pnt.vlty[[i]]]<-names(pnt.vlty[i])
		if(length(pnt.lty)) for(i in 1:length(pnt.lty)) eltp[which(duplicated(data.frame(rbind(do.call(cbind,el),pnt.lty[[i]])),fromLast=T))]<-names(pnt.lty[i])

		lay<-list()
		for(i in names(layout)) lay[[i]]<-layout[[i]](x=network.size(pbmnet))

		plot(pbmnet
			,vertex.sides=c(rep(m1vsides,om1l),rep(m2vsides,om2l))
			,vertex.col=vcol
			,edge.col=ecolp
			,vertex.cex=vcx
			,vertex.lwd=vlw
			,edge.lwd=elw
			,edge.lty=eltp
			,vertex.lty=vltp
			,layout.par=lay
		)
		box()
		mtext(paste("Original",bmain),cex=cex)
	}

	vcol<-c(rep(m1col,m1l),rep(m2col,m2l))
	el<-data.frame(as.edgelist(bmnet))
	if(!is.directed(bmnet)) el<-data.frame(t(apply(el,1,sort)),stringsAsFactors=F)
	ecolp<-rep(ecol,nrow(el))
	eltp<-rep(elt,nrow(el))
	vltp<-rep(vlt,network.size(bmnet))

	if(length(pnt.v)) for(i in 1:length(pnt.v)) vcol[network.vertex.names(bmnet)%in%pnt.v[[i]]]<-names(pnt.v[i])
	if(length(pnt.e)) for(i in 1:length(pnt.e)) ecolp[which(duplicated(rbind(do.call(cbind,el),pnt.e[[i]]),fromLast=T))]<-names(pnt.e[i])
	if(length(pnt.vlty)) for(i in 1:length(pnt.vlty)) vltp[network.vertex.names(bmnet)%in%pnt.vlty[[i]]]<-names(pnt.vlty[i])
	if(length(pnt.lty)) for(i in 1:length(pnt.lty)) eltp[which(duplicated(data.frame(rbind(do.call(cbind,el),pnt.lty[[i]])),fromLast=T))]<-names(pnt.lty[i])

	lay<-list()
	for(i in names(layout)) lay[[i]]<-layout[[i]](x=network.size(bmnet))

	plot(bmnet
		,vertex.sides=c(rep(m1vsides,m1l),rep(m2vsides,m2l))
		,vertex.col=vcol
		,edge.col=ecolp
		,vertex.cex=vcx
		,vertex.lwd=vlw
		,edge.lwd=elw
		,edge.lty=eltp
		,vertex.lty=vltp
		,layout.par=lay
	)
	box()
	mtext(ifelse(trim,paste(bmain,"w/o Pendants"),bmain),cex=cex)

	vcol<-rep(m1col,m1l)
	el<-data.frame(as.edgelist(m1net))
	if(!is.directed(m1net)) el<-data.frame(t(apply(el,1,sort)),stringsAsFactors=F)
	ecolp<-rep(ecol,nrow(el))
	eltp<-rep(elt,nrow(el))
	vltp<-rep(vlt,network.size(m1net))

	if(length(pnt.v)) for(i in 1:length(pnt.v)) vcol[network.vertex.names(m1net)%in%pnt.v[[i]]]<-names(pnt.v[i])
	if(length(pnt.e)) for(i in 1:length(pnt.e)) ecolp[which(duplicated(rbind(do.call(cbind,el),pnt.e[[i]]),fromLast=T))]<-names(pnt.e[i])
	if(length(pnt.vlty)) for(i in 1:length(pnt.vlty)) vltp[network.vertex.names(m1net)%in%pnt.vlty[[i]]]<-names(pnt.vlty[i])
	if(length(pnt.lty)) for(i in 1:length(pnt.lty)) eltp[which(duplicated(data.frame(rbind(do.call(cbind,el),pnt.lty[[i]])),fromLast=T))]<-names(pnt.lty[i])

	lay<-list()
	for(i in names(layout)) lay[[i]]<-layout[[i]](x=network.size(m1net))

	plot(m1net
		,attrname="ew"
		,loop.cex=loopcx
		,vertex.sides=rep(m1vsides,m1l)
		,vertex.col=vcol
		,edge.col=ecolp
		,vertex.cex=vcx
		,vertex.lwd=vlw
		,edge.lwd=elw
		,edge.lty=eltp
		,vertex.lty=vltp
		,layout.par=lay
	)
	box()
	mtext(m1main,cex=cex)

	vcol<-rep(m2col,m2l)
	el<-data.frame(as.edgelist(m2net))
	if(!is.directed(m2net)) el<-data.frame(t(apply(el,1,sort)),stringsAsFactors=F)
	ecolp<-rep(ecol,nrow(el))
	eltp<-rep(elt,nrow(el))
	vltp<-rep(vlt,network.size(m2net))

	if(length(pnt.v)) for(i in 1:length(pnt.v)) vcol[network.vertex.names(m2net)%in%pnt.v[[i]]]<-names(pnt.v[i])
	if(length(pnt.e)) for(i in 1:length(pnt.e)) ecolp[which(duplicated(rbind(do.call(cbind,el),pnt.e[[i]]),fromLast=T))]<-names(pnt.e[i])
	if(length(pnt.vlty)) for(i in 1:length(pnt.vlty)) vltp[network.vertex.names(m2net)%in%pnt.vlty[[i]]]<-names(pnt.vlty[i])
	if(length(pnt.lty)) for(i in 1:length(pnt.lty)) eltp[which(duplicated(data.frame(rbind(do.call(cbind,el),pnt.lty[[i]])),fromLast=T))]<-names(pnt.lty[i])

	lay<-list()
	for(i in names(layout)) lay[[i]]<-layout[[i]](x=network.size(m2net))

	plot(
		m2net
		,attrname="ew"
		,loop.cex=loopcx
		,vertex.sides=rep(m2vsides,m2l)
		,edge.col=ecolp
		,vertex.col=vcol
		,vertex.cex=vcx
		,vertex.lwd=vlw
		,edge.lwd=elw
		,edge.lty=eltp
		,vertex.lty=vltp
		,layout.par=lay
	)
	box()
	mtext(m2main,cex=cex)
	dev.off()
	cat("\nPlot saved to",out)
	ret<-list(bmnet=bmnet,m1net=m1net,m2net=m2net)
	if(trim) ret<-c(pbmnet=list(pbmnet),ret)
	ret
}

thatgirlis.f<-function(
	n
	,ew="ew"
	,s=1000
	,plot=F
)
{
	require(network)
	ns<-network.size(n)
	edist<-table(n%e%ew)
	z<-(ns*(ns-1)/2)-sum(edist)
	o<-names(edist)
	edist<-c(z,edist)
	names(edist)<-c("0",o)
	edist<-as.table(edist)
	print(cbind(Freq=edist,Prop=round(edist/sum(edist),digits=4)))
	cat("Tot:",sum(edist),"\n\n")

	perm<-list()
	maxcombo<-ns*(ns-1)/2
	combos<-1:maxcombo
	choices<-sum(n%e%ew)
	t1<-Sys.time()
	for(i in 1:s){
		cat("\r",i,"\t",sep="")
		rdist<-sample(combos,size=choices,replace=T)
		rdist<-table(rdist)
		z<-maxcombo-length(rdist)
		rdist<-table(rdist)
		o<-names(rdist)
		rdist<-c(z,rdist)
		names(rdist)<-c("0",o)
		rdist<-as.table(rdist)
		perm[[i]]<-rdist
	}
	t2<-Sys.time()
	cat(":",round((t2-t1)/60,2),"minutes to permute\n")
	maxcount<-max(sapply(perm,length))
	for(i in 1:length(perm)) perm[[i]]<-c(perm[[i]],rep(0,maxcount-length(perm[[i]])))
	perm<-do.call(rbind,perm)
	colnames(perm)<-0:(maxcount-1)
	if(dim(perm)[2]>4) {perm<-cbind(perm,apply(perm[,4:dim(perm)[2]],1,sum));colnames(perm)[dim(perm)[2]]<-">=3"}

	if(plot) for(i in which(edist>0)) {hist(perm[,i],breaks=(floor(min(perm[,i]))-.5):(ceiling(max(perm[,i]))+.5),freq=F,main=paste("Count:",i),xlab="",xlim=range(c(perm[,i],edist[i])));abline(v=edist[i],lty=2)}

	edist<-c(edist,rep(0,maxcount-length(edist)))
	names(edist)<-0:(maxcount-1)
	if(length(edist)>4) {edist<-c(edist,sum(edist[4:length(edist)]));names(edist)[length(edist)]<-">=3"}
	flush.console()
	edist<-cbind(observed=edist,expected=round(apply(perm,2,mean),1),sd=round(apply(perm,2,sd),3),t=round((edist-apply(perm,2,mean))/apply(perm,2,sd),3),"p o<=e"=round(apply(perm<=edist,2,mean),4),"p o>=e"=round(apply(perm>=edist,2,mean),4),t(round(apply(perm,2,quantile,prob=(c(seq(0,1,.1),.25,.75))),3)))
	edist
}

plotpois<-function(
	pois
	,year
	,jour
	,count
	,q1="0%"
	,q2="50%"
	,q3="100%"
)
{
	xlim<-range(year)
	plt<-do.call(c,pois[as.character(year),jour])
	plt<-plt[!is.na(plt)]
	for(i in 1:length(plt)) plt[[i]]<-cbind(year=as.numeric(names(plt[i])),plt[[i]])
	plt<-do.call(rbind,plt)
	plt<-plt[rownames(plt)==as.character(count),]

	ylim=range(plt[,c("observed",q1,q2,q3)])
	plot.new()
	plot.window(xlim=xlim,ylim=ylim)
	axis(side=1,lab=as.character(xlim[1]:xlim[2]),at=xlim[1]:xlim[2])
	axis(side=2)
	title(main=jour,xlab="Year",ylab=paste("Count of ",count,"'s",sep=""))
	w<-which(!!diff(plt[,"year"])-1)
	b<-1
	for(i in unique(c(w,dim(plt)[1]))){
		if(!length(w)) {r<-1:dim(plt)[1]} else {r<-b:i;b<-i+1}
		lines(plt[r,"year"],y=plt[r,q2],col="red",lty=1)
		lines(plt[r,"year"],plt[r,q1],lty=3,col="red")
		lines(plt[r,"year"],plt[r,q3],lty=3,col="red")
		lines(plt[r,"year"],plt[r,"observed"],lty=1,lwd=1)
		if(!length(w)) break
	}
}

subnet<-function(
	db2bel=stop("Supply original db2bel object",call.=F)
	,set=stop("set list(cr=incl cr,ut=incl ut,ncr=excl cr,nut=excl ut) and unused to NULL",call.=F)
	,source=stop("Supply source",call.=F)
)
{
	require(network)
	source(source)
	#if(!all(c("bel2mel.f","mel2net.f")%in%ls())) stop("Load correct source")
	#w<-!!sapply(lapply(db2bel$bel,"%in%",set),any)
	#if(!sum(w)) stop("Not a subset of either mode of this edgelist")
	#s<-set%in%db2bel$bel[,w]
	#if(!all(s)) {
	#	warning(paste(sum(!s),"or",round(sum(!s)/length(s)*100,1),"% of edges from subset are not in bel. First 10 excluded:"))
	#	print(head(set[s],10))
	#}
	sub<-rep(T,dim(db2bel$bel)[1])
	if(length(set$cr)) sub<-sub&db2bel$bel$cr%in%set$cr
	if(length(set$ut)) sub<-sub&db2bel$bel$ut%in%set$ut
	if(length(set$ncr)) sub<-sub&!db2bel$bel$cr%in%set$ncr
	if(length(set$nut)) sub<-sub&!db2bel$bel$ut%in%set$nut
	db2bel$bel<-db2bel$bel[sub,]
	if("pend"%in%names(db2bel)) db2bel$pend<-db2bel$pend[sub]
	cat("\n",nrow(db2bel$bel)," edges, ",length(unique(db2bel$bel$ut))," uts, and ",length(unique(db2bel$bel$cr)) ," crs fed to bel2mel.\n",sep="")
	bel2mel<-bel2mel.f(db2bel,out=getwd())
	mel2net<-mel2net.f(bel2mel)
	mel2net
}

dbl2w.f<-function(
	wok2db
	,out=stop("Specify output directory")
	,field=stop("field=c(\"field1\",\"field2\",...)")
	,variations=NULL
	,recode=NULL
)
{
require(data.table)
wok2db<-data.table(wok2db)
setnames(wok2db,c("ut","field","b"))
setkey(wok2db,field)
wok2db<-wok2db[field]
setkey(wok2db,field)

field<-tolower(field)
field<-field[field!="ut"]

l<-list()
for(i in field){
	l[[i]]<-wok2db[i=toupper(i)][j=list(ut,b)];setkey(l[[i]],ut);setnames(l[[i]],2,i)
	if(i%in%c("cr","af")) l[[i]][,i:=factor(toupper(sub(", DOI .+","",l[[i]][[i]]))),with=F]
	if(i=="af") l[[i]]<-l[[i]][i=!grepl("ANONYMOUS",l[[i]][[i]])]
	l[[i]][,c(i):=type.convert(as.character(l[[i]][[i]]))]
}
rm(wok2db)

if(is.null(variations)) {
	cat("\nvariations=list(field1=db2bel_sets1,field2=db2bel_sets2,...)")
}
else
{
	pre<-sort(unlist(c(0:9,strsplit("! \" # $ % & ' ( ) * + , - . / : ; < = > ? @ [ \\ ] ^ _ ` { | } ~"," "),letters,LETTERS)))[1]
	for(i in tolower(names(variations))) for(j in 1:length(variations[[i]])) levels(l[[i]][[i]])[levels(l[[i]][[i]])%in%variations[[i]][[j]]]<-paste(pre,"z",toupper(i),formatC(j,width=nchar(as.character(length(variations[[i]]))),format="d",flag="0"),sep="")
}

if(is.null(recode)) {
	cat("\nrecode=list(\n\tfield1=list(\n\t\t\"recode1\"=c(\"code1\",\"code2\",...)\n\t\t,\"recode2\"=c(\"code1\",\"code2\",...)\n\t)\n\t,field2=...\n)")
}
else
{
	for(i in tolower(names(recode))) for(j in names(recode[[i]]))  levels(l[[i]][[i]])[levels(l[[i]][[i]])%in%recode[[i]][[j]]]<-j
}

if("cr"%in%field) {
	npen<-l$cr[,.N,by="cr"]
	npen<-as.character(npen$cr[npen$N>1])
	l$nr<-l$cr[,.N,by="ut"]
	setnames(l$nr,c("ut","nr"))
	setkey(l$cr,"cr")
	l$nrtp<-l$cr[npen,.N,keyby="ut"]
	setnames(l$nrtp,c("ut","nrtp"))
	setkey(l$cr,"ut")
}

if("so"%in%field) {
	### improve later for custom coding of natural text
	rc<-data.frame(c=c("soci","[ck]ono","anth","poli"),r=c("Sociology","Economics","Anthropology","Political Science"))
	l$ds<-copy(l$so)
	setnames(l$ds,2,"ds")
	for(i in 1:nrow(rc)) l$ds[grep(rc$c[i],l$so$so,ignore.case=T),ds:=rc$r[i]]
	l$ds$ds<-factor(l$ds$ds)
}

### merge
dbl2w<-copy(l[[1]])
for(i in names(l)[-1]) dbl2w<-merge(dbl2w,l[[i]],all=T,allow.cartesian=TRUE)
rm(l)

### code selection effect
if("cr"%in%field) {
	dbl2w[is.na(dbl2w$nr),nr:=0]
	dbl2w[is.na(dbl2w$nrtp),nrtp:=0]
	dbl2w[,sel:=as.integer(nrtp>0)]
	dbl2w[,rej1:=as.integer(nr==0)] #document rejected if no citations
	dbl2w[,rej2:=as.integer(nrtp==0&nr!=0)] #document rejected if no citations after selection
}

w<-unique(c("ut",field[field%in%c("af","cr")]))
lvs<-lapply(as.list(1:length(w)),FUN=function(x) apply(combn(w,x),2,paste,sep=""))
lvs[[1]]<-matrix(lvs[[1]],ncol=length(w))
for(i in 1:length(lvs)) for(j in 1:ncol(lvs[[i]])) {
	levs<-lvs[[i]][,j]
	setkeyv(dbl2w,levs)

	samp<-list()
	for(k in levs) samp[[k]]<-!is.na(dbl2w[[k]])
	cat("Proportion missing:\n")
	print(round(sapply(samp,FUN=function(x) sum(!x))/nrow(dbl2w),3))
	samp<-do.call(cbind,samp)
	if(ncol(samp)>1) for(i in 2:ncol(samp)) samp[,1]<-samp[,1]&samp[,i]
	samp<-as.vector(samp[,1])
	tp<-!!dbl2w$sel

	nm<-paste(c("w",levs),collapse="")

	dbl2w<-dbl2w[i=samp,1/.N,keyby=c(levs),][dbl2w]
	dbl2w[(!samp)|is.na(dbl2w$V1),V1:=0]
	setnames(dbl2w,"V1",nm)

	setkeyv(dbl2w,levs)
	dbl2w<-dbl2w[i=tp&samp,1/.N,keyby=c(levs)][dbl2w]
	dbl2w[(!tp)|is.na(dbl2w$V1),V1:=0]
	setnames(dbl2w,"V1",paste(nm,"tp",sep=""))
	cat(nm,"\n")
	rep<-rbind(c=c(nrow(dbl2w[i=samp,NA,keyby=levs]),nrow(dbl2w[i=(!!dbl2w$sel)&samp,NA,keyby=levs])),w=c(sum(dbl2w[[nm]]),sum(dbl2w[[paste(nm,"tp",sep="")]])))
	colnames(rep)<-c("unsel","sel")
	print(rep)
}
w<-unlist(sapply(lvs,FUN=function(x) apply(x,2,FUN=function(y) paste(c("w",y),sep="",collapse=""))))
w<-c(w,paste(w,"tp",sep=""))
dbl2w<-dbl2w[j=order(names(dbl2w)%in%w),with=F]
if(!is.null(out)) save(dbl2w,file=paste(out,"dbl2w.RData",sep=.Platform$file.sep))
setkeyv(dbl2w,levs)
dbl2w
}

w2tab.f<-function(
	dbl2w
	,out=NULL
	,hg #higher (smaller) group classification
	,lg #lower (larger) group classification
	,sort #level by which all tables should be sorted in descending order
	,s.e=F
	,reps=10
	,decreasing=T
	,addlev=NULL
)
{
	require(data.table)
	w<-grep("w",names(dbl2w),value=T) #automatically processes any "weight" variables with a w- prefix
	nw<-nchar(w)
	lvar<-w[which.max(nw)]
	lvar<-gsub("(^w)|(tp$)","",lvar)
	nw<-nchar(lvar)
	lvar<-mapply(substr,lvar,seq(1,nw,2),seq(2,nw,2),USE.NAMES=F)
	lvar<-unique(c(lvar,addlev))
	dbl2w<-dbl2w[j=c(lvar,hg,lg,w),with=F]
	setkeyv(x=dbl2w,cols=lvar)

	t2prop<-function(x) if(!(is.character(x)|is.factor(x))) {prop.table(x)*100} else {x}

	tl.f<-function(dbl2w,wvar,levs,hg,lg,sort=c("key","order","given"),decreasing=T,subset=NULL){
	if(!is.null(subset)) {
		setkeyv(dbl2w,levs)
		dbl2w<-dbl2w[i=subset]
	}
	tl<-dbl2w[j=list(hc=sum(get(wvar))),keyby=hg]
	tl<-tl[dbl2w[j=list(lc=sum(get(wvar))),keyby=c(hg,lg)]]
	setkeyv(tl,c(hg,lg))
	if(sort[1]=="order") {o<-order(tl[[2]],tl[[4]],decreasing=decreasing);tl<-list(tl=tl,o=o)}
	if(sort[1]=="given") tl<-tl[i=sort,]
	tl
	}

	init<-tl.f(dbl2w,wvar=sort,sort="order",hg=hg,lg=lg)
	tl<-list(init$tl)
	names(tl)<-sort

	wo<-w[w!=sort]

	for(i in wo) tl[[i]]<-tl.f(dbl2w,wvar=i,hg=hg,lg=lg)
	tl<-tl[w]
	dup<-!duplicated(tl[[sort]][init$o,get(hg)])
	lh<-cumsum(dup)+1:nrow(tl[[sort]])
	h<-which(!(1:max(lh))%in%lh)
	hlo<-order(c(h,lh))
	ow<-grep("tp$",w,value=T,invert=T)
	if(s.e){
		se<-list()
		for(i in ow){
			wvar<-sub("tp$","",i)
			samplev<-unlist(sapply(as.list(addlev),FUN=function(x) sub(x,"",grep(x,wvar,value=T))))
			if(length(samplev)) {wvartp<-paste(samplev,"tp",sep="")} else {wvartp<-paste(wvar,"tp",sep="")}
			lvar<-gsub("(^w)|(tp$)","",wvartp) # level varnams cannot have w or tp in them, and must be exactly two characters, e.g. ut, af, cr
			levs<-mapply(substr,lvar,seq(1,nchar(lvar),2),seq(2,nchar(lvar),2),USE.NAMES=F)
			samp<-list()
			for(j in levs) samp[[j]]<-dbl2w[j=!is.na(get(j))]
			samp<-do.call(cbind,samp)
			if(ncol(samp)>1) for(j in 2:ncol(samp)) samp[,1]<-samp[,1]&samp[,j]
			samp<-as.vector(samp[,1])
			samp<-dbl2w[i=samp,NA,keyby=levs][,levs,with=F]
			cat(i,"\n")
			se[[i]]<-replicate(reps
				,tl.f(dbl2w,wvar=wvar
					,levs=levs
					,subset=samp[i=sample(1:nrow(samp),sum(dbl2w[[wvartp]]))]
					,hg=hg
					,lg=lg
				)
			,simplify=F)

			J<-se[[i]][[1]][j=c(1,3,2,4),with=F][tl[[i]][j=c(1,3),with=F]]
			for(j in 2:length(se[[i]])) J[,paste(c("hc","lc"),j,sep=""):=as.list(se[[i]][[j]][j=c(1,3,2,4),with=F][tl[[i]][j=c(1,3,2,4),with=F]][,c("hc","lc"),with=F])]
			J[is.na(J)]<-0
			J<-J[init$o,]
			se[[i]]<-cbind(
				unlist(c(J[dup,j=hg,with=F],J[j=lg,with=F]))
				,rbindlist(list(
					J[i=dup,j=grep("hc",names(J),value=),with=F]
					,J[j=grep("lc",names(J),value=),with=F]))
				)[hlo]


			se[[i]]<-list(f=se[[i]],p=copy(se[[i]]))
			col<-colnames(se[[i]]$p)[-1]
			se[[i]]$p[h,col:=lapply(se[[i]]$p[h,col,with=F],t2prop),with=F]
			se[[i]]$p[!h,col:=lapply(se[[i]]$p[!h,col,with=F],t2prop),with=F]
			ssr<-function(x) sum(round(x,4)^2)
			sr<-function(x) sum(round(x,4))
			if(length(h)==nrow(se[[i]]$f)) {
				se[[i]]$p<-data.frame(level=c(as.character(se[[i]]$p[[1]]),"High Total","High H Index"),rbind(
					se[[i]]$p[h,!1,with=F]
					,se[[i]]$p[h,lapply(.SD, sr),.SDcols=-1]
					,se[[i]]$p[h,lapply(.SD, ssr),.SDcols=-1]
				))
				se[[i]]$f<-data.frame(level=c(as.character(se[[i]]$f[[1]]),"High Total"),rbind(
					se[[i]]$f[h,!1,with=F]
					,se[[i]]$f[h,lapply(.SD, sr),.SDcols=-1]
				))
			} else {
				se[[i]]$p<-data.frame(level=c(as.character(se[[i]]$p[[1]]),"High Total","Low Total","High H Index","Low H Index"),rbind(
					se[[i]]$p[,!1,with=F]
					,se[[i]]$p[h,lapply(.SD, sr),.SDcols=-1]
					,se[[i]]$p[!h,lapply(.SD, sr),.SDcols=-1]
					,se[[i]]$p[h,lapply(.SD, ssr),.SDcols=-1]
					,se[[i]]$p[!h,lapply(.SD, ssr),.SDcols=-1]
				))
				se[[i]]$f<-data.frame(level=c(as.character(se[[i]]$f[[1]]),"High Total","Low Total"),rbind(
					se[[i]]$f[,!1,with=F]
					,se[[i]]$f[h,lapply(.SD, sr),.SDcols=-1]
					,se[[i]]$f[!h,lapply(.SD, sr),.SDcols=-1]
				))
			}
			se[[i]]$sht.p<-try(do.call(rbind,apply(se[[i]]$p[,-1],1,FUN=function(x) shapiro.test(x)[1:2]))) #shapiro test of normality of p table
			for(j in c("f","p")) {se[[i]][[j]]<-data.frame(level=se[[i]][[j]][[1]],s.e.=apply(se[[i]][[j]][,-1],1,sd));se[[i]][[j]][[2]]<-round(se[[i]][[j]][[2]],4)}
		}
	}

	t2prop<-function(x) if(!(is.character(x)|is.factor(x))) {round(prop.table(x)*100,3)} else {x}
	for(j in names(tl)) {
		tl[[j]]<-tl[[j]][init$o]
		tl[[j]]<-cbind(
				unlist(c(tl[[j]][dup,j=hg,with=F],tl[[j]][j=lg,with=F]))
				,rbindlist(list(
					tl[[j]][i=dup,j=grep("hc",names(tl[[j]]),value=),with=F]
					,tl[[j]][j=grep("lc",names(tl[[j]]),value=),with=F]))
				)[hlo]

				tl[[j]]<-list(f=tl[[j]],p=copy(tl[[j]]))
				col<-colnames(tl[[j]]$p)[-1]
				tl[[j]]$p[h,col:=lapply(tl[[j]]$p[h,col,with=F],t2prop),with=F]
				tl[[j]]$p[!h,col:=lapply(tl[[j]]$p[!h,col,with=F],t2prop),with=F]

				ssr<-function(x) sum(round(x,3)^2)
				sr<-function(x) sum(round(x,3))
				if(length(h)==nrow(tl[[j]]$f)) {
					tl[[j]]$p<-data.frame(level=c(as.character(tl[[j]]$p[[1]]),"High Total","High H Index"),rbind(
						tl[[j]]$p[h,!1,with=F]
						,tl[[j]]$p[h,lapply(.SD, sr),.SDcols=-1]
						,tl[[j]]$p[h,lapply(.SD, ssr),.SDcols=-1]
					))
					tl[[j]]$f<-data.frame(level=c(as.character(tl[[j]]$f[[1]]),"High Total"),rbind(
						tl[[j]]$f[h,!1,with=F]
						,tl[[j]]$f[h,lapply(.SD, sr),.SDcols=-1]
					))
				} else {
					tl[[j]]$p<-data.frame(level=c(as.character(tl[[j]]$p[[1]]),"High Total","Low Total","High H Index","Low H Index"),rbind(
						tl[[j]]$p[,!1,with=F]
						,tl[[j]]$p[h,lapply(.SD, sr),.SDcols=-1]
						,tl[[j]]$p[!h,lapply(.SD, sr),.SDcols=-1]
						,tl[[j]]$p[h,lapply(.SD, ssr),.SDcols=-1]
						,tl[[j]]$p[!h,lapply(.SD, ssr),.SDcols=-1]
					))
					tl[[j]]$f<-data.frame(level=c(as.character(tl[[j]]$f[[1]]),"High Total","Low Total"),rbind(
						tl[[j]]$f[,!1,with=F]
						,tl[[j]]$f[h,lapply(.SD, sr),.SDcols=-1]
						,tl[[j]]$f[!h,lapply(.SD, sr),.SDcols=-1]
					))
				}

		}

	tab<-list()
	ow<-names(sort(dbl2w[,lapply(.SD,sum),.SDcols=ow]))
	for(j in ow){
		tab[[j]]<-data.frame(tl[[j]]$p,tl[[paste(j,"tp",sep="")]]$p)[,-3]
		tab[[j]]<-data.frame(tab[[j]],apply(tab[[j]][,-1],1,diff))
		if(s.e) tab[[j]]<-data.frame(tab[[j]],se[[j]]$p[[2]],tab[[j]][[4]]/se[[j]]$p[[2]])
		tab[[j]][,-1]<-round(tab[[j]][,-1],3)
		if(s.e) {colnames(tab[[j]])<-c("l","o","s","Δ","se","t")} else {colnames(tab[[j]])<-c("l","o","s","Δ")}
		wt<-grep("Total",tl[[j]]$f$level)
		tot<-tl[[paste(j,"tp",sep="")]]$f[wt,]
		tot[[1]]<-sub("Total","N",tot[[1]])
		if(s.e) {tot<-data.frame(tot[[1]],tl[[j]]$f[wt,2],tot[[2]],round((tot[[2]]-tl[[j]]$f[wt,2])/tl[[j]]$f[wt,2]*100,3),rep(NA,length(wt)),rep(NA,length(wt)))
		} else{tot<-data.frame(tot[[1]],tl[[j]]$f[wt,2],tot[[2]],round((tot[[2]]-tl[[j]]$f[wt,2])/tl[[j]]$f[wt,2]*100,3))}
		if(s.e) {colnames(tab[[j]])<-c("l","o","s","Δ","se","t")} else {colnames(tab[[j]])<-c("l","o","s","Δ")}
		if(s.e) {colnames(tot)<-c("l","o","s","Δ","se","t")} else {colnames(tot)<-c("l","o","s","Δ")}
		tab[[j]]<-rbind(tab[[j]],tot)
	}
	tab<-tab[ow]
	tab<-data.frame(tab[[1]],lapply(tab[-1],FUN=function(x) {x[[1]]<-NULL;x}))
	cn<-grep("\\.",names(tab),value=T,invert=T)
	colnames(tab)[colnames(tab)%in%cn]<-paste(ow[1],cn,sep=".")
	colnames(tab)<-sub("^w","",colnames(tab))
	colnames(tab)[1]<-"level"
	ro<-grepl("\\.o$",colnames(tab))
	tab<-data.frame(tab,replicate(sum(ro),rep(NA,nrow(tab)),simplify=F))[,order(c(1:ncol(tab),which(ro)-.5))]
	ro<-grepl("\\.o$",colnames(tab))
	colnames(tab)[which(ro)-1]<-""
	tab[grep("Total",tab[,1]),grep("(\\.se$)|(\\.t$)",colnames(tab))]<-NA
	if(!is.null(out)) write.table(tab,file=paste(out,"selection_table.tab",sep=""),sep="\t",quote=F,row.names=F,na="")
	w2tab
}

.ls.objects <- function (
	pos = 1
	, pattern
	, order.by
	, decreasing=FALSE
	, head=FALSE
	, n=5
)
{
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.prettysize <- napply(names, function(x) {
                           capture.output(print(object.size(x), units = "auto")) })
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
    names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
# shorthand
lsos <- function(
	...
	, n=10
)
{
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

if(F){
library(linkcomm)
load("bel2mel.RData")
n<-names(bel2mel)
n<-n[grep("el$",n)]
mel2lc<-list()
for(i in n){
	pdf(paste("mel2lc_",i,".pdf",sep=""))
	mel2lc[[i]]<-getLinkCommunities(bel2mel[[i]][,-(3:4)],removetrivial=F)
	dev.off()
	save(mel2lc,file="mel2lc_rt-f_uw.RData")
}
}

rect.dendrogram<-function (tree, k = NULL, which = NULL, x = NULL, h = NULL, border = 2,
													 cluster = NULL, horiz = FALSE, density = NULL, angle = 45,
													 text = NULL, text_cex = 1, text_col = 1, xpd = TRUE, lower_rect, hpk=NULL,
													 ...
)
{
	if (!is.dendrogram(tree))
		stop("x is not a dendrogram object.")
	if (length(h) > 1L | length(k) > 1L)
		stop("'k' and 'h' must be a scalar(i.e.: of length 1)")
	if(is.null(hpk)) {tree_heights <- heights_per_k.dendrogram(tree)[-1]} else {tree_heights<-hpk[-1]}
	tree_order <- order.dendrogram(tree)
	if (!is.null(h)) {
		if (!is.null(k))
			stop("specify exactly one of 'k' and 'h'")
		ss_ks <- tree_heights < h
		k <- min(as.numeric(names(ss_ks))[ss_ks])
		k <- max(k, 2)
	}
	else if (is.null(k))
		stop("specify exactly one of 'k' and 'h'")
	if (k < 2 | k > length(tree_heights))
		stop(gettextf("k must be between 2 and %d", length(tree_heights)),
				 domain = NA)
	if (is.null(cluster))
		cluster <- cutree(tree, k = k)
	clustab <- table(cluster)[unique(cluster[tree_order])]
	m <- c(0, cumsum(clustab))
	if (!is.null(x)) {
		if (!is.null(which))
			stop("specify exactly one of 'which' and 'x'")
		which <- x
		for (n in seq_along(x)) which[n] <- max(which(m < x[n]))
	}
	else if (is.null(which))
		which <- 1L:k
	if (any(which > k))
		stop(gettextf("all elements of 'which' must be between 1 and %d",
									k), domain = NA)
	border <- rep_len(border, length(which))
	retval <- list()
	old_xpd <- par()["xpd"]
	par(xpd = xpd)
	for (n in seq_along(which)) {
		if (!horiz) {
			xleft = m[which[n]] + 0.66
			if (missing(lower_rect))
				lower_rect <- par("usr")[3L] - strheight("W") *
				(max(nchar(labels(tree))) + 1)
			ybottom = lower_rect
			xright = m[which[n] + 1] + 0.33
			ytop = tree_heights[names(tree_heights) == k]
		}
		else {
			ybottom = m[which[n]] + 0.66
			if (missing(lower_rect))
				lower_rect <- par("usr")[2L] + strwidth("X") *
				(max(nchar(labels(tree))) + 1)
			xright = lower_rect
			ytop = m[which[n] + 1] + 0.33
			xleft = tree_heights[names(tree_heights) == k]
		}
		rect(xleft, ybottom, xright, ytop, border = border[n],
				 density = density, angle = angle, ...)
		if (!is.null(text))
			text((m[which[n]] + m[which[n] + 1] + 1)/2, grconvertY(grconvertY(par("usr")[3L],
																																				"user", "ndc") + 0.02, "ndc", "user"), text[n],
					 cex = text_cex, col = text_col)
		retval[[n]] <- which(cluster == as.integer(names(clustab)[which[n]]))
	}
	par(xpd = old_xpd)
	invisible(retval)
}

strdist.dend.picker<-function(
	hd # must be dendrogram
	,s1=c(0,.075)
	,s2=c(0,1)
	,lightest.val=c(">"=.5,"="=.8)
	,maxlist=50
	,must.be.this.short.to.ride=.075
	,out=NULL
	,instruct=F
)
{

	require(dendextend)
	require(data.table)

	lintran<-function(x,s1=c(0,1),s2=c(0,1)) {a=diff(s2)/diff(s1);b=s2[1]-a*s1[1];return(a*x+b)}
	hs<-get_nodes_attr(hd,"height",simplify=T)
	l<-which(hs==0)
	n<-which(hs!=0)
	oc<-sapply(as.list(l),function(x) hs[n[which.max(which(n<x))]])
	names(oc)<-labels(hd)
	c<-lintran(oc,s1=s1,s2=s2)
	c[c>lightest.val[">"]]<-lightest.val["="]
	labels_colors(hd)<-gray(c)

	hk<-heights_per_k.dendrogram(hd)

	for(i in 1:length(hk)) if(max(table(cutree_1k.dendrogram(hd,k=as.integer(names(hk[i])),use_labels_not_values=FALSE,dend_heights_per_k=hk)))<=maxlist) break
	#if(i==1) maxlist<-as.integer(names(hk[i]))
	hdc<-cutree_1k.dendrogram(hd,k=as.integer(names(hk[i])),use_labels_not_values=T,dend_heights_per_k=hk)
	if(length(hdc)==2) hdc[2]<-1
	j<-i-1
	hkh<-hk[ifelse(j!=0,j,1)]

	df<-data.frame(get_nodes_xy(hd),get_nodes_attr(hd,"label"))
	names(df)<-c("x","y","lab")
	df<-df[!is.na(df$lab),]
	rownames(df)<-df$lab
	df[names(oc),"lh"]<-oc
	df[names(hdc),"hdc"]<-hdc
	rownames(df)<-NULL

	df<-data.table(df,key="lab")

	hdc<-split(hdc,hdc)
	hd<-color_branches(hd,k=length(hdc))
	lhdc<-length(hdc)

	graphics.off()
	sets<-list()
	if(instruct) cat("\rLeft-click to left of name to choose members of a set.\nWhen set is finished, right-click once to start a new set, or twice to advance the page.\nIf you want to accept the entire page, left-click five times, then right click twice.\nIf you make a mistake before starting a new set, select the same name three times, right-click to start a new set, and begin again.")
	c<-0

	for(i in 1:length(hdc)){
		sub<-names(hdc[[i]])
		if(min(df[sub,y])>must.be.this.short.to.ride) next
		yl<-range(df[sub,x],max(df[sub,x])-maxlist)
		#cx<-strheight(sub[1])*1.5*length(sub)/diff(yl)*10
		cx<-.7
		par(mai=c(.75,0,0,strwidth(sub[which.max(nchar(sub))],units="inches")*cx+1),cex=cx)
		plot(hd,"triangle",horiz=T,ylim=yl)
		te<-try(rect.dendrogram(hd,k=lhdc,which=lhdc-i+1,horiz=T,border="red",lty="dotted",hpk=hk),silent=T)
		if(!inherits(te,"try-error")) rect.dendrogram(hd,k=lhdc,which=lhdc-i+1,horiz=T,border="red",lty="dotted",hpk=hk)
		loc<-NA
		e<-1
		while(!is.null(loc)) {
			ind<-length(sets)+1
			nind<-rev(unlist(strsplit(as.character(ind),"")))[1]
			loc<-locator(type="p",pch=nind,col="red",cex=.5)$y
			if(is.null(loc)) break
			sets[[ind]]<-labels(hd)[round(loc)]
			if(any(table(sets[[ind]])%in%3:4)) sets[[ind]]<-sub
			segments(x0=0,y0=which(labels(hd)%in%sets[[ind]])
							 ,x1=0 #max(strwidth(sets[[ind]]))
							 ,lwd=10,col=hsv(.1*e,1,1,.5))
			e<-e+1
		}
		if(!is.null(out)) {
			c<-c+1
			dev.print(pdf,paste(out,c,file="picker_log.pdf",sep=""))
		}
	}
	sets[sapply(sets,function(x) any(table(x)>=5))]<-NULL
	sets<-lapply(sets,unique)
	sets[sapply(sets,function(x) any(sapply(setdiff(sets,list(x)),function(y) all(x%in%y))))]<-NULL
	if(!is.null(out)) save(sets,file=paste(out,file="picker_sets.RData",sep=""))
	sets
}


write.ergmm<-function(
	where=stop("\"where\" = folder filepath for scripts",call.=F)
	,dat=stop("\"dat\" = name of .RData file containing (stat)nets, w. ext",call.=F)
	,net=stop("\"net\" = named list of (stat)nets",call.=F)
	,dimensions=2
	,groups=1
	,verbosity=2
	,mcmc.size=10000
	,mcmc.burn=30000
	,mcmc.inter=10
)
{
	mod<-paste("d",dimensions,"G",groups,sep="")
	writeLines(c(
		"getwd()"
		,"rm(list=ls())"
		,"library(statnet)"
		,"library(latentnet)"
		,paste("dat<-\"",dat,"\"",sep="")
		,paste("load(dat)",sep="")
		,"hoff2fit<-list()"
		,"time<-round(unclass(Sys.time()))"
		,paste("for(i in rev(names(",net,"))){",sep="")
		,"cat(\"\\n<<<<<<<<< \",i,\" >>>>>>>>>\\n\\n\",sep=\"\")"
		,paste("hoff2fit[[paste(i,\"",mod,"\",sep=\"\")]]<-try(ergmm(",net,"[[i]]~euclidean(d=",dimensions,",G=",groups,"),verbose=",verbosity,",control=control.ergmm(\nsample.size=",mcmc.size,"\n,burnin=",mcmc.burn,"\n,interval=",mcmc.inter,"\n)))",sep="")
		,paste("save(hoff2fit,file=paste(\"hoff2fit_\",sub(\".RData\",\"\",dat),\"_\",\"",mod,"\",\"_\",time,\".RData\",sep=\"\"))",sep="")
		,"}"
	),con=paste(where,.Platform$file.sep,"write.ergmm_",sub(".RData","",dat),"_",mod,".R",sep=""))
}

cull.f<-function(
	a
	,b,cull=.2
	,noanon=F)
{
	if(any(is.na(b))) {b<-na.omit(b);attributes(b)<-NULL}
	if(noanon) {
		if(grepl("\\[ANONYMOUS\\],? ?",a)) stop("k = \"",a,"\"\n",call.=F)
		b<-sub("\\[ANONYMOUS\\],? ?","",b)
	}
	b<-b[nchar(b)>13]
	if(!is.null(cull)) b<-b[stringdist(a=a,b=b,method="jw",p=.1)<cull]
	b
}

unicr<-function(
	index
)
{
	yu<-list()
	cbeg<-list()
	cend<-list()
	cat("    ")
	for(i in 1:dim(index)[1]){
		d<-dimnames(index)[[1]][i]
		cat("\b\b\b\b",d)
		flush.console()
		yu[[d]]<-length(levels(droplevels(comdb[J(unlist(index[i,]),"CR")]$b)))
		cbeg[[d]]<-length(levels(droplevels(comdb[J(unlist(index[1:i,]),"CR")]$b)))
		cend[[d]]<-length(levels(droplevels(comdb[J(unlist(index[i:dim(index)[1],]),"CR")]$b)))
	}
	u<-data.frame(yu=unlist(yu),cbeg=unlist(cbeg),cend=unlist(cend))
	u
}

poisson_permute<-function()
{
	#permute the feckin distribution
	perm<-list()
	maxcombo<-1216*1215/2
	combos<-1:maxcombo
	choices<-sum(lc2mel$ew)
	t1<-Sys.time()
	for(i in 1:1000){
		cat("\r",i,sep="")
		rdist<-sample(combos,size=choices,replace=T)
		rdist<-table(rdist)
		z<-maxcombo-length(rdist)
		rdist<-table(rdist)
		n<-names(rdist)
		rdist<-c(z,rdist)
		names(rdist)<-c("0",n)
		rdist<-as.table(rdist)
		perm[[i]]<-rdist
	}
	t2<-Sys.time()
	cat(t2-t1,"minutes")

	cn<-sort(unique(unlist(lapply(perm,names))))
	permdb<-data.frame(matrix(0,nrow=length(perm),ncol=length(cn)))
	permdb<-data.frame(permdb,matrix(0,nrow=length(perm),ncol=length(edist)-length(cn)))
	names(permdb)<-names(edist)
	for(i in 1:length(perm)) permdb[i,names(perm[[i]])]<-perm[[i]]

	dif<-matrix(edist,nrow=dim(permdb)[1],ncol=dim(permdb)[2],byrow=T)-permdb
	md<-apply(dif,2,mean)
	sdd<-apply(dif,2,sd)
	cid<-apply(dif,2,quantile,prob=c(.05,.95))
	tad<-apply(dif,2,table)
	num<-lapply(tad,names)
	num<-lapply(num,as.numeric)

	dens<-1-(permdb$`0`/maxcombo)
	mean((1-edist["0"]/maxcombo)-dens)*100
	sd((1-edist["0"]/maxcombo)-dens)*100
	mean((1-edist["0"]/maxcombo)-dens)/sd((1-edist["0"]/maxcombo)-dens)

	apply(t(apply(permdb,1,"*",as.numeric(colnames(permdb)))),sum)/maxcombo

	round(cbind(md,sdd,t=md/sdd,p=0,t(cid)),digits=3)

	round(cbind(md=md/maxcombo,sdd=sdd/maxcombo,t=md/sdd,p=0),digits=4)[1:7,]
	sum(edist[-(1:7)])/maxcombo

	plot(as.table((md/sdd)[1:7]),type="l",ylim=range((md/sdd)[1:7]),lwd=3)
	abline(h=0,lty="dotted",lwd=3)

	round(cbind(md,sdd,t=md/sdd,p=0),digits=1)[1:7,]
	sum(edist[-(1:7)])

	mp<-apply(dif/maxcombo,2,mean)
	sdp<-apply(dif/maxcombo,2,sd)
	cip<-apply(dif/maxcombo,2,quantile,prob=c(.05,.95))

	round(cbind(mp,sdp,p=0,t(cip)),digits=5)[1:7,]
	sum(edist[-(1:7)])

	mvd<-mean((sum(((as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2)*edist)/maxcombo)-(apply(apply(permdb,1,"*",(as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2),2,sum)/maxcombo))
	sddv<-sd((sum(((as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2)*edist)/maxcombo)-(apply(apply(permdb,1,"*",(as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2),2,sum)/maxcombo))
	cidv<-quantile((sum(((as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2)*edist)/maxcombo)-(apply(apply(permdb,1,"*",(as.numeric(names(edist))-sum(lc2mel$ew)/maxcombo)^2),2,sum)/maxcombo),prob=c(.05,.95))
	round(c(mvd=mvd,sddv=sddv,tvd=mvd/sddv,p=0,cidv=cidv),digits=4)

	m<-apply(permdb,2,mean)
	sd<-apply(permdb,2,sd)
	ci<-apply(permdb,2,quantile,prob=c(.05,.95))
	cbind(m,sd,ci)

	zscores<-(edist[1:7]-apply(permdb,2,mean))/apply(permdb,2,sd)


	apply(matrix(edist,nrow=dim(permdb)[1],ncol=dim(permdb)[2],byrow=T)-permdb,2,mean)
	apply(matrix(edist,nrow=dim(permdb)[1],ncol=dim(permdb)[2],byrow=T)-permdb,2,sd)
}

extract.components<-function(
	graph
	,quantile=0
)
{
	components<-component.dist(graph)
	subgraphs<-list()
	j<-1
	for(i in order(components$csize[which(components$csize>=quantile(components$csize, c(quantile)))],decreasing=T)){
		subgraphs[[j]]<-graph%s%which(components$membership==i)
		j<-j+1
	}
	cat("\nComponent size distribution:")
	print(table(components$csize))
	return(subgraphs)
}


igraph2statnet<-function(
	net
)
{
	require(igraph)
	require(network)
	require(data.table)

	if(class(net)=='igraph'){
		el<-data.table(get.edgelist(net))
		dir<-!igraph::is.directed(net)
		if(dir){
			w<-el[,V1>V2]
			r<-el[w,V1]
			el[w,V1:=V2]
			el[w,V2:=r]
		}
		setkey(el,V1,V2)
		bf<-nrow(el)
		el<-unique(el)
		cat('\n',bf-nrow(el),'duplicate edges deleted\n')
		el<-as.matrix(el)
		net<-network::network(el,matrix.type="edgelist",directed=ifelse(dir,T,F))
	}

	if(class(net)=='network'){

	}
	net
}
