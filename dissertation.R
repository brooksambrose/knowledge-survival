cat('\014') # clear console
#############################
############ TOC ############
#############################

# SECTION 1
## 1.1 load database
## 1.2 extract bimodal edge list
## 1.3 convert to two unimodal edge lists
# SECTION 2

###############################
############ SETUP ############
###############################

rm(list=ls()) #clear memory
setwd('~/Dropbox/GitHub/knowledge-survival') # work in Git repo
source('dissertation_source.R') # load source
rm(list=ls()[!ls()%in%c(
'wok2db.f'
,'wok2db'
,'db2bel.f'
,'db2bel'
,'fuzzy.sets'
,'bel2mel.f'
,'bel2mel'
,'.ls.objects'
,'lsos'
)]) # functions we're using so far
tc<-1 # a counter for which table we are working on
fc<-1 # a counter for which figure we are working on
library(data.table)

###################################
############ SECTION 1 ############
###################################

############ 1.1 load database ############
if(F){wok2db<-wok2db.f(
	dir='in'
	,out='out'
	,sample.batches=F
	,sample.size=100
	,save=F
)
save(wok2db,file='wok2db.RData')
}

## 1.1.1 Export list of journals
library(stringdist)
library(igraph)
library(data.table)
load('wok2db.RData')
setkey(wok2db,field,id,val)
wok<-wok2db[list('SO'),.N,by=val]$val
#wok<-do.call(c,strsplit(wok,'-'))
rm(wok2db)
jstoro<-read.delim('/Users/bambrose/Dropbox/2013-2014/forLynne/WOKsamp.txt')
match<-data.frame(abbr=NA,jstor=jstoro$jstor[amatch(tolower(wok),tolower(jstoro$ssci))],wok=wok)
match<-match[order(match$wok,match$jstor),]
if(F) write.table(match,file='supplemental/1.1.1 JSTOR-WOK/JSTOR-samp-abbr.tab',quote=F,row.names=F,na='',sep='\t')
match<-read.table('supplemental/1.1.1 JSTOR-WOK/JSTOR-samp-abbr.tab',sep='\t',header=T)
match$abbr<-sub('^.+=','',match$abbr)

#these can be pasted into a browser and data requests can be made manually to dfr.jstor.org
cat(c('',paste('http://dfr.jstor.org/fsearch/submitrequest?cs=jcode%3A',match$abbr[match$abbr!=''],'&cc=ty%3Afla',sep=''),'\n',sep=''),sep='\n')

#these links will download csvs of annual record counts for each journal
lapply(
	c('fla','brv')
	,function(x) cat(
		paste('http://dfr.jstor.org/fsearch/csv?cs=jcode%3A',as.character(match$abbr[match$abbr!='']),'%5E1.0%7Cty%3A',x,'%5E1.0&fs=rtm1%3Ayrm1%3Atym1%3Asnm1&view=text&&csv=yr&fmt=csv',sep='')
		,sep='\n'))

"http://dfr.jstor.org/fsearch/csv?cs=jcode%3Aamerjarch%5E1.0%7Cty%3Abrv%5E1.0&fs=rtm1%3Ayrm1%3Atym1%3Asnm1&view=text&&csv=yr&fmt=csv"

"http://dfr.jstor.org/fsearch/csv?cs=jcode%3Aamerjarch%5E1.0%7Cty%3Abrv%5E1.0&fs=rtm1%3Ayrm1%3Atym1%3Asnm1&view=text&&csv=yr&fmt=csv"

#tabulate data sources across web of knowledge and jstor
match<-data.table(match)
setkey(match,wok)
setkey(wok2db,field,val)
wok.so<-wok2db['SO',.N,by=val]
setnames(wok.so,'val','wok')
wok.so[match]

rm(list=ls())
gc()
cat('\014')
t0<-proc.time()
# after downloading

#go to dissertation_source and look for jstor2db.
match<-read.table('supplemental/1.1.1 JSTOR-WOK/JSTOR-samp-abbr.tab',sep='\t',header=T)
missing<-match$jstor[!match$jstor%in%dfr.jour$journaltitle]
(missing<-as.character(missing[missing!='']))

#venn diagram for comparing jstor and wok samples
require(venneuler)
v <- venneuler(c(A=450, B=1800, "A&B"=230))
plot(v)

# just to 1920
ret[,y:=as.integer(substr(pubdate,1,4))]
setkey(ret,y)
ret00.41<-ret[list(1900:1941)]
rm(ret)
save(ret00.41,file='ret00.41.RData')
# topic model
require(stm)
vocab<-table(unlist(ret$bow))
vocab<-sort(unique(unlist(docs))) # stm expects as input a list of matrices where each element of the list is a document and where the first row of each matrix is the index of a word and the second row is the count of the word. This is a more memory efficient form since zeros are not stored.
stm.docs<-list()
for(i in names(docs)){
	t<-table(docs[[i]])
	stm.docs[[i]]<-rbind(vocab.index=which(vocab%in%names(t)),frequency=t)
}

pre2stm<-list()
pre2stm$model<-stm(documents=stm.docs,vocab=vocab,K=k,control=list(alpha=alpha))
pre2stm$top.word.phi.beta<-sapply(data.frame(pre2stm$model$beta$logbeta),function(x) sapply(x,function(y) ifelse(is.infinite(y),.Machine$double.eps,exp(y)))) # called beta by stm, epsilon closest thing to zero the machine can represent, necessary to prevent error
colnames(pre2stm$top.word.phi.beta)<-pre2stm$model$vocab
pre2stm$doc.top.theta<-pre2stm$model$theta
rownames(pre2stm$doc.top.theta)<-names(docs)
pre2stm$doc.length<-sapply(docs,length)
pre2stm$vocab<-pre2stm$model$vocab
tn<-table(unlist(docs))
pre2stm$term.frequency<-as.integer(tn)
names(pre2stm$term.frequency)<-names(tn)

save(pre2stm,file=sub(paste(rep(.Platform$file.sep,2),collapse=""),.Platform$file.sep,paste(out.dir,paste("stm-model-k",k,"-alpha",round(alpha,3),".RData",sep=""),sep=.Platform$file.sep)))


if(F){
	jstor<-as.character(read.csv('supplemental/1.1.1 JSTOR-WOK/JSTOR-all-journals.csv',skip=2)$Journal)
	jstor<-do.call(c,strsplit(jstor,' / '))
		jstor[jstor$ssci!='',]
	jnet<-stringdistmatrix(tolower(c(wok,jstor)),tolower(c(wok,jstor)),useNames=T,method='jw',p=.1)
	jnet[1:length(wok),1:length(wok)]<-1
	match<-data.frame(wok,jstor=apply(jnet[1:length(wok),],1,function(x) jstor[which.min(x)-length(wok)]))
	rownames(match)<-NUL
	# jnet<-jnet/max(jnet)
	# jnet<-graph.adjacency(jnet,mode='upper',weighted=T,diag=F)
	# jcom<-spinglass.community(jnet,spins=300)
	# (jcom<-split(names(membership(jcom)),membership(jcom)))
}


############ 1.2 resolve identity uncertainty for citations ############

## 1.2.1 Greedy string comparison (already performed via hoffman)

## 1.2.2 pairwise comparisons within greedy sets

## 1.2.3 community detection to cluster pairwise comparisons


############ 1.3 compile bimodal edge list ############

## 1.3.1 pre-process edge list

if(F){
load('wok2db.RData')
setkey(wok2db,id)
system.time(db2bel<-db2bel.f(
	wok2db=wok2db
	,out='out'
	,saved_recode=fuzzy.sets
))
}

## 1.2.4 recode original, statistics on the difference it made
if(F){
	load('fuzzy-sets.RData')
	lsos()
	setkey(wok2db,field)
	years<-as.data.frame(wok2db['PY',list(val)][,.N,keyby=val])
	years$perc<-round(prop.table(years$N)*100,2)
	print(years)

	setkey(wok2db,field)
	db2cr<-wok2db['CR',list(id,val)]
}


############ 1.4 convert to two unimodal edge lists ############
load('db2bel.RData')
lsos()
db2bel<-db2bel[which(!(zpend|zdup)),list(ut,zcr)] # use only recoded non-pendant
lsos()
db2bel[,`:=`(ut=droplevels(ut),zcr=droplevels(zcr))]
lsos()

system.time(bel2mel<-bel2mel.f(
	db2bel=db2bel
	,out='out'
))



#export to CFinder

for(i in names(bel2mel)){
setwd('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out')
long<-factor(c(as.character(bel2mel[[i]][[1]]),as.character(bel2mel[[i]][[2]])))
longlev<-levels(long)
writeLines(longlev,paste(i,'_1941_CFinder_levels.txt',sep=''))
write.table(cbind(matrix(as.numeric(long),ncol=2),bel2mel[[i]][['ew']]),file=paste(i,'_1941_CFinder.txt',sep=''),quote=F,na='',row.names=F,col.names=F)
}

#reading CFinder output
#import k-clique community relatedness graph
rm(list=ls())
db<-'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/CFinder1941-exact'
db<-list.files(db,recursive=T,full.names=T)
goc<-grep('graph_of_communities',db)
com<-grep('/communities$',db)
lk<-grep('links',db)
mem<-grep('membership_distribution',db)



#just trying to get networks within communities here
l<-list() #a list with each undifferentiated link record in string form
for(i in lk) l[[length(l)+1]]<-readLines(db[i])
for(i in 1:length(l)) l[[i]]<-l[[i]][8:length(l[[i]])]
m<-list() #a list with cutpoints for each link record
for(i in 1:length(l)) {m[[i]]<-grep(':',l[[i]]);m[[i]]<-c(m[[i]],length(l[[i]]))}


f<-list()
for(i in com) f[[length(f)+1]]<-readLines(db[i])
for(i in 1:length(f)) f[[i]]<-f[[i]][8:length(f[[i]])]
f<-lapply(f,strsplit,split=' ')
for(i in 1:length(f)) for(j in 1:length(f[[i]])) f[[i]][[j]]<-as.numeric(f[[i]][[j]][-1])
sdis<-list()
for(i in 1:length(f)) sdis[[i]]<-sapply(f[[i]],length)

g<-list()
for(i in goc) {
	g[[length(g)+1]]<-readLines(db[i])
}

q<-list()
for(i in mem) {
	q[[length(q)+1]]<-readLines(db[i])
}


k<-list()
for(i in 1:length(g)) {
	k[[i]]<-as.integer(sub('^.+=','',g[[i]][1]))
}
o<-order(unlist(k))
g<-g[o]
l<-l[o]
m<-m[o]
q<-q[o]
f<-f[o]
k<-sort(unlist(k))
names(g)<-paste('k',k,sep='')
names(l)<-names(g)
names(q)<-names(g)
names(m)<-names(g)
names(f)<-names(g)
p<-list() #a list of k levels with a list of each k-clique edgelist within
for(i in names(m)) {
	nam<-paste('l',0:(length(m[[i]])-1),sep='')
	for(j in 1:(length(nam)-1)) p[[i]][[nam[j]]]<-lapply(strsplit(l[[i]][(m[[i]][j]+1):(m[[i]][j+1]-1)],' '),as.numeric)
	}

names(p)<-names(g)
{g<-g[!duplicated(g,fromLast=T)];attr(g,'source')<-'graph_of_communities'}
{l<-l[!duplicated(l,fromLast=T)];attr(l,'source')<-'communities_links'}
{m<-m[!duplicated(m,fromLast=T)];attr(m,'source')<-'communities_links'}
{q<-q[!duplicated(q,fromLast=T)];attr(q,'source')<-'membership_distribution'}
{f<-f[!duplicated(f,fromLast=T)];attr(f,'source')<-'communities'}
{p<-p[!duplicated(p,fromLast=T)];attr(p,'source')<-'communities_links'}
for(i in names(p)) for(j in names(p[[i]])) p[[i]][[j]]<-do.call(rbind,p[[i]][[j]])
for(i in names(f)) names(f[[i]])<-paste('l',0:(length(f[[i]])-1),sep='')
gfpq<-list(g=g,f=f,p=p,q=q)
save(gfpq,file='gfpq.RData')


#### Hierarchical clustering of communities
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/gfpq.RData')
f<-gfpq$f
rm(gfpq)
f<-do.call(c,f)
f<-lapply(f,sort)
f<-f[!duplicated(f,fromLast=T)]
cn<-list()
for(i in names(f)) cn[[i]]<-sapply(lapply(f,'%in%',f[[i]]),all) # for every community which communities does it completely contain?
cn<-do.call(data.frame,cn)
diag(cn)<-F # each community contains itself, so ignore
cn<-sapply(cn,which) #becomes a list of containers and contained
nonest<-!sapply(cn,length)
cn<-cn[!nonest] # get rid of communities who contain nooone else
for(i in names(cn)) cn[[i]]<-names(f)[cn[[i]]] ### provides a list of containers and all of the communities that they completely contain
cf<-f
for(i in names(cn)) cf[[i]]<-setdiff(f[[i]],unique(unlist(f[cn[[i]]])))
for(i in names(cn)) cn[[i]]<-setdiff(cn[[i]],unique(unlist(cn[cn[[i]]]))) ### so communities are only contained by the highest k
#odd that a k3 community contains another k3 comunity

xml<-lapply(cn,as.list)
for(i in rev(names(xml))) {
	names(xml[[i]])<-cn[[i]]
	for(j in names(xml[[i]])) xml[[i]][[j]]<-''
}

ix<-rep(names(cn),sapply(cn,length))
jx<-sub('.+\\.([0-9])','\\1',make.unique(ix))
jx[grep('^k',jx)]<-'0'
jx<-type.convert(jx)+1
hx<-rep(1:length(cn),sapply(cn,length))
kx<-unlist(cn,use.names=F)
for(i in rev(names(xml))) if(i%in%kx) {
	j<-which(kx%in%i)
	xml[[hx[j]]][[jx[j]]]<-xml[[i]]
	xml<-xml[!names(xml)%in%i]
}
nm<-names(unlist(xml))
if(F){
nk<-lapply(strsplit(gsub('[^0-9]+','-',nm),'-'),as.integer)
srt<-matrix(Inf,nrow=length(nm),ncol=max(sapply(nk,length)))
for(i in 1:length(nk)) srt[i,1:length(nk[[i]])]<-nk[[i]]
srt<-srt[,seq(2,dim(srt)[2],2)]
so<-as.list(data.frame(sapply(nk,length),srt))
so<-so[c(2,1,3:length(so))]
so<-do.call(order,so)
nm<-nm[so]
nk<-nk[so]
ha<-do.call(rbind,strsplit(sub('.+\\.(k[^k]+)\\.(k[^k]+)$','\\1-\\2',nm),'-'))
### not the right way to do it, should just alphebatize so that each level is taken care of
}

nm<-sort(nm)
nm<-gsub('\\.k','-k',nm)
nt<-strsplit(nm,'-')
### if a community is size (tt-1) or smaller, add it into the category it was derived from. This creates a condition in which small spans are not limited simply by the size of the category. Leaving only communities of size (tt-1) ameliorates this potential, and this is fair because these straglers are actually members of that community.
temp<-rep(names(cn),sapply(cn,length))
names(temp)<-unlist(cn)
ak<-round(sapply(lapply(lapply(nt,sub,pattern='k([0-9]+).+','\\1',replacement='\\1'),as.integer),mean),1)
cf<-scf
sf<-sapply(cf,length)
tt<-2 # tiny threshold, will retain this number and higher. must be higher than 1 or ut=1 will perfectly predict failure
tiny<-names(cf[sf>0&sf<tt])
run<-1

if(T){ while(length(tiny)){
	cat('\n\n########### Run',run,'###########')
	cat('\n#### Length:',length(tiny))
	run<-run+1
	for(i in tiny) {
	cat('\n\n',i,sep='')
	pkw<-sapply(nt,grep,pattern=paste(i,'$',sep=''))
	pkx<-which(!!sapply(pkw,length))
	pkw<-unlist(pkw)
	cat('pkw:',pkw,'pkx:',pkx,'\n')
	cat('',nm[grep(i,nm)],sep='\n')
	inh<-mapply(FUN=function(x,w) x[w+ifelse(unique(pkw)==1,1,-1)],x=nt[pkx],w=pkw) #inheritor
	ak1<-table(ak[pkx])
	if(ak1[[as.character(max(ak[pkx]))]]==1) {inh<-inh[which.max(ak[pkx])];cf[[inh]]<-union(cf[[inh]],cf[[i]])} #if there are multiple inheritors, add to least cohesive. if tied, add nowhere and delete
	cf[[i]]<-numeric(0)
	}
sf<-sapply(cf,length)
tiny<-setdiff(names(cf[sf>0&sf<tt]),ifelse(sf[[inh]]>0&sf[[inh]]<tt,inh,''))
}}

#eliminate empty rings
sfut<-lapply(lapply(nt,FUN=function(x) rev(sf[x])),diff)
sfut<-data.frame(
	kl=unlist(lapply(lapply(nt,FUN=rev),FUN=function(x) x[-1]))
	,diff=unlist(sfut)
)
sfut<-aggregate(sfut$diff,by=sfut['kl'],sum)
sfut<-as.character(sfut$kl[sfut$x==0]) # categories that are empty and followed by an empty set
e<-lapply(nt,FUN=function(x) x[!sf[x]]) # categories that are empty and have only one child
e<-table(unlist(e))
sfut<-union(sfut,names(e[e==1]))
nt<-lapply(nt,setdiff,y=sfut)

nu<-unique(unlist(nt))
un<-table(unlist(nt))
un<-un[nu]
nut<-nu
tl<-sapply(nt,tail,1)
for(i in 1:length(nut)) if(any(tl==nut[i])) nut[i]<-nt[tl==nut[i]]

nutl<-list() ## list of when categories open and close
for(i in nu) nutl[[i]]<-range(which(sapply(lapply(nut,grepl,pattern=paste(i,'$',sep='')),any)))
nutl<-do.call(rbind,nutl)
nutl<-cbind(nutl,!apply(nutl,1,diff))
cl<-table(nutl[!nutl[,3],2])

binc<-as.integer(sub('k([0-9]+).+','\\1',nu))
binc<-c(nmax=-max(binc),med=round(quantile(binc,.25)),zero=0,avg=mean(binc))

nh<-lapply(nu,FUN=function(x){
	y<-nt[[grep(x,nm)[1]]]
	w<-which(y==x)
	w<-try(y[0:(w-1)],silent=T)
	if(class(w)=='try-error'|!length(w)) w<-NA_character_
	w
	})
names(nh)<-nu

hkl<-list(memb=cf[!!sapply(cf,length)],nest=nh,thresh=tt)

#### Discover components in bel, isolates in mel (dropped by bel2mel)
levs<-readLines('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/tpcrel_1941_CFinder_levels.txt')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
utcrp<-aggregate(!db2bel$pend,by=db2bel$bel['ut'],mean) #proportion of cr for each ut that are pendant
utcrp<-data.frame(utcrp[1],round(utcrp[2],3))
utcrr<-aggregate(db2bel$bel$cr[db2bel$pend],by=list(db2bel$bel[db2bel$pend,'ut']),length) #number of cr remaining after pendants dropped
colnames(utcrr)<-c('ut','crr')
colnames(utcrp)<-c('ut','crp')

iso<-merge(db2bel$bel[db2bel$pend,],utcrr,all.y=T)
iso<-aggregate(iso$crr,iso['cr'],FUN=function(x) all(x==1))
iso<-iso[iso$x,'cr'] #the CR's lost because isolate in unimodal projection; would be components all UT attached to one CR in bimodal
isow<-table(db2bel$bel$cr[db2bel$bel$cr%in%iso])
isot<-table(isow)
hkl$iso<-isow

#add isolate communities, but only if (tt-1) or more UT or else ut==1 perfectly predicts failure
hkl$isoin<-names(isow[isow>=hkl$thresh])
a<-length(hkl$isoin)
b<-length(levs)
hkl$mlev<-c(levs,hkl$isoin)
c<-(b+1):(b+a)
c<-as.list(c)
names(c)<-paste('i',1:a,sep='')
hkl$memb<-c(hkl$memb,c)


#export to xml for visualization
if(T){
library(XML)
treeviz<-suppressWarnings(xmlTree('k0',attrs=c(name='unmarked',size='0',k='0')))
if('nunstr'%in%ls()) treeviz$addNode('k1',attrs=c(name='unclustered',size=nunstr,k=1),close=T)
for(i in 1:length(nu)) {
	treeviz$addNode(nu[i]
		,attrs=c(
			name=nu[i]
			,size=ifelse(nutl[[i,3]],sf[[nu[i]]],0)
			,k=ifelse(nutl[[i,3]],as.integer(sub('k([0-9]+).+','\\1',nu[i])),binc['zero']))
		,close=ifelse(nutl[i,3],T,F))
	if(!nutl[i,3]&sf[[nu[i]]]>0) treeviz$addNode(nu[i],attrs=c(name=nu[i],size=sf[[nu[i]]],k=as.integer(sub('k([0-9]+).+','\\1',nu[i]))),close=T)
	if(i%in%names(cl)) replicate(cl[as.character(i)],treeviz$closeTag())
}
for(i in setdiff(names(f)[!names(f)%in%unlist(nt)],sfut)) treeviz$addNode(i,attrs=c(name=i,size=sf[[i]],k=as.integer(sub('k([0-9]+).+','\\1',i))),close=T)
for(i in grep('^i',names(hkl$memb),value=T)) treeviz$addNode(i,attrs=c(name=i,size=1,k='i'),close=T)
treeviz$closeTag()
cat(saveXML(treeviz),'\n')
saveXML(treeviz, file='/Users/bambrose/Downloads/test.xml')
}

save(hkl,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')

if(F){
#inspect subsets here
source('/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')
library(statnet)

#view size of ut and cr dimension
utgr<-data.frame(ut=aggregate(d$ut,d['idn'],sum)$x,k=aggregate(d$k,d['idn'],unique)$x)
rownames(utgr)<-unique(d$id)

tut<-50 # thresh for viewing size of UT
plot(utgr[utgr$ut<=tut,],cex=0)
text(utgr[utgr$ut<=tut,],labels=rownames(utgr[utgr$ut<=tut,]),cex=.5)
utgr<-utgr[order(utgr$ut,utgr$k,decreasing=T),]

t1<-proc.time()
net<-lapply(rownames(utgr),FUN=function(x) {cat(x,'\n',sep='\n');subnet(db2bel=db2bel,set=list(cr=hkl$mlev[hkl$memb[[x]]],ut=NULL,nut=NULL,ncr=NULL),source='/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R')})
print(proc.time()-t1)
names(net)<-rownames(utgr)
net[!sapply(net,length)]<-
net<-do.call(rbind,net)
save(net,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/utknets_41.RData')
utnet<-na.omit(net[,'tputel'])
sapply(net,gden)

netn.f<-function(net){
	require(statnet)
	netn<-rep('',network.size(net))
	netn[which.max(degree(net))]<-network.vertex.names(net)[which.max(degree(net))]
	netn
}
pdf('utknets_41.pdf',h=(8.5-2)/2,w=(8.5-2)/2)

mapply(n=,m=,FUN=function(n,m) plot(n,main=m,displaylabels=T,label=netn,label.cex=.5,label.bg='lightgray', boxed.labels=F,label.lwd=0,edge.lwd=.1,vertex.cex=.5,edge.col=gray(.5,.5),vertex.lwd=.1))
dev.off()


#ergmm
library(igraph)
library(statnet)
inet<-graph.adjacency(net[[2]][,],mode='undirected',weighted=NULL)
ic<-optimal.community(inet)
icm<-ic$membership
ergmm.prior(Z.pK=icm)
ich<-multilevel.community(inet)
comr<-apply(ich$memberships,1,max)
ichmax<-ich$memberships[which.max(comr),]
ichmin<-ich$memberships[which.min(comr),]
#par(mfrow=c(1,1),mar=rep(.5,4))
plot(create.communities(ichmin),inet,vertex.label=NA,vertex.size=5,main='Minimal Clusters')
plot(ic,inet,vertex.label=NA,vertex.size=5,main='Optimal Clusters')
plot(create.communities(ichmax),inet,vertex.label=NA,vertex.size=5,main='Maximal Clusters')

G<-do.call(seq,as.list(comr))
names(G)<-paste('d2',G,sep='G')
G<-rev(G)

fitlsm<-function(G,Gprior,net) {
	require(statnet)
	cat('\n\nG = ',G,'\n',sep='')
	t1<-proc.time()
		ergmm(net~euclidean(d=2,G=G)+rsociality()
		,prior=ergmm.prior(Z.K=Gprior)
		,verbose=2
		,control=control.ergmm(burnin=100000,interval=20,sample.size=5000))
	t2<-proc.time()
	print(t2-t1)
	}
# Gprior = [1] 3 5 8 2 4 1 1 2 3 4 2 5 4 5 5 6 2 7 2 6 6 6 7 5 8 2 2 5 6 8 6 4
f<-list()
for(i in rev(names(G))) f[[i]]<-try(fitlsm(G=G[[i]],Gprior=ichmax,net=net[[2]]))

plots<-list(f1=plot(fits[[1]],pie=T,suppress.axes=T))
plots<-c(plots,lapply(fits[-1],FUN=function(x) plot(x,pie=T,suppress.axes=T,Z.ref=plots$f1)))
}

#test model failure
{graphics.off();t<-network(pn[sample(1:nrow(pn),round(.75*nrow(pn))),1:2],directed=F);plot(t);tfit<-ergmm(t~euclidean(d=2,G=4)+rsociality(),verbose=2,control=control.ergmm(burnin=30000));plot(tfit,pie=T);mcmc.diagnostics(tfit)}
}
#### table of all cr's with data attached

######## Match UT and PY DT to each community membership
#rm(list=ls())


save(hkl,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')

####################### Exprot duration database to Stata?
#rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
tcr<-length(unique(db2bel$bel$cr))
cat(round(length(hkl$mlev)/tcr*100,2),'% of references contribute to community structure.',sep='')
dis<-c(
	sum(!db2bel$pend)
	,sum(hkl$iso<hkl$thresh)
	,length(unique(unlist(hkl$memb[grep('^k',names(hkl$memb))])))
	,length(unique(unlist(hkl$memb[grep('^i',names(hkl$memb))])))
)
names(dis)<-c('pendants',paste('iso.th',hkl$thresh,sep='<'),'structural')
dis<-cbind(dis,'%'=round(dis/sum(dis)*100,1))
write.table(dis,sep='\t',file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/cr_selection.tab')
print(dis)

cr<-unique(db2bel$bel$cr[db2bel$pend])
db2bel<-db2bel$bel[db2bel$pend,]
cr<-split(db2bel$cr,f=db2bel$ut)
(nunstr<-length(unique(db2bel$cr))-length(hkl$mlev)) #
rm(db2bel)
for(i in 1:length(cr)) cr[[i]]<-which(hkl$mlev%in%cr[[i]]) #convert CRs to CFinder codes

load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')
rownames(wok2db_41)<-NULL
py<-droplevels(wok2db_41[wok2db_41$field=='PY',-2])
colnames(py)<-c('ut','py')
dt<-droplevels(wok2db_41[wok2db_41$field=='DT',-2])
colnames(dt)<-c('ut','dt')
au<-droplevels(wok2db_41[wok2db_41$field=='AU',-2])
colnames(au)<-c('ut','au')
au<-cbind(ut=as.character(sort(unique(au$ut))),au=split(as.character(au$au),f=au$ut,drop=T))
af<-droplevels(wok2db_41[wok2db_41$field=='AF',-2])
colnames(af)<-c('ut','af')
af<-cbind(ut=as.character(sort(unique(af$ut))),af=split(as.character(af$af),f=af$ut,drop=T))
pu<-droplevels(wok2db_41[wok2db_41$field=='PU',-2])
colnames(pu)<-c('ut','pu')
so<-droplevels(wok2db_41[wok2db_41$field=='SO',-2])
colnames(so)<-c('ut','so')
sc<-droplevels(wok2db_41[wok2db_41$field=='SC',-2])
colnames(sc)<-c('ut','sc')
wc<-droplevels(wok2db_41[wok2db_41$field=='WC',-2])
colnames(wc)<-c('ut','wc')
wc<-cbind(ut=as.character(sort(unique(wc$ut))),wc=split(as.character(wc$wc),f=wc$ut,drop=T))
nr<-droplevels(wok2db_41[wok2db_41$field=='NR',-2])
colnames(nr)<-c('ut','nr')
tc<-droplevels(wok2db_41[wok2db_41$field=='TC',-2])
colnames(tc)<-c('ut','tc')
pg<-droplevels(wok2db_41[wok2db_41$field=='PG',-2])
colnames(pg)<-c('ut','pg')

db<-merge(py,dt,all=T)
db<-merge(db,au,all=T)
db<-merge(db,af,all=T)
db<-merge(db,pu,all=T)
db<-merge(db,so,all=T)
db<-merge(db,sc,all=T)
db<-merge(db,wc,all=T)
db<-merge(db,nr,all=T)
db<-merge(db,tc,all=T)
db<-merge(db,pg,all=T)

rownames(db)<-db$ut
rm(wok2db_41,py,dt,au,af,pu,so,sc,wc,nr,tc,pg)

hkl[['dbfull']]<-droplevels(db)

dbo<-!rownames(db)%in%names(cr) # just the ut's that have references
dbout<-db[dbo,]
db<-db[!dbo,]
db<-db[order(rownames(db)),]
cr<-cr[order(names(cr))]
db<-data.frame(db[rep(names(cr),sapply(cr,length)),],cr=unlist(cr))

hkl[['db']]<-droplevels(db)
rownames(hkl$db)<-NULL

hkl$dbfull<-merge(merge(hkl$dbfull,utcrp,all=T),utcrr,all=T)
hkl$db<-merge(merge(hkl$db,utcrp,all=T),utcrr,all.x=T)

save(hkl,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')



#split into spans
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')
s<-lapply(hkl$memb,FUN=function(x) droplevels(merge(matrix(x,ncol=1,dimnames=list(NULL,'cr')),y=hkl$db,by='cr')))
for(i in names(s)) s[[i]]<-data.frame(id=rep(i,nrow(s[[i]])),s[[i]])
s<-lapply(s,FUN=function(x) {
	x<-data.frame(lapply(x,FUN=function(y) type.convert(as.character(y))))
	x<-x[order(x$py,x$ut,x$cr),]
	rownames(x)<-NULL
	x<-data.frame(x,span=1+cumsum(c(1,diff(x$py))>1))
	x<-split(x,f=x$span,drop=T)
	names(x)<-paste('s',1:length(x),sep='')
	x
})

#custom aggregate data to duration format
s<-do.call(c,s)
d<-lapply(s,FUN=function(x) {
	n<-aggregate(x$ut,x[c('py','ut')],length)
	n<-n[order(n$py),]$x
	n<-rep(n,n)
	data.frame(
		id=aggregate(x$id,x['py'],unique)$x
		,py=aggregate(x$py,x['py'],unique)$x
		,fail=aggregate(x$py,x['py'],function(y) as.integer(unique(y==max(x$py))))$x
		,span=aggregate(x$span,x['py'],function(y) unique(y))$x
		,start=rep(min(x$py),length(unique(x$py)))
		,end=rep(max(x$py),length(unique(x$py)))
		,ut=aggregate(x$ut,x['py'],function(y) length(unique(y)))$x
		,au=aggregate(x$au,x['py'],function(y) length(unique(unlist(y))))$x
		,cr=aggregate(x$cr,x['py'],function(y) length(unique(y)))$x
		,crr=aggregate(x$crr/n,x['py'],sum)$x
		,nr=aggregate(x$nr/n,x['py'],sum)$x
		,tc=aggregate(x$tc/n,x['py'],sum)$x
		,pg=aggregate(x$pg/n,x['py'],sum)$x
		,dta=aggregate(as.integer(grepl('Article',x$dt))/n,x['py'],sum)$x
		,dtb=aggregate(as.integer(grepl('Book Review',x$dt))/n,x['py'],sum)$x
		,dtr=aggregate(as.integer(grepl('^Review$',x$dt))/n,x['py'],sum)$x
		,puwb=aggregate(as.integer(grepl('WILEY',x$pu))/n,x['py'],sum)$x
		,puaea=aggregate(as.integer(grepl('AMER ECONOMIC ASSOC',x$pu))/n,x['py'],sum)$x
		,puucp=aggregate(as.integer(grepl('UNIV CHICAGO PRESS',x$pu))/n,x['py'],sum)$x
		,puasa=aggregate(as.integer(grepl('AMER SOCIOLOGICAL ASSOC',x$pu))/n,x['py'],sum)$x
		,pusag=aggregate(as.integer(grepl('SAGE PUBLICATIONS INC',x$pu))/n,x['py'],sum)$x
		,puaaea=aggregate(as.integer(grepl('AGRICULTURAL',x$pu))/n,x['py'],sum)$x
		,puunc=aggregate(as.integer(grepl('CAROLINA',x$pu))/n,x['py'],sum)$x
		,pumit=aggregate(as.integer(grepl('M ?I ?T PRESS',x$pu))/n,x['py'],sum)$x
		,puox=aggregate(as.integer(grepl('OXFORD',x$pu))/n,x['py'],sum)$x
		,puspr=aggregate(as.integer(grepl('SPRINGER',x$pu))/n,x['py'],sum)$x
		,soaa=aggregate(as.integer(grepl('AMERICAN ANTHROPOLOGIST',x$so))/n,x['py'],sum)$x
		,soaer=aggregate(as.integer(grepl('AMERICAN ECONOMIC REVIEW',x$so))/n,x['py'],sum)$x
		,soajs=aggregate(as.integer(grepl('AMERICAN JOURNAL OF SOCIOLOGY',x$so))/n,x['py'],sum)$x
		,soasr=aggregate(as.integer(grepl('AMERICAN SOCIOLOGICAL REVIEW',x$so))/n,x['py'],sum)$x
		,soehr=aggregate(as.integer(grepl('ECONOMIC HISTORY REVIEW',x$so))/n,x['py'],sum)$x
		,soej=aggregate(as.integer(grepl('^ECONOMIC JOURNAL$',x$so))/n,x['py'],sum)$x
		,soe=aggregate(as.integer(grepl('ECONOMICA',x$so))/n,x['py'],sum)$x
		,sojes=aggregate(as.integer(grepl('JOURNAL OF EDUCATIONAL SOCIOLOGY',x$so))/n,x['py'],sum)$x
		,sojfe=aggregate(as.integer(grepl('JOURNAL OF FARM ECONOMICS',x$so))/n,x['py'],sum)$x
		,sojpe=aggregate(as.integer(grepl('JOURNAL OF POLITICAL ECONOMY',x$so))/n,x['py'],sum)$x
		,sojsf=aggregate(as.integer(grepl('JOURNAL OF SOCIAL FORCES',x$so))/n,x['py'],sum)$x
		,soqje=aggregate(as.integer(grepl('QUARTERLY JOURNAL OF ECONOMICS',x$so))/n,x['py'],sum)$x
		,sores=aggregate(as.integer(grepl('REVIEW OF ECONOMIC STUDIES',x$so))/n,x['py'],sum)$x
		,sosf=aggregate(as.integer(grepl('SOCIAL FORCES',x$so))/n,x['py'],sum)$x
		,sos=aggregate(as.integer(grepl('SOCIOMETRY',x$so))/n,x['py'],sum)$x
		,sozfn=aggregate(as.integer(grepl('ZEITSCHRIFT FUR NATIONALOKONOMIE',x$so))/n,x['py'],sum)$x
		,sca=aggregate(as.integer(grepl('Anthropology',x$sc))/n,x['py'],sum)$x
		,scp=aggregate(as.integer(grepl('^Psychology$',x$sc))/n,x['py'],sum)$x
		,scs=aggregate(as.integer(grepl('Sociology',x$sc))/n,x['py'],sum)$x
		,scbe=aggregate(as.integer(grepl('^Business . Economics',x$sc))/n,x['py'],sum)$x
		,sce=aggregate(as.integer(grepl('Education',x$sc))/n,x['py'],sum)$x
		,scag=aggregate(as.integer(grepl('Agriculture',x$sc))/n,x['py'],sum)$x
		,scm=aggregate(as.integer(grepl('Math',x$sc))/n,x['py'],sum)$x
		)
})

d<-do.call(rbind,d)
rownames(d)<-NULL
d<-data.frame(d[,1:4],gap=NA,d[,5:ncol(d)])
w<-c(0,diff(d$span))&!c(0,diff(d$id))
for(i in which(w)) d$gap[d$id==d$id[i]&d$span==d$span[i]]<-(c(1,diff(d$py))-1)[i]
d<-d[d$py!=max(d$py),]
d<-data.frame(idn=as.numeric(d$id),d)
w<-which(colnames(d)=='end')
w2<-which(colnames(d)=='dtr')
depth<-sapply(hkl$nest,length)
d<-data.frame(
	sid=paste(d$id,d$span,sep='.s')
	,d[,1:w]
	,iso=ifelse(grepl('^i',d$id),1,0)
	,k=as.integer(sub('k([0-9]+).+','\\1',d$id))
	,depth=0
	,d[,(w+1):w2]
	,dtac=unlist(aggregate(d$dta,d['id'],cumsum)$x) # cumulative sum of reviews
	,dtbc=unlist(aggregate(d$dtb,d['id'],cumsum)$x) # cumulative sum of reviews
	,dtrc=unlist(aggregate(d$dtr,d['id'],cumsum)$x) # cumulative sum of reviews
	,d[,(w2+1):ncol(d)])
d$depth[d$id%in%names(depth)]<-depth[as.character(d$id[d$id%in%names(depth)])]

n<-table(unlist(hkl$nest))
nn<-na.omit(unique(unlist(hkl$nest)))
nn<-nn[order(as.integer(sub('^k([0-9]+).+','',nn),as.integer(sub('^.+l([0-9]+)','',nn))))]
names(nn)<-nn
hkl$nest<-hkl$nest[!is.na(hkl$nest)]
d<-data.frame(d
	,lapply(as.list(nn)
	,FUN=function(x) sapply(hkl$nest[d$id]
		,FUN=function(y) ifelse(is.null(y),yes=0,no=as.integer(any(x==y)))
		)
)
)
save(d,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/duration.RData')
write.table(d
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/stata/duration.tab'
	,quote=F
	,sep='\t'
	,na='.'
	,row.names=F
)

############ DESCRIPTIVES ############
##### Illustrate network projections that are small versions of my distribution

if(F){
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
source('/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/edge_sample.RData')
library(statnet)

narts<-list()
for(i in 1:10){
	es<-list(bel=db2bel$bel[db2bel$pend,][sample(1:sum(db2bel$pend),2000),])
	es$bel<-es$bel[order(es$bel$ut,es$bel$cr),]
	tcr<-table(es$bel$cr)
	es$pend<-es$bel$cr%in%names(tcr)[tcr>1]
	#save(es,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/edge_sample.RData')
	utd<-table(es$bel$ut)
	utq<-quantile(utd,prob=seq(.01,1,.01))
	utd<-table(utd)

	crd<-table(es$bel$cr)
	crq<-quantile(crd,prob=seq(.01,1,.01))
	crd<-table(crd)

	cat('\n\ni =',i)
	cat('\n')
	print(cbind(utd,round(prop.table(utd),3)))
	cat('\n')
	print(cbind(crd,round(prop.table(crd),3)))
	cat('\n')
	print(tail(rbind(cbind(utq=utq,crq=crq),tot=apply(cbind(utq=utq,crq=crq),2,sum)),11))
	fxygt<-function(x,y,gt=T) if(sum(x)>sum(y)) x else y
	fxylt<-function(x,y,gt=T) if(sum(x)<sum(y)) x else y

	longer<-round(fxygt(utq,crq))
	shorter<-round(fxylt(utq,crq))
	trunc<-longer[(1+which(cumsum(longer)==abs(sum(longer)-sum(shorter)))):length(longer)]
	for(j in 1:10){
	sim<-data.frame(
		m1=sample(rep(
			paste('m1.',sapply(mapply(rep,0,times=3-(nchar(1:length(trunc)))),paste,collapse=''),1:length(trunc),sep='')
			,trunc))
		,m2=sample(rep(
			paste('m2.',sapply(mapply(rep,0,times=3-(nchar(1:length(shorter)))),paste,collapse=''),1:length(shorter),sep='')
			,shorter))
	,stringsAsFactors=F)
	sim<-sim[order(sim$m1,sim$m2),]
	sim<-sim[!duplicated(sim),]
	rownames(sim)<-NULL
narts[[paste(i,j,sep='_')]]<-plot.mode.projection(sim
	,out=paste('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/',i,'_',j,'_',sep='')
	,trim.pendants='m2',ecol='black',elt=0,vlw=.5,m1vsides=20)
}
}
save(narts,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/narts.RData')

load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/narts.RData')
source('/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R')
el<-as.edgelist(narts$'6_7'$pbmnet)
el<-t(apply(el,1,sort))
el<-el[!duplicated(el),]
plot.mode.projection(el
	,out='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/painted'
	,trim.pendants='m2',ecol='darkgray',elt='blank',vcx=2,elw=5,vlw=.5,m1vsides=20,m1col='white',m2col='white'
	,pnt.v=list(
		darkgray='m1.085'
		,black='m2.100'
	)
	,pnt.e=list(black=cbind('m2.100','m2.100'))
	,pnt.vlty=list(
		solid='m1.085'
		,solid='m2.100'
	)
	,layout=list(
		repulse.rad=function(x) x^2*log(x)/2
		,niter=function(x) 500
	)
)

load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/narts.RData')
library(igraph)
ig<-graph.adjacency(narts$'6_7'$bmnet[,])
cig<-optimal.community(ig)
pdf('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/igraph_opt.pdf')
plot(cig,ig)
dev.off()

}

#####Some summary stats for Scrivener
if(F){
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')

wok2db_41<-wok2db_41[wok2db_41$field%in%c('PY','PU','SO','DT','PG','AU','AF'),]
colnames(wok2db_41)<-c('ut','field','b')
rownames(wok2db_41)<-NULL
q<-droplevels(wok2db_41[wok2db_41$field=='AU',c('ut','b')])
colnames(q)[2]<-'au'
r<-droplevels(wok2db_41[wok2db_41$field=='AF',c('ut','b')])
colnames(r)[2]<-'af'
s<-droplevels(wok2db_41[wok2db_41$field=='PG',c('ut','b')])
colnames(s)[2]<-'pg'
t<-droplevels(wok2db_41[wok2db_41$field=='DT',c('ut','b')])
colnames(t)[2]<-'dt'
u<-droplevels(wok2db_41[wok2db_41$field=='PY',c('ut','b')])
colnames(u)[2]<-'py'
v<-droplevels(wok2db_41[wok2db_41$field=='PU',c('ut','b')])
colnames(v)[2]<-'pu'
w<-droplevels(wok2db_41[wok2db_41$field=='SO',c('ut','b')])
colnames(w)[2]<-'so'
x<-droplevels(aggregate(db2bel$bel$cr,db2bel$bel['ut'],length))
colnames(x)<-c('ut','nr')
y<-droplevels(aggregate(db2bel$bel$cr[db2bel$pend],list(db2bel$bel[db2bel$pend,'ut']),length))
colnames(y)<-c('ut','tpnr')

cbind(table(v$pu))
v$pu[v$pu=='M I T PRESS']<-'MIT PRESS'
v$pu[v$pu=='WILEY-BLACKWELL PUBLISHING, INC']<-'WILEY-BLACKWELL'
v$pu[v$pu=='OXFORD UNIV PRESS INC']<-'OXFORD UNIV PRESS'
v<-droplevels(v)
cbind(table(v$pu))

cbind(table(w$so))
levels(w$so)[levels(w$so)=='ECONOMIC HISTORY REVIEW-FIRST SERIES']<-'ECONOMIC HISTORY REVIEW'
w$so[w$so=='ECONOMICA-NEW SERIES']<-'ECONOMICA'
w$so[w$so=='JOURNAL OF SOCIAL FORCES']<-'SOCIAL FORCES'
w<-droplevels(w)
cbind(table(w$so))

soup<-merge(u,v,all=T)
soup<-merge(soup,w,all=T)
soup<-merge(soup,x,all=T)
soup<-merge(soup,y,all=T)
soup<-merge(soup,t,all=T)
soup<-merge(soup,s,all=T)
soup[is.na(soup)]<-0
rm(wok2db_41)

soup$sel<-as.integer(soup$tpnr>0)
soup$rej1<-as.integer(soup$nr==0) #document rejected if no citations
soup$rej2<-as.integer(soup$tpnr==0&soup$nr!=0) #document rejected if only reference was pendant
soup$pg<-as.integer(as.character(soup$pg))

soup$dis<-NA
soup$dis[grep('soci',soup$so,ignore.case=T)]<-'Sociology'
soup$dis[grep('[ck]ono',soup$so,ignore.case=T)]<-'Economics'
soup$dis[grep('anth',soup$so,ignore.case=T)]<-'Anthropology'
soup$dis[grep('poli',soup$so,ignore.case=T)]<-'Political Science'
soup$dis<-factor(soup$dis)

head(soup)
save(soup,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/soup.RData')
}

#### Data for Pivot Tables of Different Populations and Selection effect
if(F){
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/soup.RData')
round(prop.table(table(soup$rej1))*100,3)
round(prop.table(table(soup$rej2))*100,3)
round(prop.table(table(soup$rej1+soup$rej2))*100,3)
round(prop.table(table(soup$sel))*100,3)

seldt<-soup[,c('dt','nr')]
seldt$nr[seldt$nr>1]<-2
round(prop.table(table(seldt),2)*100,3)

if(F) write.table(soup[,-1]
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/journal_tables.tab'
	,quote=F
	,sep='\t'
	,na='.'
	,row.names=F
)

z<-merge(db2bel$bel,v,all.x=T)
z<-merge(z,w)
z<-merge(z,t)
tz<-1/table(z$cr)
z$w<-tz[z$cr]
head(z)

if(F) write.table(z[,3:6]
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/puso_record_table.tab'
	,quote=F
	,sep='\t'
	,na=''
	,row.names=F
)

tpz<-merge(db2bel$bel[db2bel$pend,],v,all.x=T)
tpz<-merge(tpz,w)
tpz<-merge(tpz,t)
ttpz<-1/table(tpz$cr)
tpz$w<-ttpz[tpz$cr]
head(tpz)

if(F) write.table(tpz[,3:5]
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/pusotp_record_table.tab'
	,quote=F
	,sep='\t'
	,na=''
	,row.names=F
)

dis<-soup[,c('dis','tpnr')]
dis$tpnr[dis$tpnr>1]<-1
dis<-table(dis)
dis<-cbind(apply(dis,1,sum),dis[,2])
dis<-rbind(round(prop.table(dis,2)*100,2),Tot=c(NA,NA),HH=c(NA,NA),N=apply(dis,2,sum))
dis['Tot',]<-apply(dis[1:4,],2,sum)
dis['HH',]<-round(apply(dis[1:4,],2,FUN=function(x) sum(x^2)),2)
dis<-cbind(dis,round(apply(dis,1,diff),2))
dis[length(dis)]<-round(dis[length(dis)]/dis[nrow(dis),1],2)
dimnames(dis)<-list(rownames(dis),Records=c('Original','Selected','Δ'))
dis

dt<-soup[,c('dt','tpnr')]
dt$tpnr[dt$tpnr>1]<-1
dt<-table(dt)
dt<-cbind(apply(dt,1,sum),dt[,2])
dt<-rbind(round(prop.table(dt,2)*100,2),Tot=c(NA,NA),HH=c(NA,NA),N=apply(dt,2,sum))
dt['Tot',]<-apply(dt[1:3,],2,sum)
dt['HH',]<-round(apply(dt[1:4,],2,FUN=function(x) sum(x^2)),2)
dt<-cbind(dt,round(apply(dt,1,diff),2))
dt[length(dt)]<-round(dt[length(dt)]/dt[nrow(dt),1],2)
dimnames(dt)<-list(rownames(dt),Records=c('Original','Selected','Δ'))
dt

dtp<-cbind(xtabs(pg~dt,data=soup[,c('dt','pg')])
,xtabs(pg~dt,data=soup[soup$tpnr!=0,c('dt','pg')]))
dtp<-rbind(round(prop.table(dtp,2)*100,2),Tot=c(NA,NA),HH=c(NA,NA),N=apply(dtp,2,sum))
dtp['Tot',]<-apply(dtp[1:3,],2,sum)
dtp['HH',]<-round(apply(dtp[1:4,],2,FUN=function(x) sum(x^2)),2)
dtp<-cbind(dtp,round(apply(dtp,1,diff),2))
dtp[length(dtp)]<-round(dtp[length(dtp)]/dtp[nrow(dtp),1],2)
dimnames(dtp)<-list(rownames(dtp),Pages=c('Original','Selected','Δ'))
dtp

dtz<-cbind(xtabs(w~dt,data=z)
,xtabs(w~dt,data=tpz))
dtz<-rbind(round(prop.table(dtz,2)*100,2),Tot=c(NA,NA),HH=c(NA,NA),N=apply(dtz,2,sum))
dtz['Tot',]<-apply(dtz[1:3,],2,sum)
dtz['HH',]<-round(apply(dtz[1:4,],2,FUN=function(x) sum(x^2)),2)
dtz<-cbind(dtz,round(apply(dtz,1,diff),2))
dtz[length(dtz)]<-round(dtz[length(dtz)]/dtz[nrow(dtz),1],2)
dimnames(dtz)<-list(rownames(dtz),References=c('Original','Selected','Δ'))
dtz

dtzz<-cbind(table(z$dt)
,table(tpz$dt))
dtzz<-rbind(round(prop.table(dtzz,2)*100,2),Tot=c(NA,NA),HH=c(NA,NA),N=apply(dtzz,2,sum))
dtzz['Tot',]<-apply(dtzz[1:3,],2,sum)
dtzz['HH',]<-round(apply(dtzz[1:4,],2,FUN=function(x) sum(x^2)),2)
dtzz<-cbind(dtzz,round(apply(dtzz,1,diff),2))
dtzz[length(dtzz)]<-round(dtzz[length(dtzz)]/dtzz[nrow(dtzz),1],2)
dimnames(dtzz)<-list(rownames(dtzz),Citations=c('Original','Selected','Δ'))
dtzz
}

#handle cr and af pops differently because more than one category per unit, weight required
if(F){
z1<-1/table(z$cr)
z1<-merge(z,data.frame(cr=names(z1),w=as.numeric(z1)))[,2:3]
z1<-aggregate(z1$w,z1['pu'],sum)
colnames(z1)[2]<-'cr'

zz<-droplevels(merge(db2bel$bel,w,all.x=T)[,c('cr','so')])
zz<-zz[!duplicated(z),]
zz1<-1/table(zz$cr)
zz1<-merge(zz,data.frame(cr=names(zz1),w=as.numeric(zz1)))[,2:3]
zz1<-aggregate(zz1$w,zz1['so'],sum)
colnames(zz1)[2]<-'cr'

tpz<-droplevels(merge(db2bel$bel[db2bel$pend,],v,all.x=T)[,c('cr','pu')])
tpz<-tpz[!duplicated(tpz),]
tpz1<-1/table(tpz$cr)
tpz1<-merge(tpz,data.frame(cr=names(tpz1),w=as.numeric(tpz1)))[,2:3]
tpz1<-aggregate(tpz1$w,tpz1['pu'],sum)
colnames(tpz1)[2]<-'tpcr'

tpzz<-droplevels(merge(db2bel$bel[db2bel$pend,],w,all.x=T)[,c('cr','so')])
tpzz<-tpzz[!duplicated(tpz),]
tpzz1<-1/table(tpzz$cr)
tpzz1<-merge(tpzz,data.frame(cr=names(tpzz1),w=as.numeric(tpzz1)))[,2:3]
tpzz1<-aggregate(tpzz1$w,tpzz1['so'],sum)
colnames(tpzz1)[2]<-'tpcr'

(soup2<-list(pucr=z1,socr=zz1,putpcr=tpz1,sotpcr=tpzz1))
for(i in names(soup2)) {
	rn<-as.character(soup2[[i]][,1])
	soup2[[i]]<-cbind(soup2[[i]][,-1])
	p<-round(prop.table(soup2[[i]])*100,2)
	p<-rbind(p,sum(p),sum(soup2[[i]]))
	rownames(p)<-c(rn,'Tot','N')
	colnames(p)<-i
	soup2[[i]]<-p
}
soup2


for(i in names(soup2)){
write.table(soup2[[i]]
	,file=paste('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/',i,'_record_table.tab',sep='')
	,quote=F
	,sep='\t'
	,na=''
	,row.names=T
)
}
}

#weights for journals under multiple publishers
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
zzz<-droplevels(wok2db_41[wok2db_41$field=='CR',c('ut','b')])
colnames(zzz)[2]<-'cr'

zzz<-merge(z,zz)
zzz<-droplevels(zzz[!duplicated(zzz),])
zzzt<-1/table(zzz$cr)
zzz$w<-zzzt[zzz$cr]

write.table(zzz[,-1]
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/puso_citation_table.tab'
	,quote=F
	,sep='\t'
	,na=''
	,row.names=F
)

if(F){
usp<-aggregate(hkl$db[,c('so','pu')],hkl$db['ut'],unique)[,-1]
u<-table(usp)
t<-which(!!u,T)
write.table(data.frame(so=rownames(t),pu=colnames(u)[t[,2]],c=u[t])
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/descriptive.tab'
	,quote=F
	,sep='\t'
	,na='.'
	,row.names=F
)
}




# now we could add authors as a column in the population statistics
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/AFout/db2bel_AF_41.RData')
#follow same procedure as CR
az<-merge(db2bel$bel[,-2],v,all.x=T)
az<-merge(az,w)
az<-merge(az,t)
atz<-1/table(az$zcr)
az$w<-atz[az$zcr]
head(az)

write.table(az[,3:6]
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/AF_record_table.tab'
	,quote=F
	,sep='\t'
	,na=''
	,row.names=F
)

tpaz<-merge(db2bel$bel[db2bel$pend,-2],v,all.x=T)
tpaz<-merge(tpaz,w)
tpaz<-merge(tpaz,t)
ttpaz<-1/table(tpaz$zcr)
tpaz$w<-ttpaz[tpaz$zcr]
head(tpaz)

write.table(tpaz[,3:6]
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/AFtp_record_table.tab'
	,quote=F
	,sep='\t'
	,na=''
	,row.names=F
)


dtaz<-cbind(xtabs(w~dt,data=az)
,xtabs(w~dt,data=tpaz))
dtaz<-rbind(round(prop.table(dtaz,2)*100,2),Tot=c(NA,NA),HH=c(NA,NA),N=apply(dtaz,2,sum))
dtaz['Tot',]<-apply(dtaz[1:3,],2,sum)
dtaz['HH',]<-round(apply(dtaz[1:4,],2,FUN=function(x) sum(x^2)),2)
dtaz<-cbind(dtaz,round(apply(dtaz,1,diff),2))
dtaz[length(dtaz)]<-round(dtaz[length(dtaz)]/dtaz[nrow(dtaz),1],2)
dimnames(dtaz)<-list(rownames(dtaz),Citations=c('Original','Selected','Δ'))
dtaz

write.table(cbind(dtaz,rep(NA,nrow(dt)),dt,rep(NA,nrow(dt)),dtp,rep(NA,nrow(dt)),dtz,rep(NA,nrow(dt)),dtzz)
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/dtAF_record_table.tab'
	,quote=F
	,sep='\t'
	,na=''
	,row.names=T
)

### back to hoffman for string matching
if(F){
for(i in c('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_1.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_2.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_3.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_4.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_5.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_6.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_7.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1980out/wok2db_80.RData')){
	cat('\n',i,sep='')
	load(i)
	cr<-get(grep('wok2db',ls(),value=T))
	rm(list=ls()[!ls()%in%c('i','cr')])
	cr<-levels(droplevels(cr$b[cr$field=='CR']))
	cr<-toupper(cr)
	cr<-sub(', DOI .+','',cr)
	save(cr,file=paste(sub('^(.+/).+','\\1',i),sub('^.+db_([^.]+).+$','cr\\1',i),'.RData',sep=''))
	print(ls())
}

bigcr<-NULL
for(i in c('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_1.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_2.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_3.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_4.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_5.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_6.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/cr15_7.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/cr41.RData'
,'/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1980out/cr80.RData'
)){
	cat('\n',i,sep='')
	load(i)
	bigcr<-unique(c(bigcr,cr))
}
rm(cr)
bigcr<-sort(bigcr)
save(bigcr,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/bigcr.RData')
lsos()

###establish a reasonable baseline similarity threshold
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/bigcr.RData')
split<-strsplit(bigcr,'')
split<-unlist(split)
split<-table(split)
l<-names(split)
p<-split/sum(split)
n<-as.list(sample(nchar(bigcr),1000,replace=F))
blcr<-sapply(n,FUN=function(x) paste(sample(l,size=x,p=p),collapse=''))
library(parallel)
library(stringdist)
s<-sample(bigcr,100)
test<-list()
if(F) {for(i in c('osa', 'lv', 'dl', 'hamming', 'lcs', 'qgram', 'cosine', 'jaccard', 'jw', 'soundex')) {cat('\n',i,sep='');test$t[[i]]<-system.time(test$m[[i]]<-try(stringdistmatrix(s,s,method=i,p=.1)));print(test$m[[i]][1:10,1:10])}
cbind(sort(do.call(rbind,test$t)[,3]))}

blmat<-stringdistmatrix(blcr,blcr,method='jw',p=.1,ncores=detectCores())
blmat<-as.dist(blmat)
(rthresh<-round(quantile(blmat,.025),2))
# 97.5% rounded to 2 is 0.50
act<-sample(bigcr,1000)
s1<-sample(bigcr,1)
rm(bigcr)
act<-stringdistmatrix(act,act,method='jw',p=.1,ncores=detectCores())
act<-as.dist(act)
plot(density(act))
(athresh<-round(quantile(act,.025),2))
# 97.5% rounded to 2 is 0.64

load('/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/1900-1917/db2bel_sets.RData')
sts<-lapply(sets,FUN=function(x) as.dist(stringdistmatrix(x,x,method='jw',p=.1)))
sts<-unlist(sts)
cd<-list(rp=density(blmat),ra=density(act),ma=density(sts))

plot.new()
plot.window(xlim=range(unlist(lapply(cd,FUN=function(x) x$x))),ylim=range(unlist(lapply(cd,FUN=function(x) x$y))))
axis(1)
axis(2)
lines(cd$rp,lty='dotted')
lines(cd$ra,lty='dashed')
lines(cd$ma,lty='solid')

round(mean(blmat<=.2)*100,3)
round(mean(act<=.2)*100,3)
round(mean(sts<=.2)*100,3)
#



s1<-sample(bigcr,1)
if(F) system.time(real<-stringdist(s1,bigcr,method='jw',p=.1))
nc<-nchar(bigcr)
tnc<-table(nc)
for(i in 1:30) {
	cat('\n')
	print(tnc[i])
	tst<-nc==as.integer(names(tnc[i]))
	stst<-sum(tst)
	print(cbind(sample(bigcr[tst],ifelse(stst>5,5,stst))))
}
### criterion, must have at least 3 coded fields (two commas), a year
s<-sample(bigcr,20)
data.frame(s=s,m=regmatches(s,regexpr('(^| )((((17)|(18)|(19))[0-9]{2})|(20[01][0-9]))(,|$)',s)))

date<-grepl('(^| )((((17)|(18)|(19))[0-9]{2})|(20[01][0-9]))(,|$)',bigcr)
commas<-gregexpr(',',bigcr)
w0<-sapply(commas,FUN=function(x) x[1]==-1)
commas<-sapply(commas,length)
commas[w0]<-0
cbind(round(prop.table(table(commas))*100,2))
chco<-data.frame(ch=nc,cm=commas,d=date)
chco$ch[chco$ch>20]<-99
chco$cm[chco$cm>3]<-99
(tchco<-table(chco))
ptchco<-prop.table(tchco)
sum(ptchco[,,'FALSE'],ptchco[-dim(ptchco)[1],1:2,'TRUE'])
nc99<-do.call(seq,as.list(range(quantile(nc,p=c(.005,.995)))))
cull<-!(date==T&nc%in%nc99&commas>1)
mean(cull)
mean(cull)*length(bigcr)
mchco<-cbind(round(prop.table(tchco)*100,5))
mchco<-cbind(mchco,Tot=apply(mchco,1,sum))
(mchco<-rbind(mchco,Tot=apply(mchco,2,sum)))

###culled
bigcr<-bigcr[!cull]
save(bigcr,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/cbigcr.RData')

rm(list=ls())
library(lineprof)
source('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/memtest_bigcr.R')
prof <- lineprof(memtest('max',portions=3))
shine(prof)

rm(list=ls())
library(lineprof)
source('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/memtest_bigcr_.R')
prof <- lineprof(memtest())
shine(prof)

##ooh boy
hoft<-data.frame(matrix(
	c(166.924,119,8,4,16
	,151.334,120,8,4,16
	,175.467,119,11,4,16
	,160.975,120,11,4,16
	,161.404,119,10,4,16
	,148.251,120,10,4,16
	,192.086,119,14,4,16
	,188.702,120,14,4,16
	,176.364,119,13,4,16
	,167.875,120,13,4,16
	,172.729,119,12,4,16
	,156.226,120,12,4,16
	,163.344,119,9,4,16
	,151.283,120,9,4,16 )
,byrow=T,ncol=5))
colnames(hoft)<-c('s','n','p','c','cd')
hoft<-hoft[order(hoft$n,hoft$p),c('n','s','p','c','cd')]
hoft$dc<-factor(apply(hoft[,c('c','cd')],1,paste,collapse='-'))
fit<-lm(s~n+dc+p,hoft)
pl<-split(hoft,hoft[,c('p','dc')],drop=T)
legc<-sort(unique(unlist(hoft['dc'])))
legp<-sort(unique(unlist(hoft['p'])))
pcc<-'black'
#pcc<-rainbow(length(legc),.8)
pcp<-rainbow(length(legp),.8)
#pcp<-gray(1-(1/1:length(legp)),.8)
snakes<-aggregate(hoft$s,hoft[c('n','dc','p')],mean)
snakes<-split(snakes,snakes[c('dc','p')],drop=T)
plot(hoft[,c('n','s')],cex=0)
for(i in 1:length(snakes)) {pl[[i]]<-pl[[i]][order(pl[[i]]$n),];snakes[[i]]<-snakes[[i]][order(snakes[[i]]$n),]}
for(i in 1:length(pl)) {
	points(pl[[i]][,c('n','s')],pch=22,col=pcc[legc==unique(pl[[i]]$dc)],bg=pcp[legp==unique(pl[[i]]$p)],cex=2,lwd=.5)
	lines(snakes[[i]][,c('n','x')],lwd=5,col=pcp[legp==unique(snakes[[i]]$p)],ljoin='mitre',lend='square')
	lines(snakes[[i]][,c('n','x')],lwd=5,col=pcc[legc==unique(snakes[[i]]$dc)],lty='dotted',ljoin='mitre',lend='square')
}
legend('topright',c(paste('P',legp),paste('C',legc)),col=c(pcp,pcc),lwd=5,bty='n')
#text(x=159,y=sapply(pl,FUN=function(x) x[1,2]),labels=names(pl),col=pc)
points(x=hoft$n,y=predict(fit),pch=21,cex=1,col='white',bg=pcc[as.numeric(hoft$dc)])
}

####ooooh doggies
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hoffstringdist/biglist.RData')
cat(paste('c(',apply(cbind(seq(1,length(biglist)-1,(length(biglist)-1)/100),c((seq(1,length(biglist)-1,(length(biglist)-1)/100)-1)[-1],(length(biglist)-1))),1,paste,collapse=':'),')',sep=''),sep='\n')


##### binning bigcr for smooth sailing
if(F){
nb<-nchar(bigcr)
bigcr<-bigcr[order(nb,decreasing=T)]
lbig<-length(bigcr)

bigbox<-array(c(1:lbig,rep(NA,ceiling(lbig^(1/3))^3-lbig)),dim=rep(ceiling(lbig^(1/3)),3))
bigbox[1:10,1:10,1:10]
bigbox[111:120,111:120,111:120]
biglist<-apply(bigbox,3,c)
biglist<-apply(biglist,1,list)
biglist<-lapply(biglist,FUN=function(x) {x<-na.omit(unlist(x));attributes(x)<-NULL;x})
bignc<-sapply(biglist,FUN=function(x) nchar(bigcr[x]))
bigsd<-sapply(bignc,sd)
plot(1:length(bigsd),bigsd,type='l')
bigm<-sapply(bignc,mean)
plot(1:length(bigm),bigm,type='l')
bigt<-sapply(bignc,sum)
plot(1:length(bigt),bigt,type='l')
bigl<-sapply(bignc,length)
cbind(table(bigl))

biglist<-rev(biglist)
names(biglist)<-paste('batch',1:length(biglist),sep='')
biglist<-c(cr=list(bigcr),biglist)
save(biglist,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hoffstringdist/biglist.RData')

bigpor<-list(cr=bigcr)
for(i in 8:12) {
	j<-paste('por',i,sep='')
	bigpor[[j]]<-matrix(c(1:lbig,rep(NA,(i*ceiling(lbig/i))-lbig)),ncol=i,byrow=T)
	bigpor[[j]]<-apply(bigpor[[j]],2,list)
	bigpor[[j]]<-lapply(bigpor[[j]],FUN=function(x) {x<-na.omit(unlist(x));attributes(x)<-NULL;x})
}
save(bigpor,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/hoffstringdist/bigpor.RData')
}

#hrm, fuzzy matching broken
if(F){
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')

wok2db_41<-wok2db_41[wok2db_41$field%in%c('AF'),]
colnames(wok2db_41)<-c('ut','field','b')
rownames(wok2db_41)<-NULL
r<-droplevels(wok2db_41[wok2db_41$field=='AF',c('ut','b')])
colnames(r)[2]<-'af'
rm(wok2db_41)

require(stringdist)
#sd<-unique(head(toupper(r$af),100))
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=1))
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=3))
#sd<-unique(head(toupper(r$af),1000))
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=1))
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=3))
#sd<-toupper(unique(r$af))
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=1)) # 60 seconds through terminal
#system.time(mat<-adist(sd)) # really long apparently # 335 seconds through terminal
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=3)) # 50 seconds
#system.time(mat<-stringdistmatrix(sd,sd,method='jw',p=0.1,ncores=4)) # 32 seconds!!
#save(mat,file='mat.RData')
}

### plot of selection effect over time
if(F){
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/soup.RData')
sef<-array(0,dim=c(length(unique(soup$py)),length(unique(soup$dt)),2),dimnames=list(sort(levels(soup$py)),unique(soup$dt),c('Original','Selected')))
tsf<-table(soup[,c('py','dt')])
sef[rownames(tsf),colnames(tsf),'Original']<-tsf
tsf<-table(soup[soup$tpnr!=0,c('py','dt')])
sef[rownames(tsf),colnames(tsf),'Selected']<-tsf

all<-apply(sef[,,'Selected'],1,sum)/apply(sef[,,'Original'],1,sum)
xv<-as.integer(names(all))
xl<-range(xv)

sot<-aggregate(as.integer(as.character(soup$py)),soup['so'],min)
colnames(sot)[1:2]<-c('so','beg')
sot<-data.frame(sot,end=aggregate(as.integer(as.character(soup$py)),soup['so'],max)$x,y=all[as.character(sot$beg)])
.simpleCap <- function(x) {
    s <- strsplit(x, ' ')[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = '', collapse = ' ')
}
sot$so<-as.character(sot$so)
sot$so<-sub(' Fur ',' fur ',sub(' Of ',' of ',sapply(as.list(sot$so),FUN=function(x) .simpleCap(tolower(x)))))
sot<-aggregate(sot,sot['y'],unique)
sot$so<-sapply(sot$so,paste,collapse='; ')
cx<-.5
pdf('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/selection.pdf',h=(11-2)/3,w=8.5-2)
par(omi=rep(0,4),mai=c(.5,.5,.5,0))
plot.new()
plot.window(ylim=range(all),xlim=xl)
sot$pos<-c(2,4)[as.integer((strwidth(sot$so,cex=cx)+sot$beg)<xl[2])+1]
#lines(xl,sef[,'Article','Selected']/sef[,'Article','Original'])
#rev<-sef[,'Review','Selected']/sef[,'Review','Original']
#rev<-rev[!is.na(rev)]
#lines(as.integer(names(rev)),rev)
#lines(xl,sef[,'Book Review','Selected']/sef[,'Book Review','Original'])
for(i in 1:nrow(sot)) text(sot$beg[i],sot$y[i],labels=sot$so[i],col='darkgray',pos=sot$pos[i],cex=cx)
points(sot$beg,sot$y,pch=19,col='black')
lines(xv,all)
axis(1,at=seq(xl[1],xl[2],by=5),labels=seq(xl[1],xl[2],by=5))
axis(2,las=1)
title(main='Selection Effect Over Time',xlab='',ylab='')
dev.off()
}

### random Herfindahl indices
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/soup.RData')
rherf<-replicate(1000,sum((table(sample(unique(soup$pu),10000,replace=T))/100)^2))

### UT permutations
##publisher
permPUUT<-replicate(1000,prop.table(table(soup$pu[sample(1:43736,14016)]))*100)
psdpuut<-apply(permPUUT,1,sd)
psdpuut<-psdpuut[order(apply(permPUUT,1,mean),decreasing=T)]
cat(round(psdpuut,3),sep='\n')
#test normality
puUTnorm<-apply(permPUUT,1,FUN=function(x) {qqnorm(x);qqline(x,col=2);shapiro.test(x)})
#random Herfindahl
rherfpuut<-apply(permPUUT,2,FUN=function(x) sum(x^2))
sd(rherfpuut)
{qqnorm(rherfpuut);qqline(rherfpuut,col=2);shapiro.test(rherfpuut)}
##journal
##publisher
permSOUT<-replicate(1000,prop.table(table(soup$so[sample(1:43736,14016)]))*100)
psdsout<-apply(permSOUT,1,sd)
psdsout<-psdsout[order(apply(permSOUT,1,mean),decreasing=T)]
cat(round(psdsout,3),sep='\n')
#test normality
soUTnorm<-apply(permSOUT,1,FUN=function(x) {qqnorm(x);qqline(x,col=2);shapiro.test(x)})
#random Herfindahl
rherfsout<-apply(permSOUT,2,FUN=function(x) sum(x^2))
sd(rherfsout)
{qqnorm(rherfsout);qqline(rherfsout,col=2);shapiro.test(rherfsout)}



sum((^2))
95/sd(rherfUT)
### at 1000 replications, difference is 7.104764 sds above the random deviations from the original distribution

##### degree distribution of each mode
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/mel2net_41.RData')
library(network)


uth<-5 #set upper threshold for summary table
ddt<-list()
ddt$bmut<-c('0'=sum(wok2db_41$b[wok2db_41$field=='NR']=='0'),table(table(db2bel$bel$ut)))
ddt$tput<-table(table(db2bel$bel$ut[db2bel$pend]))
ddt$ut<-table(degree(mel2net$tput,gmode='graph'))
ddt$bmcr<-table(table(db2bel$bel$cr))
ddt$tpcr<-table(table(db2bel$bel$cr[db2bel$pend]))
ddt$cr<-table(degree(mel2net$tpcr,gmode='graph'))
br<-paste('>=',uth+1,sep='')
dd<-matrix(nrow=uth+2,ncol=6,dimnames=list(degree=c(0:uth,br),mode=c('bmut','tput','ut','bmcr','tpcr','cr')))
for(i in names(ddt)) {
	nms<-names(ddt[[i]])[names(ddt[[i]])%in%as.character(0:uth)]
	dd[nms,i]<-ddt[[i]][nms]
	dd[br,i]<-sum(ddt[[i]][as.integer(names(ddt[[i]]))>uth])
}
dd[is.na(dd)]<-0
dd<-round(prop.table(dd,2)*100,2)
(dd<-rbind(dd,Tot=apply(dd,2,sum),N=sapply(ddt,sum),sapply(ddt,function(x) tail(rep(names(x),x),5))))
write.table(dd
	,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/descriptive/degreedist.tab'
	,quote=F
	,sep='\t'
	,na='.'
	,row.names=T
)
rm(db2bel,mel2net,wok2db_41)

##### which clusters are interesting
ksig<-readLines('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/stata/ksig.tab')
ksig<-grep('^ +k[0-9]',ksig,value=T)
ksig<-do.call(rbind,strsplit(ksig,'[ |]+'))
nks<-ksig[,2]
ksig<-data.frame(lapply(data.frame(matrix(
	cbind(sub('^k([0-9]+).+','\\1',ksig[,2]),ksig[,3:8])
	,dimnames=list(ksig[,2],c('k','hr','se','z','p','l95','r95'))
	,nrow=nrow(ksig)
	),stringsAsFactors=F),type.convert))
rownames(ksig)<-nks

pdf('ksig.pdf')
plot(ksig[,c('z','hr')],cex=5,col=hsv(1,1,1,.5))
text(ksig[,c('z','hr')],labels=rownames(ksig),cex=.1,col=gray(0,.5))
abline(h=1,v=0)
dev.off()

#####effects of network selection
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hkl.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/db2bel_41.RData')
sel<-list(ut0=NULL,ut1=NULL,ut2=NULL,ut3=NULL,ut4=NULL,cr0=NULL,cr1=NULL,cr2=NULL,cr3=NULL,cr4=NULL)

#criterion 0: original coding, out because journal not selected, and because coverage of journal is spotty
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')
sel$ut0<-levels(droplevels(wok2db_41$id.))
rm(wok2db_41)
sel$cr0<-unique(db2bel$bel$cr)

##levs are all of the CRs that were fed into CFinder as exported from the single mode CR node UT edge network
#so why isn't lev the length of the unique CRs in bel2mel? Because isolates were dropped. Where does an isolate come from? A ut with only one reference left after pendants are dropped.

#criterion 1: out because of projections as bimodal UT to CR list, ut's with zero references excluded
sel$ut1<-unique(db2bel$bel$ut)
sel$cr1<-unique(db2bel$bel$cr)

#criterion 2: out because of cr size threshold

sel$cr2<-hkl$mlev[unique(unlist(hkl$memb))]
sel$ut2<-unique(db2bel$bel$ut[db2bel$bel$cr%in%sel$cr2])

#criterion 3: out because of unimodal projection: cr only exists in UT where it is the only reference
sel$cr3<-readLines('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/tpcrel_1941_CFinder_levels.txt')
sel$ut3<-unique(db2bel$bel$ut[db2bel$bel$cr%in%sel$cr3])

#criterion 4: reintroduce isolates
sel$cr4<-hkl$mlev
sel$ut4<-unique(db2bel$bel$ut[db2bel$bel$cr%in%sel$cr4])


#criterion 5: out because cr not cohesive, member of 3-clique
sel$cr5<-hkl$mlev[unique(hkl$db$cr)]
sel$ut5<-unique(db2bel$bel$ut[db2bel$bel$cr%in%sel$cr5]) # dbout has every ut from criterion 1 and 2
#note that nothing is eliminated

if(T){
cat('\nVerify zero overlap between selection categories:')
com<-combn(names(sel),2)
for(i in 1:ncol(com)) cat('\n',length(intersect(sel[[com[1,i]]],sel[[com[2,i]]])),' b/w ',com[1,i],' and ',com[2,i],sep='')
com<-matrix(sapply(sel,length),ncol=2)
colnames(com)<-c('ut','cr')
rownames(com)<-c('C0','C1','C2','C3','C4','C5')
cat('\n\n')
com<-data.frame(com,round(rbind(apply(com,2,diff),c(NA,NA))/com*100,1)
,'note'=c('# cut 0 cr'
,'# cut pendants'
,paste('# cut k-clique<',hkl$thresh,sep='')
,paste('# cut because not in k>=',min(na.omit(as.integer(sub('(k([0-9]*))?.*','\\2',names(hkl$memb))))),'-clique',sep='')
,'# final selection'))
print(com)
sapply(sel,FUN=function(x) sum(table(x)>1))

cat('\n\nExclusive sets\n\n')
esel<-sel
esel$ut0<-setdiff(esel$ut0,esel$ut1)
esel$ut1<-setdiff(esel$ut1,esel$ut2)
esel$ut2<-setdiff(esel$ut2,esel$ut3)
esel$ut3<-setdiff(esel$ut3,esel$ut4)
esel$cr0<-setdiff(esel$cr0,esel$cr1)
esel$cr1<-setdiff(esel$cr1,esel$cr2)
esel$cr2<-setdiff(esel$cr2,esel$cr3)
esel$cr3<-setdiff(esel$cr3,esel$cr4)

cat('\nVerify zero overlap between selection categories:')
ecom<-combn(names(esel),2)
for(i in 1:ncol(ecom)) cat('\n',length(intersect(esel[[ecom[1,i]]],esel[[ecom[2,i]]])),' b/w ',ecom[1,i],' and ',ecom[2,i],sep='')
ecom<-matrix(sapply(esel,length),ncol=2)
colnames(ecom)<-c('ut','cr')
rownames(ecom)<-c('C0','C1','C2','C3','C4','C5')
cat('\n\n')
print(data.frame(ecom,apply(ecom,2,cumsum),'note'=c('# first cut for 0 references','# cut when pendants dropped','# cut as isolates in unimodal','# cut because lower than triad density','# reintroduce isolates','# final selection')))
table(table(unlist(esel)))
}


########## Modeling ##########
rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/duration.RData')

#####Examine high correlations
hdm1<-hclust(as.dist(1-abs(cor(d[sapply(d,class)!='factor'&colnames(d)!='gap']))))
pdf('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/stata/corr.pdf')
plot(hdm1,xlab='',main='1-abs(corr()) Clusters',cex=.75)
hd<-hclust(as.dist(abs(cor(d[sapply(d,class)!='factor'&colnames(d)!='gap']))))
plot(hd,xlab='',main='abs(corr()) Clusters',cex=.75)
dev.off()

#####ooooh! from stata haz ratio estimates for years ending and years beginning
start<-read.delim('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/stata/start.txt',header=F)
end<-read.delim('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/stata/end.txt',header=F)

y<-range(c(start[,-1],end[,-1]))
plot(start[,c(1,2)],ylim=y,type='l',lwd=2,ylab='Haz Ratio',xlab='',main='Communities Starting in Year')
lines(start[,c(1,3)],col='black',lty=2,lwd=2)
lines(start[,c(1,4)],col='black',lty=2,lwd=2)
abline(h=1)

plot(end[,c(1,2)],ylim=y,col='red',type='l',lty=1,lwd=2,ylab='Haz Ratio',xlab='',main='Communities Failing in Year')
lines(end[,c(1,3)],col='red',lty=2,lwd=2)
lines(end[,c(1,4)],col='red',lty=2,lwd=2)
abline(h=1)

###########
g<-gfpq$g
el<-list()
for(i in 1:length(g)) {
	el[[i]]<-matrix(do.call(rbind,lapply(strsplit(g[[i]][8:length(g[[i]])],' '),as.numeric)),ncol=3)
	if(!is.na(el[[i]])) el[[i]]<-el[[i]][order(el[[i]][,1],el[[i]][,2]),]
}
names(el)<-names(g)
dup<-duplicated(el,fromLast=T)
el<-el[!dup]
sdis<-sdis[!dup]
library(network)
n<-list()

for(i in 1:length(el)) if(class(el[[i]])=='numeric') {n[[i]]<-network(matrix(el[[i]][1:2],ncol=2),matrix.type='edgelist',directed=F)} else {n[[i]]<-network(cbind(as.character(el[[i]][,1]),as.character(el[[i]][,2])),matrix.type='edgelist',directed=F)}

pdf('k-cliques.pdf')
	for(i in 1:(length(n)-1)) plot(n[[i]],edge.lwd=data.frame(el[[i]])[[3]],vertex.cex=1)
dev.off()

##############################
#### BACK TO BEL2MEL NETS ####
##############################

##nah

##see how k-qliques are nested
load('gfpq.RData')
p<-gfpq$p

for(i in 1:length(p)) for(j in 1:length(p[[i]])) p[[i]][[j]]<-lapply(p[[i]][[j]],sort)
nest<-list()
for(i in length(p):2) for(j in 1:length(p[[i-1]])){
	x<-paste(names(p[i]),names(p[i-1]),sep='->')
	y<-paste('->l',j-1,sep='')
	cat('\r',x,y,'\t\t')
	nest[[x]][[y]]<-sapply(lapply(p[[i]],'%in%',p[[i-1]][[j]]),mean) #proportion of the smaller, higher level cluster that is contained in the larger, lower level clusters. Should also reveal boundaries within high level clusters.
	names(nest[[x]][[y]])<-paste(names(nest[[x]][[y]]),'->',sep='')
}
for(i in 1:length(nest)) {
	if(!is.list(nest[[i]])) nest[[i]]<-as.list(nest[[i]])
	nest[[i]]<-do.call(rbind,nest[[i]])
}
nest_41<-nest #these are interk community relationships by edgelist
save(nest_41,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/nest_41.RData')

for(i in 1:length(nest)) colnames(nest[[i]])<-paste('l',0:(dim(nest[[i]])[2]-1),'->',sep='')
for(i in 1:length(f)) names(f[[i]])<-paste('l',0:(length(f[[i]])-1),sep='')

rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/gfpq.RData')
f<-gfpq$f

fnest<-list() #interk community relationships by nodelist
fnest$down<-list()
fnest$up<-list()

#DOWN
for(i in length(f):2) for(j in 1:length(f[[i-1]])){
	x<-paste(names(f[i]),names(f[i-1]),sep='->')
	y<-paste('->l',j-1,sep='')
	cat('\r',x,y,'\t\t')
	fnest$down[[x]][[y]]<-sapply(lapply(f[[i]],'%in%',f[[i-1]][[j]]),mean) #proportion of the larger, lower level cluster that is contained in the smaller, higher level clusters. Should also reveal boundaries within high level clusters.
	names(fnest$down[[x]][[y]])<-paste(names(fnest$down[[x]][[y]]),'->',sep='')
}
for(i in 1:length(fnest$down)) {
	if(!is.list(fnest$down[[i]])) fnest$down[[i]]<-as.list(fnest$down[[i]])
	fnest$down[[i]]<-do.call(rbind,fnest$down[[i]])
	colnames(fnest$down[[i]])<-paste('l',0:(dim(fnest$down[[i]])[2]-1),'->',sep='')
}

#UP
for(i in (length(f)-1):1) for(j in 1:length(f[[i+1]])){
	x<-paste(names(f[i]),names(f[i+1]),sep='<-')
	y<-paste('<-l',j-1,sep='')
	cat('\r',x,y,'\t\t')
	fnest$up[[x]][[y]]<-sapply(lapply(f[[i+1]][[j]],'%in%',f[[i]]),mean) #proportion of the smaller, higher level cluster that is contained in the larger, lower level clusters. Should also reveal boundaries within high level clusters.
	names(fnest$up[[x]][[y]])<-paste(names(fnest$up[[x]][[y]]),'<-',sep='')
}
for(i in 1:length(fnest$up)) {
	if(!is.list(fnest$up[[i]])) fnest$up[[i]]<-as.list(fnest$up[[i]])
	fnest$up[[i]]<-do.call(rbind,fnest$up[[i]])
	colnames(fnest$up[[i]])<-paste('l',0:(dim(fnest$up[[i]])[2]-1),'->',sep='')
}

fnest_41<-fnest #these are interk community relationships by edgelist
save(fnest_41,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/fnest_41.RData')

fnest_41<-lapply(fnest_41,round,3)
colnames(fnest_41[[1]])<-'l0->'
kel<-list()
for(i in 1:length(fnest_41)){
	n<-unlist(strsplit(names(fnest_41[i]),'->'))
	t<-as.logical(fnest_41[[i]])
	a<-sub('->','',paste(n[1],colnames(fnest_41[[i]]),sep=''))
	b<-sub('->','',paste(n[2],rownames(fnest_41[[i]]),sep=''))
	kel[[i]]<-cbind(expand.grid(b,a,KEEP.OUT.ATTRS=F,stringsAsFactors=F),ew=unlist(fnest_41[i]))[t,c(2,1,3)]
}
kel<-do.call(rbind,kel)
rownames(kel)<-NULL
kel<-kel[order(as.character(kel[[1]]),as.character(kel[[2]])),]
kelt<-network(kel[,1:2],matrix.type='edgelist')
set.edge.attribute(kelt,'ew',kel[[3]])
kelt%e%'ew'<-kel[[3]]

#how to visualize nesting
lapply(f,sapply,length) #gives size of each community (many nodes counted multiple times)


if(F){

	load('1941out/by_year/thatgirlis_41_by.RData')
	trend<-list()
	for(i in names(thatgirlis_41_by)) for(j in c('0','1','2')) trend[[j]][[i]]<-c(year=as.integer(i),t=thatgirlis_41_by[[i]][j,'t'])
	for(j in c('0','1','2')) trend[[j]]<-do.call(rbind,trend[[j]])
	plot.new()
	plot.window(
		xlim=range(do.call(rbind,trend)[,'year'])
		,ylim=range(do.call(rbind,trend)[,'t'])
	)
	lines(trend[['0']],type='l',lty=3)
	lines(trend[['1']],type='l',lty=2)
	lines(trend[['2']],type='l',lty=1)
	abline(h=0,col='red')
	axis(1)
	axis(2)
}



#gini test

rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')
cr<-droplevels(wok2db_41[wok2db_41$field=='CR',c('id.','b')])
py<-droplevels(wok2db_41[wok2db_41$field=='PY',c('id.','b')])
crpy<-merge(cr,py,by='id.')
names(crpy)<-c('ut','cr','py')
crpyby<-split(crpy$cr,crpy$py,drop=T)
crpyby<-lapply(crpyby,droplevels)
crpyby<-crpyby[order(names(crpyby))]
crpyby<-lapply(crpyby,table)
count<-sapply(crpyby,sum)
ucount<-sapply(crpyby,length)
library(ineq)
gini<-data.frame(py=as.numeric(names(crpyby)),g=sapply(crpyby,ineq,type='Gini'),c=count,uc=ucount)
for(i in names(crpyby)) crpyby[[i]]<-crpyby[[i]][crpyby[[i]]!=1]
gini_1<-data.frame(py=as.numeric(names(crpyby)),g=sapply(crpyby,ineq,type='Gini'),c=count,uc=ucount)
g41<-list(gini=gini,gini_1=gini_1)
save(g41,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/g41.RData')

rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1980out/wok2db_80.RData')
cr<-droplevels(wok2db_80[wok2db_80$field=='CR',c('id.','b')])
py<-droplevels(wok2db_80[wok2db_80$field=='PY',c('id.','b')])
crpy<-merge(cr,py,by='id.')
names(crpy)<-c('ut','cr','py')
crpyby<-split(crpy$cr,crpy$py,drop=T)
crpyby<-lapply(crpyby,droplevels)
crpyby<-crpyby[order(names(crpyby))]
crpyby<-lapply(crpyby,table)
count<-sapply(crpyby,sum)
ucount<-sapply(crpyby,length)
library(ineq)
gini<-data.frame(py=as.numeric(names(crpyby)),g=sapply(crpyby,ineq,type='Gini'),c=count,uc=ucount)
for(i in names(crpyby)) crpyby[[i]]<-crpyby[[i]][crpyby[[i]]!=1]
gini_1<-data.frame(py=as.numeric(names(crpyby)),g=sapply(crpyby,ineq,type='Gini'),c=count,uc=ucount)
g80<-list(gini=gini,gini_1=gini_1)
save(g80,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1980out/g80.RData')

rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_1.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_2.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_3.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_4.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_5.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_6.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_7.RData')

wok2db_15<-rbind(wok2db_15_1,wok2db_15_2,wok2db_15_3,wok2db_15_4,wok2db_15_5,wok2db_15_6,wok2db_15_7)
rm(list=ls()[grep('wok2db_15_',ls())])
cr<-droplevels(wok2db_15[wok2db_15$field=='CR',c('id.','b')])
py<-droplevels(wok2db_15[wok2db_15$field=='PY',c('id.','b')])
crpy<-merge(cr,py,by='id.')
names(crpy)<-c('ut','cr','py')
crpyby<-split(crpy$cr,crpy$py,drop=T)
crpyby<-lapply(crpyby,droplevels)
crpyby<-crpyby[order(names(crpyby))]
crpyby<-lapply(crpyby,table)
count<-sapply(crpyby,sum)
ucount<-sapply(crpyby,length)
library(ineq)
gini<-data.frame(py=as.numeric(names(crpyby)),g=sapply(crpyby,ineq,type='Gini'),c=count,uc=ucount)
for(i in names(crpyby)) crpyby[[i]]<-crpyby[[i]][crpyby[[i]]!=1]
gini_1<-data.frame(py=as.numeric(names(crpyby)),g=sapply(crpyby,ineq,type='Gini'),c=count,uc=ucount)
g15<-list(gini=gini,gini_1=gini_1)
save(g15,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/g15.RData')


rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/g41.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1980out/g80.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/g15.RData')
gini<-rbind(g41$gini,g80$gini,g15$gini)[,c('py','g')]
#gini_1<-rbind(g41$gini_1,g80$gini_1,g15$gini_1) #[,c('py','g')]

cor(gini)
#cor(gini_1)

require(graphics)

plot(gini)
#plot(gini_1)
abline(v=seq(1900,2015,10),col='lightgray')

#points(gini,pch=19)

#(z <- lm(g~py,gini_1))
(x <- lm(g~py,gini))
#lz<-lowess(gini_1,f=1/10)
lx<-lowess(gini,f=1/10)

#abline(coef(z),lty=2)
abline(coef(x),lty=1)
#lines(lz,lty=2)
lines(lx,lty=1)

## Tukey-Anscombe Plot :
#plot(residuals(z) ~ fitted(z), main = deparse(z$call))


library(rdd)
mod<-list()
for(i in seq(1900.5,2000,1)) mod[[as.character(i)]]<-try(RDestimate(g~py,data=gini,cutpoint=i,verbose=T))
for(i in which(sapply(mod,class)=='RD')) try(plot(mod[[i]],range=c('min','max')))

library(scatterplot3d)
scatterplot3d(gini$py,gini$g,gini$uc)
library(rgl)
plot3d(gini$py,gini$uc,gini$g)

ggplotRegression <- function (fit) {

require(ggplot2)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
  geom_point() +
  stat_smooth(method = 'lm', col = 'red') +
  opts(title = paste('Adj R2 = ',signif(summary(fit)$adj.r.squared, 5),
                     '; Intercept =',signif(fit$coef[[1]],5 ),
                     '; Slope =',signif(fit$coef[[2]], 5),
                     '; P =',signif(summary(fit)$coef[2,4], 5)))
}

#batch db
source('/u/home/b/bambrose/clique/dissertation_source.R')
inp<-list.files('/u/home/b/bambrose/clique/1981-2015',recursive=T,full.names=T)
tin<-table(sub('(.+/)(.+)(/.+$)','\\2',inp))
b<-cumsum(tin)%/%round(sum(tin)/7)
b<-split(names(b),b)

rm(list=ls())
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1980out/wok2db_80.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_1.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_2.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_3.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_4.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_5.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_6.RData')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/2015out/wok2db_15_7.RData')
allcr1900_2015<-rbind(wok2db_15_1,wok2db_15_2,wok2db_15_3,wok2db_15_4,wok2db_15_5,wok2db_15_6,wok2db_15_7,wok2db_41,wok2db_80)

if(F){
	rm(list=ls())
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/glmfpq.RData")

	g<-glmfpq$g
	g<-g[as.logical(sapply(g,length)!=7)]
	for(i in 1:length(g)) g[[i]]<-g[[i]][-(1:7)]
	g<-g[!duplicated(g,fromLast=T)]
	g<-lapply(lapply(lapply(g,strsplit," "),sapply,rbind),t)

	q<-glmfpq$q
	q<-q[as.logical(sapply(q,length)!=7)]
	for(i in 1:length(q)) q[[i]]<-q[[i]][-(1:7)]
	q<-q[!duplicated(q,fromLast=T)]
	for(i in names(q)) q[[i]]<-q[[i]][-length(q[[i]])]
	q<-lapply(lapply(lapply(q,strsplit," "),sapply,rbind),t)

	f<-glmfpq$f

	library(statnet)
	library(latentnet)

	n<-list()
	c<-list() #count distributions
	for(i in names(g)) {
		n[[i]]<-network(g[[i]][,1:2],matrix.type="edgelist",directed=F)
		set.edge.attribute(n[[i]],attrname="ew",value=as.integer(g[[i]][,3]))
		c[[i]]<-c("0"=sum(!n[[i]][,]),table(n[[i]]%e%"ew"))
		add.isolates(n[[i]],n=length(f[[i]])-network.size(n[[i]]))
		c[[i]]<-list(
			"o"=c[[i]]
			,"e.mean"=table(rpois(sum(c[[i]]),mean(rep(c[[i]],names(c[[i]])))))
			,"e.var"=table(rpois(sum(c[[i]]),var(c[[i]]))))
	}
	save(n,file="h.com.nets1941.RData")

	pdf("com_graph.pdf")
	for(i in rev(names(n))){
		plot(n[[i]]
				 ,displaylabels=T
				 ,label.pos=5
				 ,label.cex=.5
				 ,vertex.col="white"
				 ,vertex.border="gray"
				 ,vertex.cex=1.5
				 ,edge.lwd="ew"
				 ,main=i)
	}
	dev.off()

	## must load write.ergmm

	for(i in 1:5){
		write.ergmm(where="/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hoff"
								,dat="h.com.nets1941.RData"
								,net="n"
								,groups=i
								,dimensions=1
								,verbosity=4
		)}
	setwd("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hoff")
	n<-n[paste("k",20:23,sep="")]

	###report out logs
	setwd("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/hoff")
	fls<-grep("write.ergmm.+out",list.files(),value=T)
	r<-list()
	for(i in fls){
		r$raw[[i]]<-readLines(i)
	}
	names(r$raw)<-sub(".+_([^.]+).+","\\1",names(r$raw))
	r$beg<-lapply(r$raw,grep,pattern="^<")
	names(r$beg)<-names(r$raw)
	r$end<-lapply(r$beg,"-",2)
	for(i in names(r$end)) r$end[[i]]<-c(r$end[[i]][-1],length(r$raw[[i]])-4)
	r$sub<-list()
	for(i in names(r$raw)) for(j in 1:length(r$beg[[i]])) r$sub[[i]][[sub("<+ ([^ ]+) >+","\\1",r$raw[[i]][r$beg[[i]]][j])]]<-r$raw[[i]][r$beg[[i]][j]:r$end[[i]][j]]

	s<-array(NA,dim=c(length(names(r$sub[[1]])),length(names(r$sub)),2),dimnames=list(k=names(r$sub[[1]]),G=names(r$sub),att=c("fin","bic")))
	s[,,"fin"]<-grepl("Finished.",unlist(r$raw)[unlist(r$end)])

	#get bic from hoffman!
	load(tail(grep("fit2bic.+RData",list.files(),value=T),1))
	s[,,"bic"]<-as.numeric(unlist(fit2bic))
	save(s,file="s.RData")

	#plot bics
	load("s.RData")
	pdf("bic.pdf",h=22,w=4)
	plot.new()
	plot.window(xlim=c(1,5), ylim=log(range(s[,,"bic"],na.rm=T)))
	text(labels=rownames(s)[!is.na(s[,1,"bic"])],x=.925,y=log(s[!is.na(s[,1,"bic"]),1,"bic"]),cex=.5)
	abline(v=1:5,col="lightgray")
	apply(log(s[!is.na(s[,1,"bic"]),,"bic"]),1,lines)
	axis(1)
	axis(2)
	title(main="BIC for Groups in k3:k50")
	title(xlab="G")
	title(ylab="BIC")
	box()
	dev.off()

	#choosing cross section based on low variance convention
	(plt<-cbind("#"=sapply(glmfp[["f"]],length),sd=sapply(sapply(glmfp[["f"]],sapply,length),sd)))
	pdf("community_sizebysd.pdf")
	plot(plt,xlab="# communites",ylab="SD community size",type="n",main=c("Variation in Community Size","by Number of Communities"))
	lines(plt,col="lightgray")
	text(plt,labels=names(glmfp[["f"]]),cex=.45,offset=0,pos=4)
	dev.off()

	(k8<-table(sapply(glmfp[["f"]]$k8,length)))
	(k9<-table(sapply(glmfp[["f"]]$k9,length)))
	nam<-sort(unique(as.integer(c(names(k8),names(k9)))))
	kcom<-matrix(nrow=length(nam),ncol=2,dimnames=list(nam,c("k8","k9")))
	kcom[names(k8),"k8"]<-k8
	kcom[names(k9),"k9"]<-k9
	cbind(kcom,apply(kcom,1,diff))

	d8<-density(sapply(glmfp[["f"]]$k8,length),bw=1.822)
	d9<-density(sapply(glmfp[["f"]]$k9,length),bw=1.822)
	plot(d8,ylim=range(d8$y,d9$y),xlim=range(d8$x,d9$x),lty=2)
	lines(d9)

	#k9 appears to be better because in k8 variance is balanced by inclusion of one big community and a bunch of minimum communities


	###########################
	# match dates to clusters #
	###########################
	#require loaded bel2mel, wok2db
	rm(list=ls())
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData")
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/bel2mel_41.RData")
	library(data.table)
	py<-data.table(droplevels(wok2db_41[wok2db_41$field=="PY",-2]),key="b.ind.")
	m<-list()
	system.time(for(i in 1:dim(bel2mel$tpcrel)[1]) m[i]<-as.integer(as.character(py[bel2mel$tpcrel$x[[i]]]$b)))
	bel2mel$tpcrel$py<-m
	rm(m)

	##match cluster ids
	load("/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/glmfpq.RData")
	l<-glmfpq$l[["k9"]]
	rm(glmfpq)

}

### hoffZCRcompiler.R
if(F){

	if(F){ rm(list=ls())
				 require(data.table)
				 long<-data.table(k=list())
	}

	all<-c(
		grep("similars",list.files("/Volumes/Iomega_HDD/go",recursive=T,full.name=T),value=T)
		,grep("similars",list.files("/Volumes/Iomega_HDD/miss",recursive=T,full.name=T),value=T)
	)
	bt<-table(type.convert(sub(".+_Sbatch([^_]+)_.+","\\1",all)))
	btn<-as.integer(names(bt))
	which(!(1:max(btn))%in%btn)
	bt[bt!=10]

	cat("\n")
	t0<-proc.time()

	for(h in 1:max(btn)){
		batch<-grep(paste("_Sbatch",h,"_",sep=""),all,value=T)
		### first take all partitions and combine horizontally

		load(batch[1])
		dat<-data.table(k=savedsim)
		for(i in 2:length(batch)) {
			load(batch[i])
			dat[,paste("P",i,sep=""):=savedsim]
		}

		dat<-apply(dat,1,FUN=function(x) {x<-unlist(list(x));names(x)<-NULL;x})
		names(dat)<-names(savedsim)
		dat<-lapply(dat)
		### eliminate zero length strings

		long[[length(long)+1]]<-dat
		cat("\rFinished with",h,"or",round(h/max(btn)*100,3),"% complete")
	}
	t1<-proc.time()
	t1-t0
	long<-do.call(c,long)
	long<-data.table(k=names(long),set=long,key="k")

	#save(long,file="/Volumes/Iomega_HDD/long.RData")

	### remove "[ANONYMOUS]"
	t0<-proc.time()

	t1<-proc.time()
	t1-t0

	anon<-grepl("\\[ANONYMOUS\\],? ?",long[,k])
	alarm()
	long[anon,set:=lapply(set,FUN=function(x) sub("\\[ANONYMOUS\\],? ?","",x))]
	long[anon,set:=lapply(set,FUN=function(x) x[nchar(x)>13])]
	long[anon,k:=sub("\\[ANONYMOUS\\],? ?","",k)]
	t2<-proc.time()
	t2-t1

	cull<-.1

	#load cull.f

	if(F){t<-data.table(a=replicate(20,paste0(sample(letters,5),collapse="")))
				t[,b:=replicate(20,sample(t$a,6),simplify=F)]
	}
	long[anon,set:=mapply(cull.f,a=k,set,SIMPLIFY=F)]
	t3<-proc.time()
	t3-t2

	setkey(long,k)
	if(F) system.time(save(long,file="/Users/bambrose/long.RData"))
	### combine overlapping sets
	#if(T) {
	library(data.table)
	library(stringdist)
	system.time(load("/Users/bambrose/long.RData"))
	alarm()
	cull<-.1
	cull.f<-function(a,b,cull=.2,noanon=F) {
		if(any(is.na(b))) {b<-na.omit(b);attributes(b)<-NULL}
		if(noanon) {
			if(grepl("\\[ANONYMOUS\\],? ?",a)) stop("k = \"",a,"\"\n",call.=F)
			b<-sub("\\[ANONYMOUS\\],? ?","",b)
		}
		b<-b[nchar(b)>13]
		if(!is.null(cull)) b<-b[stringdist(a=a,b=b,method="jw",p=.1)<cull]
		b
	}
	#}
	grp<-list()
	t3<-proc.time()
	for(i in 1:nrow(long)){
		if(is.na(long[i,set])) {cat("\rSkipping:",i,"\t");next}
		if(nchar(long[i,k])<14) {long[i,set:=NA,nomatch=0];next}
		cat("\rOn:",i,"\b...\t")
		init<-sub("\\[ANONYMOUS\\],? ?","",unique(unlist(long[i,set,nomatch=0])))
		if(!is.null(cull)) init<-cull.f(a=long[i,k],b=init,cull=cull)
		grp[[i]]<-list(unique(unlist(long[init,list(mapply(cull.f,a=k,b=set,cull=cull,noanon=T,SIMPLIFY=F)),nomatch=0])))
		rem<-setdiff(unlist(grp[[i]]),init)
		long[init,set:=NA]
		while(length(rem)){
			init<-unlist(grp[[length(grp)]])
			grp[[length(grp)]][[2]]<-list(unique(unlist(mapply(cull.f,a=rem,b=long[rem,set,nomatch=0],cull=cull,noanon=T,SIMPLIFY=F))))
			grp[[length(grp)]]<-list(unique(unlist(grp[[length(grp)]])))
			long[rem,set:=NA]
			rem<-unique(setdiff(unlist(grp[[i]]),init))
			cat("\ri:",i,"Length rem:",length(rem),"First:",rem[1],"\t")
		}
		grp[[length(grp)]]<-sort(unlist(grp[[length(grp)]]))
		cat("\r\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tCompleted",i,"or",round(i/nrow(long)*100,3),"%\t")
	}
	t4<-proc.time()
	alarm()
	t4-t3
	save(grp,file="grp.RData")
	alarm()
	t4-t0


	### manually audit sets
}

if(F){
	source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
	library(data.table)

	setwd("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/ECON")
	load("ECONwok2db.RData")
	comdb<-wok2db
	setwd("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/SOC")
	load("SOCwok2db.RData")
	comdb<-data.table(rbind(comdb,wok2db))
	rm(wok2db)
	setkey(comdb,field,b.ind.)
	#comdb[,b:=as.character(comdb$b)]
	#setkey(comdb,field,b.ind.)
	#comdb[J("CR"),b:=gsub("(\\w)","\\U\\1",sub(", DOI .+","",comdb$b[comdb$field=="CR"]),perl=T)] #remove DOI & capitalize

	pyj9<-droplevels(data.table(py=comdb["PY",b]$b,j9=comdb["NR",b]$b,key="py"))
	(pyj9<-table(pyj9))
	pynr<-droplevels(data.table(py=comdb["PY",b]$b,nr=as.integer(comdb["NR",b]$b),key="py"))
	pynr<-pynr[,sum(nr),by=py]
	setnames(pynr,old=c("py","V1"),new=c("py","nr"))
	pynr[,cs:=cumsum(pynr$nr)]
	(pynr[,cut:=pynr$cs%/%round(pynr["1970"]$cs/20)])
	cod<-levels(droplevels(comdb["CR",b]$b))
	tu<-length(cod)

	for(i in 6:10*100) assign(paste("dis",i,sep=""),do.call(rbind,lapply(as.list(cod[1:i]),levenshteinSim,cod[1:i])))
	plot(rev(lsos(,n=20)[grep("^dis",dimnames(lsos(n=20))[[1]]),"Size"]))
	size<-rev(lsos(n=20)[grep("^dis",dimnames(lsos(n=20))[[1]]),"Size"])
	#write.table(cbind(size,size2),file="size.txt",sep="\t")
	p<-1:150000
	y<-size

	x<-1:10*100
	x2<-x^2
	sqx<-x^.5

	lm<-lm(y~x-1)
	qm<-lm(y~x+x2-1)
	sm<-lm(y~x+sqx-1)
	summary(lm)
	summary(qm)
	summary(sm)
	anova(lm,qm)
	anova(lm,sm)

	lm$coef
	sm$coef
	qm$coef

	lp<-lm$coef*p/1024/1024
	qp<-((qm$coef[1]*p)+(qm$coef[2]*p^2))/1024/1024
	sp<-((sm$coef[1]*p)+(qm$coef[2]*p^.5))/1024/1024

	max(which(lp<2048))
	max(which(qp<2048)) #15000
	max(which(sp<2048))



	dates<-sort(as.character(unique(comdb[J("PY")]$b)))
	journals<-as.character(unique(comdb[J("J9")]$b))
	index<-list()
	setkey(comdb,field,b)
	system.time(for(i in dates) for(j in journals) {cat("\r",i,j,"\t\t");index[[j]][[i]]<-intersect(comdb[J("PY",i)]$b.ind.,comdb[J("J9",j)]$b.ind.)})
	index<-do.call(cbind,index)
	save(index,file="index.RData")
	setkey(comdb,b.ind.,field)

	# load unicr

	system.time(u<-unicr(index=index))

	u<-data.frame(yu=u[-71,],gain=diff(u$cbeg),loss=diff(u$cend),net=diff(u$cbeg)+diff(u$cend),p=u$yu[-71]+diff(u$cbeg)+diff(u$cend),e=u$yu[-71]+diff(u$cbeg)+diff(u$cend)-u$yu[-1])

	batchs<-12000
	setkey(comdb,b.ind.,field,b)
	j<-1
	rang<-list()
	nun<-list()
	while(!is.na(dimnames(index)[[1]][j])){
		cbeg<-list()
		cat("\n    ")
		for(i in j:dim(index)[1]){
			d<-dimnames(index)[[1]][i]
			cat("\b\b\b\b",d)
			cbeg[[d]]<-length(unique(comdb[J(unlist(index[j:i,]),"CR")]$b))
			if(cbeg[[d]]>batchs) {cbeg[[d]]<-NULL;nun[[length(nun)+1]]<-cbeg[length(cbeg)];break}
		}
		rang[[length(rang)+1]]<-names(cbeg)
		j<-which(dimnames(index)[[1]]==tail(rang[[length(rang)]],1))+1
	}
	if(is.na(dimnames(index)[[1]][j])) nun[[length(nun)+1]]<-cbeg[length(cbeg)]
	cbind(unlist(nun))

	for(i in 1:length(rang)){
		out<-paste(getwd(),.Platform$file.sep,head(rang[[i]],1),ifelse(length(rang[[i]])==1,"",paste("-",tail(rang[[i]],1),sep="")),sep="")
		system(paste("mkdir \"",out,"\"",sep=""))
		wok2bel.f(wok2db=comdb[J(unlist(index[rang[[i]],]))],out=out,man_recode=T)
	}

	#seams: 1916:1919
	out<-paste(getwd(),.Platform$file.sep,1916,"-",1919,sep="")
	system(paste("mkdir \"",out,"\"",sep=""))
	wok2bel.f(wok2db=comdb[J(unlist(index[as.character(1916:1919),]))],out=out,man_recode=T)
	#run it again using
	source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
	wok2bel.f(wok2db=comdb[J(unlist(index[as.character(1916:1919),]))],out=out,man_recode=T,height=20,periodicals=c("MOTOR AGE","LABOUR GAZETTE","RAILWAY AGE GAZETTE","COM FIN CHRON S","OFFICIAL B","EC CIRCULAR","COMMERCE FINANC","NEW YORK J COMM","WOMENS WEAR","MONITEUR","NZ OFFICIAL YB","BANK ARCH","J COMMERCE","MONTHLY LABOR RE","COMMERCIAL AM","NEGLIGENCE COMPENSAT","ANNALIST","ECONOMIST","J DOCUMENTS HOUSE DE","SO LUMBERMAN","AM LUMBERMAN","NATIONS BUSINESS","P GLASS BOTTLE BLOWE","STONE CUTTERS J","IRON MOLDERS J","J POLITICAL EC","EC J","AM J SOCIOLOGY","Q J EC","POLITICAL SCI Q","FORUM","ANN AM ACAD","NILES REGISTER","WOMAN WORKER"))



	#research decision: tables in same statiscal or administrative series are different if from different years
	#research decision: if no year is given, different volumes are treated differently
	#research decision: articles titled part 1 and part 2, etc, can be treated as one
	#research decision: newspapers really shouldn't be treated as the same



}
