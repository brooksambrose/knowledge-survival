cat('\014')
source('/Users/bambrose/Dropbox/GitHub/knowledge-survival/dissertation_source.R')
library(data.table)
library(igraph)

system.time(cfnet<-cfinder2all.f(cf.in='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/CFinder1941-exact',which='communities'))

net<-cfnet[[1]]

strat<-rbindlist(list(net[,list(src,ew)],net[,list(tgt,ew)]),use.names=F)
strat<-strat[,list(mxk=max(ew),mnk=min(ew),k=list(sort(ew))),by=src]
strat[,`:=`(r=mxk-mnk,avk=round(sapply(k,FUN=mean),2),mdk=sapply(k,FUN=quantile,p=.5))]
strat[,amdif:=avk-mdk]
setkey(strat,amdif)
k<-strat[,.N,keyby=k]
k[,k:=k+2] # k=3 counted as 1, etc.
k[,Nk:=round(N/k,2)]

### visualize global graph

setkey(net,src,tgt,ew)
graph<-graph.edgelist(as.matrix(net[,list(src,tgt)]),directed=F)
E(graph)$weight<-net[,ew]
community<-fastgreedy.community(graph)
plot(community,graph,vertex.size=1,edge.lty='blank',vertex.label=NA)

### test idea of nestedness

# here we will try to measure each node according to whether the majority of its sending ties are to adjascent strata, conditional on size of strata. We need a predicted and expected vector.

et<-strat[,
	list(o=list(table(factor(do.call(c,k),levels=et$mxk))))
	,keyby=mxk]
et[,Nt:=strat[,list(Nt=length(do.call(c,k))),keyby=mxk]$Nt]
et[,`:=`(e=lapply(1:length(k$k),function(x) prop.table(k$N[1:x])))]

#if that stratum sent as many ties as it did each to a random alter at or below its level, then the expectation of ties to each substratum would be the conditional stratum proportions. Because the counts depend on the distribution, it is simpler to run a permutation test rather than calculate conditional probabilities.

perm.test<-function(x) {
	reps=100000
	r<-replicate(
		n=reps
		,expr=table(
			sample(
				x=factor(as.integer(names(x)))
				,size=Nt # something is broken here
				,replace=T
				,prob=prop.table(x) # these need to be conditional probabilities of choosing a node, but here they are conditional probabilities of choosing a tie proportion, which should just reproduce the original data, all fubar
			)
		)
	)
	sd<-apply(
		X=r
		,MARGIN=1
		,FUN=sd
	)
	m<-apply(
		X=r
		,MARGIN=1
		,FUN=mean
	)
	pego<-apply(
		X=r>matrix(rep(x,reps),ncol=reps)
		,MARGIN=1
		,FUN=mean
	)
	width<-function(y) y #c(y,rep(NA,nrow(et)-length(y)))

	ret<-list(
		o=width(x)
		,e=width(m)
		,sde=width(sd)
		,t=width((x-m)/sd)
		,pego=width(pego)
	)
	ret
}
dat<-lapply(et$o[-1],perm.test) # this needs to be an mapply, and we need to feed different population proportions and tie proportions
dimn<-list(ks=et$mxk,kr=et$mxk,dat=c('o','e','sde','t','pego'))
et2excel<-array(dim=sapply(dimn,length),dimnames=dimn)
for(i in 2:nrow(et2excel)) for(j in dimn$dat) et2excel[i,,j]<-dat[[i-1]][[j]]

graphics.off()
plot.new()
plot.window(
	xlim=c(1,37)
	,ylim=c(-1,1)*log(abs(range(et2excel[,,'t'][is.finite(et2excel[,,'t'])]))))
abline(h=-log(3),col='blue')
abline(h=log(3),col='blue')
abline(h=0)
axis(1)
axis(2)
box()

for(i in 2:37){
w<-(et2excel[i,,'t'])<0
w[is.na(w)]<-F
l<-log(abs(et2excel[i,,'t']))
l[w]<-l[w]*-1
lines(l,type='o',col=heat.colors(36)[i],ylog=T,bg='darkgray')
}
dev.off()

pdf('s-r-kclique.pdf')
conf<-3

t<-et2excel[,,'t']
w<-which(t<0)
t<-log(abs(t))
t[w]<--1*t[w]

bx<-function(){
boxplot(t)
abline(h=-log(conf),col='blue')
abline(h=log(conf),col='blue')
abline(h=0)
}
colr<-rainbow(nrow(t))
for(i in 2:nrow(t)) {
	bx()
	lines(t[i,],type='o',col=colr[i])
	text(x=0,y=t[i,which.min(is.na(t[i,]))],labels=i,col=colr[i],cex=.5)
}


bx()
int.case<-apply(t,MARGIN=1,function(x) mean(abs(na.omit(x))>log(conf)))
ic<-list()
qic<-quantile(int.case,na.rm=T,prob=seq(0,1,.25))
colr<-rainbow(length(qic))
names(colr)<-names(qic)
for(i in names(qic)) {
	ic[[i]]<-which.min(abs(int.case-qic[i]))
	lines(t[ic[[i]],],type='o',col=colr[i])
	text(x=0,y=t[ic[[i]],which.min(is.na(t[ic[[i]],]))],labels=paste(names(ic[i]),names(ic[[i]])),col=colr[i],cex=.5)
}


flop<-t
flop[]<-NA
for(i in 2:nrow(flop)) flop[i,(ncol(flop)-i+1):ncol(flop)]<-t[i,1:i]
colnames(flop)<-(ncol(flop)-1):0
flop<-flop[,ncol(flop):1]

bx<-function(){
	boxplot(flop)
	abline(h=-log(conf),col='blue')
	abline(h=log(conf),col='blue')
	abline(h=0)
}

colr<-rainbow(nrow(t))
for(i in 2:nrow(t)) {
	bx()
	lines(flop[i,],type='o',col=colr[i])
	text(x=0,y=flop[i,which.min(is.na(flop[i,]))],labels=i,col=colr[i],cex=.5)
}

dev.off()



