rm(list=ls())
cat('\014')
setwd('/Users/bambrose/Dropbox/GitHub/knowledge-survival')
d<-readLines('dissertation.R')
ds<-readLines('dissertation_source.R')

(gfuns<-grep('^[^ 	].+<- *function',ds,value=F))
funs<-sub('^([^<]+).+$','\\1',ds[gfuns])


gcoms<-grep("#",d)
labels3<-sub("^[^#]*(#.*)$","\\1",d[gcoms])

library(data.table)
dt<-data.table(ds=gfuns,fun=funs)
dt[,'lines':=list(lapply(funs,function(x) grep(x,d)))]
dt[,'ads':=sapply(lines,mean)]
dt<-dt[!is.na(dt$ads)]
lines<-list()
for(i in 1:nrow(dt)) lines[[dt$fun[i]]]<-cbind(dt$ads[i],dt$lines[[i]])
lines<-lines[sapply(lines,length)!=1]
lines<-data.frame(do.call(rbind,lines))
names(lines)<-c('y0','y1')
labels2<-list()
for(i in 1:nrow(dt)) labels2[[dt$fun[i]]]<-paste(dt$fun[i],dt$lines[[i]])
labels2<-do.call(c,labels2)
labels2<-grep(" [0-9]+$",labels2,value=T)
labels1<-dt[,paste(ds,fun)]
ylim<-rev(range(c(lines,gcoms)))
plot.new()
plot.window(xlim=c(0,3),ylim=ylim)
pad1<-max(strwidth(labels1))
pad2<-max(strwidth(c(labels3,labels2)))
cx<-abs(max(strheight(c(labels1,labels2,labels3))))
inpad<-max(strwidth(c(labels3,labels2),units="inches"))
incx<-max(strheight(labels1,units="inches"))
dev.off()
pdf('function_map.pdf',width=7+inpad,height=incx*max(ylim))
plot.new()
axis(2,side=4,at=seq(ylim[1],ylim[2],-100),labels=seq(ylim[1],ylim[2],-100))
plot.window(xlim=c(1-pad1-.15,2+pad2+.15),ylim=ylim)
text(x=1,y=dt$ads,label=labels1,pos=2,col=gray(0,.5),offset=.15) #
text(x=2,y=unlist(dt$lines),label=labels2,pos=4,col=gray(0,.5),offset=.15) #
text(x=2,y=gcoms,label=labels3,pos=4,col="#FF303080",offset=.15) # hashes
segments(
	x0=rep(1,nrow(lines))
	,x1=rep(2,nrow(lines))
	,y0=lines$y0
	,y1=lines$y1
	,col=gray(0,.25))
box()
dev.off()



remove.dups<-function(){
cd<-readLines('/Users/bambrose/Dropbox/2014-2015/Sprints/2/Terminal grep community detection.txt')
cd<-cd[cd!='']

library(stringdist)
system.time(sd<-stringdistmatrix(cd,cd,method='cosine'))

# remove near duplicates
pick.thresh<-function(sd,xlim=NULL){
	quants<-quantile(sd,p=seq(0,1,.05))
	print(quants)
	dsd<-density(sd)
	h<-which(diff(sign(diff(dsd$y)))==-2)+1
	humps<-dsd$x[h]
	t<-which(diff(sign(diff(-dsd$y)))==-2)+1
	troughs<-dsd$x[t]
	#par(mfrow=c(2,1))
	xlim<-if(is.null(xlim)) {range(dsd$x)} else {xlim}
	plot(dsd,xlim=xlim,axes=F)
	axis(2)
	axis(1,at=seq(from=0,to=max(xlim),by=0.1))
	abline(v=humps,col="blue")
	abline(v=troughs,col="green")
	text(y=0,x=humps,col="blue")
	text(y=0,x=troughs,col="green")
	wh<-min(which(dsd$y[h]>=sort(dsd$y[h],decreasing=T)[2]))
	wt<-min(which(dsd$x[t]>dsd$x[wh]))
	choice<-mean(c(humps[wh],troughs[wt]))
	abline(v=choice,col="red")
	h<-hist(sd,breaks=quants,xlim=xlim,axes=F)
	axis(2)
	axis(1,at=seq(from=0,to=max(xlim),by=0.1))
	abline(v=choice,col="red")
	(choice<-max(quants[quants<choice]))
	abline(v=choice,col="purple")
	text(y=max(h$density),x=choice,label=round(choice,3),pos=4,col="purple",cex=.75)
	ret<-list(choice=choice,xlim=range(dsd$x))
	ret
}
graphics.off()
pdf(paste(getwd(),'thresh.pdf',sep=.Platform$file.sep))
choice<-pick.thresh(sd)

system(paste('cd',getwd(),'; mkdir -p nids'))
for(i in 1:nrow(sd)){ #out of all nearly identical, choose 1
	nids<-sd[i,]<=pick.thresh(sd[i,],choice$xlim)$choice
	writeLines(cd[nids],con=paste('nids',.Platform$file.sep,i,'.txt',sep=''))
	nids[sample(which(nids),1)]<-FALSE
	sd<-sd[!nids,!nids]
	cd<-cd[!nids]
	if(i>=nrow(sd)) break
}
dev.off()

library(fastcluster)
library(dendextend)
dd<-as.dist(sd) %>% hclust %>% as.dendrogram

#library(pvclust)
#fit <- pvclust(sd, method.hclust='ward',method.dist='euclidean')
#plot(fit) # dendogram with p values
# add rectangles around groups highly supported by the data
#pvrect(fit, alpha=.95)

library(fpc)
system.time(pc<-pamk(as.dist(sd),5:(nrow(sd)-1)))
pc$nc

library(DendSer)
library(colorspace)
dd %>% rotate_DendSer %>% color_branches(k=pc$nc,col=c("green","blue"),groupLabels=F) %>% plot(leaflab="none",type="rectangle")

sets<-split(cd,pc$pamobject$clustering)
setwd('/Users/bambrose/Dropbox/GitHub/plagiat-scrivener-sync/community detection notes')
for(i in 1:length(sets)) writeLines(sets[[i]],paste(i,'txt',sep='.'))
}


