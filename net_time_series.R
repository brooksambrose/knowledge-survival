source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
setkey(comdb,b.ind.,fields,b)
for(i in 1){
	out<-paste(getwd(),.Platform$file.sep,head(rang[[i]],1),ifelse(length(rang[[i]])==1,"",paste("-",tail(rang[[i]],1),sep="")),sep="")
	wok2bel<-wok2bel.f(wok2db=comdb[J(unlist(index[rang[[i]],]))],out=out,man_recode=T,saved_recode=sets)
}

source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
jbymel<-list()
for(j in colnames(index)) for(i in as.character(1900:1917)) {cat(i,j,"\n");jbymel[[j]][i]<-list(bel2mel.f(wok2bel=wok2bel,subset=index[[i,j]],man_recode=T,out=getwd()))}
(jbymel<-do.call(cbind,jbymel))
save(jbymel,file="jbymel.RData")

jb3ymel<-list()
for(j in colnames(index)) for(i in as.character(1902:1917)) {cat(i,j,"\n");jb3ymel[[j]][i]<-list(bel2mel.f(wok2bel=wok2bel,subset=unlist(index[(which(rownames(index)==i)-2):which(rownames(index)==i),j]),man_recode=T,out=getwd()))}
(jb3ymel<-do.call(cbind,jb3ymel))
save(jb3ymel,file="jb3ymel.RData")


source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
jb3ynet<-list()
for(j in colnames(index)) for(i in as.character(1902:1917)) {cat(i,j,"\n");jb3ynet[[j]][i]<-list(mel2net(jb3ymel[[i,j]]))}
jb3ynet<-do.call(cbind,jb3ynet)
jb3ynet[!sapply(jb3ynet,length)]<-NA
save(jb3ynet,file="jb3ynet.RData")

source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
load(file="/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/jb3ynet.RData")
pois<-list()
for(j in colnames(jb3ynet)) for(i in rownames(jb3ynet)) if(is.na(jb3ynet[i,j])) {pois[[j]][[i]]<-NA} else {pois[[j]][[i]]<-list(thatgirlis(jb3ynet[i,j][[1]]$tpzcrel))}
pois<-do.call(cbind,pois)
pois<-pois[,!apply(is.na(pois),2,all)]
save(pois,file="pois.RData")

source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
pdf("journal_poisson.pdf")
for(j in colnames(pois)) for(k in c(0:2,">=3")) {cat(j,k,"");try(plotpois(pois,year=1902:1917,jour=j,count=k,q1="25%",q2="50%",q3="75%"))}
dev.off()

source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R")
sbel<-wok2bel.f(wok2db=soc,out="/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/SOC",man_recode=T)
ebel<-wok2bel.f(wok2db=econ,out="/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/ECON",man_recode=T)

load("wok2bel.RData")

sel<-list()
sel[["ASR"]][["1936"]]<-bel2mel.f(wok2bel=wok2bel,out="/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear",trim_pendants=T,man_recode=T)

eel<-list()
debug(bel2mel.f)
eel[["AM ECON REV"]][["1936"]]<-bel2mel.f(wok2bel=ebel,out="/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/ECON",trim_pendants=T,man_recode=T)


library(statnet)
library(latentnet)
load("bel2mel.RData")

el<-list()
el[["ASR"]][["1936"]]<-bel2mel
n<-list()
n$ASR$`1936`<-network(el$ASR$`1936`$tpzcrel[,1:2],matrix.type="edgelist",directed=F)
n$ASR$`1936`%e%"ew"<-el$ASR$`1936`[[1]]$ew
n[["AM ECON REV"]][["1936"]]<-network(eel$"AM ECON REV"$`1936`$tpzcrel[,1:2],matrix.type="edgelist",directed=F)
n[["AM ECON REV"]][["1936"]]%e%"ew"<-eel$"AM ECON REV"$`1936`[[1]]$ew


f<-ergm(n$ASR$`1936`~edges,verbose=1,response="ew",reference=~Poisson,burnin=500000)
mcmc.diagnostics(f)

t<-triad.census(n$ASR$`1936`, mode ="graph")

source("/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R", chdir = TRUE)
thatgirlis(n=n[["AM ECON REV"]][["1936"]])

#bel by year
if(F){
	db2bel_41_by<-list()
	years<-sort(as.character(unique(wok2db_41$b[wok2db_41$fields=='PY'])))
	for(i in years){
		set<-NULL
		set<-as.character(wok2db_41$id[wok2db_41$fields=='PY'&wok2db_41$b==i])
	db2bel_41_by[[i]]<-db2bel.f(wok2db=droplevels(wok2db_41[wok2db_41$id.%in%set,])
		,out='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/by_year'
		,man_recode=F
		,saved_recode=NULL
		,cut_samp_def=0
	)
}
save(db2bel_41_by,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/by_year/db2bel_41_by.RData')
}

#mel by year
if(F){
	bel2mel_41_by<-list()
	for(i in years) {bel2mel_41_by[[i]]<-bel2mel.f(db2bel=db2bel_41_by[[i]]
		,out='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/by_year'
		)
	}
	save(bel2mel_41_by,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/by_year/bel2mel_41_by.RData')

}

#net by year
if(F){
	mel2net_41_by<-list()
	for(i in years){
	mel2net_41_by[[i]]<-mel2net.f(bel2mel_41_by[[i]])
	}
	save(mel2net_41_by,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/by_year/mel2net_41_by.RData')

}

#pois by year
if(F){
	thatgirlis_41_by<-list()
	for(i in years){
		if((!length(mel2net_41_by[[i]]))|class(mel2net_41_by[[i]]$tpcrel)=='logical') {thatgirlis_41_by[[i]]<-NULL;next}
		thatgirlis_41_by[[i]]<-thatgirlis.f(mel2net_41_by[[i]]$tpcrel)
	}
	save(thatgirlis_41_by,file='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/by_year/thatgirlis_41_by.RData')
}