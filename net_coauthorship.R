if(F){
### should do same for authors, but need to do fuzzy set replacement on author names. Treat AF as if it were CR and it should work
source('/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R')
load('/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/wok2db_41.RData')
wok2db_41<-wok2db_41[wok2db_41$fields!='CR',]
rownames(wok2db_41)<-NULL
wok2db_41$fields[wok2db_41$fields=='AF']<-'CR'

##happens to make a co-authorship network!
require(parallel)
afbel<-db2bel.f(wok2db=wok2db_41,out='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out/AFout',man_recode=T,recode_cores=detectCores(),manual_audit=T)
# 3799 Anonymouses cut
#hrm, fuzzy matching fixed! and fast! parallelelelel
#coool, coauthorship network
}