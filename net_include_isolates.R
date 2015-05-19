##### new mel that includes isolates
if(F){
	source('/Users/bambrose/Dropbox/2013-2014/2013-2014_A_Fall/netBYjournalBYyear/dissertation_source.R')
	bel2mel<-bel2mel.f(db2bel=db2bel,out='/Users/bambrose/Dropbox/2013-2014/winter2014_isi_data/1941out')
	mel2net<-mel2net.f(bel2mel)
}
