
setwd("~/dat/generalismtheory/")
load('160517_generalism_keyObj.RData')
library(ape)
library(phangorn)
library(gdata)
library(picante)

# contains...
# "assoc" which is a cleaned up version of the associations.txt file you sent me (with lines containing NAs and nfs removed)
# "generalism" which is the scores generated from the full associations list you sent
# "paralist" which is the full list of parasites
# "pgd_final" which is the matrix of pairwise distances with all missing data imputed
# "ro_allhost_tree" is the upgma tree of all hosts
# "cdmat" is a "Community data matrix" which just records all the host-parasite associations
# "pd_allhost" is a dataframe with parasite name, Faith's PD, host number

# use tree length as generalism metric
# first need to sort out the indiv parasite's host pgd matrix as it needs a complete triangle
# for(i in 1:length(lowerTriangle(p1mat))){if(is.na(lowerTriangle(p1mat)[i])){lowerTriangle(p1mat)[i]<-upperTriangle(p1mat, byrow=TRUE)[i]}}
# where p1mat equals pgd_final[hlist, hlist]
# can then sum the branch lengths of tree
# sumBLen<-sum(upgma(p1mat)$edge.length) # or sum(nj(p1mat)$edge.length)

#################################################################################################
#
# to calculate Faith's Phylogenetic Diversity (pd)
#
#################################################################################################

# update to use only definitive hosts (bipdef list)
assoc <- bipdef
assoc$h_name <- gsub(" ","_",assoc$new_hname)
assoc$h_name[assoc$new_hname=="Capoetobrama kuschakewitschi kuschakewitsch"] <- "Capoetobrama_kuschakewitschi" #this one didn't match
assoc <- filter(assoc,h_name %in% colnames(pgd_final)) # drop a few more (checked not fixable)
save(assoc,file="~/dat/generalismtheory/assoc_def.RData")

paralist <- as.character(unique(assoc$new_pname))

# make the tree of all hosts (ultrametric using upgma)
allhost_tree<-upgma(pgd_final)
plot(allhost_tree, show.tip.label = FALSE, direction="downwards")
plot(allhost_tree, show.tip.label = FALSE, type="fan")

# make a community data matrix (which is just a binary matrix to record associations - parasites as rows and hosts as columns)
cdmat<-matrix(rep(0,length(paralist)*length(colnames(pgd_final))),nrow=length(paralist))
colnames(cdmat)<-colnames(pgd_final)
rownames(cdmat)<-paralist
for(paras in 1:length(paralist)){
  hlist<-assoc[assoc$new_pname==paralist[paras],3]
  cdmat[paras,intersect(hlist, colnames(pgd_final))]<-1
}
# compute phylogenetic diversity
ro_allhost_tree<-reorder(allhost_tree, "cladewise")
pd_allhost<-pd(cdmat, ro_allhost_tree)
save(pd_allhost,file="~/dat/generalismtheory/pd_allhost.RData")
#################################################################################################
#
# make my dataframe with various indices
#
#################################################################################################

# make a dataframe with nrows=number of parasites i.e. length(paralist) and colnames=c("pname", "hostN", "edges", "sumPGD", "STD", "sSTD")
generalism<-data.frame(pname=character(), hostN=integer(), edges=integer(), sumPGD=numeric(), STD=numeric(), sSTD=numeric(), stringsAsFactors=FALSE)

# compute generalism scores for each parasite based on pairwise genetic distances of all hosts
for(paras in 1:length(paralist)){
	hlist<-assoc[assoc$new_pname==paralist[paras],3]
#  hmat<-pgd_final[hlist, hlist]
#	for(i in 1:length(lowerTriangle(hmat))){if(is.na(lowerTriangle(hmat)[i])){lowerTriangle(hmat)[i]<-upperTriangle(hmat, byrow=TRUE)[i]}}
#	sumBLen_upgma<-sum(upgma(hmat)$edge.length)
  pname<-as.character(paralist[paras])
	hostN<-length(hlist)
	edges<-(hostN*(hostN-1))/2
	sumPGD<-sum(pgd_final[hlist, hlist], na.rm=TRUE)
	STD<-sumPGD/edges
	sSTD<-hostN*STD
	generalism[paras,]<-c(pname, hostN, edges, sumPGD, STD, sSTD)
}

pd_allhost$pname <- rownames(pd_allhost)
generalism2 <- full_join(generalism,pd_allhost,by="pname")
save(generalism2,file="~/dat/generalismtheory/indices29may.RData")
generalism.traits <- left_join(generalism2,par.traits,by=c("pname"="new_pname"))
save(generalism.traits,file="~/dat/generalismtheory/indices29maytraits.RData")

#########################################
## recalculate with individual GEOs
assocG <- bipdefgeo
assocG$h_name <- gsub(" ","_",assocG$new_hname)
assocG$h_name[assocG$new_hname=="Capoetobrama kuschakewitschi kuschakewitsch"] <- "Capoetobrama_kuschakewitschi" #this one didn't match
assocG <- filter(assocG,h_name %in% colnames(pgd_final)) # drop a few more (checked not fixable)
assocG2 <- unite(assocG,paraloc,new_pname,GEO,sep="_",remove=FALSE)
paralist <- unique(assocG2[,1:2])

#################
cdmat<-matrix(rep(0,nrow(paralist)*length(colnames(pgd_final))),nrow=nrow(paralist))
colnames(cdmat)<-colnames(pgd_final)
rownames(cdmat)<-paralist$paraloc
for(paras in 1:nrow(paralist)){
  hlist<-assocG2[assocG2$paraloc==paralist$paraloc[paras],5]
  cdmat[paras,intersect(hlist, colnames(pgd_final))]<-1
}

# compute phylogenetic diversity
ro_allhost_tree<-reorder(allhost_tree, "cladewise")
pd_allhost_GEO<-pd(cdmat, ro_allhost_tree)
save(pd_allhost_GEO,file="~/dat/generalismtheory/pd_allhost_GEO.RData")

###

# make a dataframe with nrows=number of parasites i.e. length(paralist) and colnames=c("pname", "hostN", "edges", "sumPGD", "STD", "sSTD")
generalism.GEO<-data.frame(pnameGEO=character(), hostN=integer(), edges=integer(), sumPGD=numeric(), STD=numeric(), sSTD=numeric(), stringsAsFactors=FALSE)

# compute generalism scores for each parasite based on pairwise genetic distances of all hosts
for(paras in 1:nrow(paralist)){
  hlist<-assocG2[assocG2$paraloc==paralist$paraloc[paras],5]
  #  hmat<-pgd_final[hlist, hlist]
  #	for(i in 1:length(lowerTriangle(hmat))){if(is.na(lowerTriangle(hmat)[i])){lowerTriangle(hmat)[i]<-upperTriangle(hmat, byrow=TRUE)[i]}}
  #	sumBLen_upgma<-sum(upgma(hmat)$edge.length)
  pnameGEO<-as.character(paralist$paraloc[paras])
  hostN<-length(hlist)
  edges<-(hostN*(hostN-1))/2
  sumPGD<-sum(pgd_final[hlist, hlist], na.rm=TRUE)
  STD<-sumPGD/edges
  sSTD<-hostN*STD
  generalism.GEO[paras,]<-c(pnameGEO, hostN, edges, sumPGD, STD, sSTD)
}

pd_allhost_GEO$pnameGEO <- rownames(pd_allhost_GEO)
generalism.GEO.2 <- full_join(generalism.GEO,pd_allhost_GEO,by="pnameGEO")
generalismGEO <- separate(generalism.GEO.2,pnameGEO,c("pname","GEO"),sep = "_")
save(generalismGEO,file="~/dat/generalismtheory/indices29may_geo.RData")
generalism.traits.GEO <- left_join(generalismGEO,par.traits,by=c("pname"="new_pname"))
save(generalism.traits.GEO,file="~/dat/generalismtheory/indices29maytraits_GEO.RData")
