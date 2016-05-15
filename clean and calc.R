# final data cleaning and calculations
# 8 May 2016

# use Stephen's cleaned host and parasite name lists

library(dplyr)
library(tidyr)

setwd("~/dat/generalismtheory")
dat <- read.csv("fishparas.csv",header=T)

# update parasite names
pnames <- read.csv("pnames_key.csv")
dat.p <- dat %>% unite(old_pname,P_Genus,P_species,sep=" ", remove=FALSE)
dat.p <- left_join(dat.p,pnames,by="old_pname")

# update host names
hnames <- read.csv("corrected_host_taxonomy.csv")
dat.ph <- dat.p %>% unite(orig_hname, hgenus, H_species,sep=" ",remove=FALSE)
dat.ph <- left_join(dat.ph,hnames,by="orig_hname")

dat.ph.drop <- select(dat.ph,-subspp,-hname,-H_species,-hgenus,-orig_hname,-H_F,-H_O,-P_Genus,-P_species,-old_pname,-P_Family)

# fix trait levels
levels(dat.ph.drop$P_Taxon) <- c("A","A","C","C","CR","M","M","N","N","T","T")
levels(dat.ph.drop$Host_stage) <- c("Both","Both","Definitive","Definitive","Intermediate",NA)
levels(dat.ph.drop$Endoparasite) <- c("Both","Ecto","Endo","Endo",NA)
levels(dat.ph.drop$Gender) <- c("Dioecious","Hermaphrodite","Hermaphrodite",NA)
levels(dat.ph.drop$Complex) <- c(NA,NA,"No","Yes","Yes_No")
levels(dat.ph.drop$Sexual) <- c(NA,"Both",NA,"Sexual","Sexual","Sexual")
levels(dat.ph.drop$Reproduction)[1:2] <- c(NA,NA)
levels(dat.ph.drop$Host_no)[c(1,8)] <- c(NA,NA)
levels(dat.ph.drop$Environment_LHS)[1:2] <- c(NA,NA)
levels(dat.ph.drop$Env_mode) <- c(NA,"napp",NA,"Water","Water")
levels(dat.ph.drop$Motile_LHS) <- c(NA,NA,"No","Yes","Yes")
levels(dat.ph.drop$Vector) <- c(NA,NA,"No","No","Yes")
levels(dat.ph.drop$Trophic) <- c(NA,NA,NA,NA,"No","Yes")
levels(dat.ph.drop$Horizontal) <- c(NA,"Yes",NA,"Yes")
levels(dat.ph.drop$GEO) <- c(NA,"AFR","ANT","AUS","IND",NA,NA,"NEA","NEO","PAL",NA,NA)

save(dat.ph.drop,file="dat.ph.drop.RData")

##
load("dat.ph.drop.RData")

# separate unique host-par pairs and calculate indices of parasites
# degree, S_TD (taxonomy), VarS_TD (taxonomy), betweenness?

# # http://stackoverflow.com/questions/17171148/non-redundant-version-of-expand-grid
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}

# write function to calculate S_TD
# only congeneric hosts will have value 1, so call one species only value 0?

get.S_TD <- function(hname,genus,family,order,class,phylum) {
  dat <- data.frame(cbind(hname,genus,family,order,class,phylum))
  dat <- distinct(dat)
  #print(dat)
  S_TD <- 0
  VarS_TD <- NA
  
  if(nrow(dat) > 1) {
    temp <- data.frame(expand.grid.unique(as.character(dat$hname),as.character(dat$hname)))
    temp$genus <- 0
    temp$family <- 0
    temp$order <- 0
    temp$class <- 0
    temp$phylum <- 0
    
    #print(temp)
    
    for (i in 1:nrow(temp)) {
      temp2 <- rbind(dat[as.character(dat$hname) == temp$X1[i],],dat[as.character(dat$hname) == temp$X2[i],])
      #print(temp2)
      temp$genus[i] <- as.numeric(identical(temp2$genus[1],temp2$genus[2]))
      temp$family[i] <- as.numeric(identical(temp2$family[1],temp2$family[2]))
      temp$order[i] <- as.numeric(identical(temp2$order[1],temp2$order[2]))
      temp$class[i] <- as.numeric(identical(temp2$class[1],temp2$class[2]))
      temp$phylum[i] <- as.numeric(identical(temp2$phylum[1],temp2$phylum[2]))
    }
    
    temp$omega <- with(temp, 6 - genus - family - order - class - phylum)
    
    S_TD <- 2*(sum(temp$omega))/(nrow(dat)*(nrow(dat)-1))
    temp$Var <- (temp$omega - S_TD)^2 
    #print(temp)
    VarS_TD <- 2*(sum(temp$Var))/(nrow(dat)*(nrow(dat)-1))
  }
  return(data.frame(S_TD=S_TD,VarS_TD=VarS_TD,degree=nrow(dat)))
  #return(paste(S_TD,VarS_TD,sep="@"))
}

index <- dat.ph.drop %>% group_by(new_pname) %>% do(get.S_TD(as.character(.$new_hname),as.character(.$genus),as.character(.$hfamily),as.character(.$horder),as.character(.$hclass),as.character(.$hphylum))) %>% ungroup()


save(index,file="~/dat/generalismtheory/index.RData")

load("~/dat/generalismtheory/index.RData")
#ind.ph2 <- left_join(dat.ph.drop,index,by="new_pname")

# # only if definitive...
# ind.def <- ind.ph2 %>% filter(Host_stage=="Definitive" | Host_stage=="Both")

# this takes a long long time to calculate...not feasible (run overnight)
###bip <- dplyr::select(dat.ph.drop,new_pname,new_hname) %>% distinct()
###write.table(bip,file="associations.txt")
# bip$group <- "a"
# library(bipartite)
# pairsweb <- frame2webs(bip,varnames=c("new_hname","new_pname","group"))
# bipindex <- specieslevel(pairsweb$a,level="higher",index="betweenness")
# save(bipindex,file="bipindex.RData")

# summarize traits - mean and sd host traits by parasite
# clean out duplicates

host.traits <- dat.ph.drop %>% select(GEO,MaxL,new_hname) %>% group_by(new_hname,GEO) %>% summarize(maxL=max(MaxL,na.rm=T))

par.traits <- dat.ph.drop %>% select(new_pname,P_Taxon,Endoparasite,Complex,Trophic,Horizontal,Host_no,Host_stage) %>%  filter(!is.na(new_pname))  %>% filter(!(new_pname=="Spiracanthus bovichthys" & P_Taxon=="M")) %>% filter(!(new_pname=="Neoechinorhynchus cylindratus" & P_Taxon=="M")) %>% arrange(new_pname,Endoparasite,Complex,Trophic,Horizontal) %>% distinct(new_pname)

summary(par.traits) # look at complex - host no. how many are complex? intermediate vs definitive? only look at...??
table(par.traits$Host_no,par.traits$Complex,par.traits$Host_stage,useNA="ifany")

par.traits$new_Complex <- par.traits$Complex
par.traits$new_Complex[as.numeric(par.traits$Host_no) > 1 & par.traits$Complex=="No"] <- NA # factor level 1 is "1"
#endo.direct <- par.traits %>% filter(Endoparasite=="Endo")
complexNA <- par.traits[which(is.na(par.traits$Complex)),]

#### HOW SHOULD I USE THIS DISTINCTION? What does NA mean? 
hp <- dat.ph.drop %>% select(new_pname,new_hname,Host_stage,GEO) %>% distinct()
hp$Definitive <- 0
hp$Definitive[hp$Host_stage=="Definitive" | hp$Host_stage=="Both"] <- 1
hp$Intermediate <- 0
hp$Intermediate[hp$Host_stage=="Intermediate" | hp$Host_stage=="Both"] <- 1

# make new index only for definitive hosts
dat.ph.2 <- dat.ph.drop
dat.ph.2$Definitive <- 0
dat.ph.2$Definitive[dat.ph.2$Host_stage=="Definitive" | dat.ph.2$Host_stage=="Both"] <- 1
dat.ph.2$Intermediate <- 0
dat.ph.2$Intermediate[dat.ph.2$Host_stage=="Intermediate" | dat.ph.2$Host_stage=="Both"] <- 1

with(dat.ph.2,table(Complex,Host_stage,useNA="ifany"))

dat.ph.def <- dat.ph.2 %>% filter(Definitive==1)

index.complex <- dat.ph.def %>% group_by(new_pname) %>% do(get.S_TD(as.character(.$new_hname),as.character(.$genus),as.character(.$hfamily),as.character(.$horder),as.character(.$hclass),as.character(.$hphylum))) %>% ungroup()

index.geo <- dat.ph.def %>% group_by(new_pname,GEO) %>% do(get.S_TD(as.character(.$new_hname),as.character(.$genus),as.character(.$hfamily),as.character(.$horder),as.character(.$hclass),as.character(.$hphylum))) %>% ungroup()

hp.htraits <- left_join(hp,host.traits,by=c("new_hname","GEO"))
hp.traits <- left_join(hp.htraits,par.traits,by=c("new_pname"))

# hp.traits.dup <- hp.traits[(duplicated(hp.traits[,c(1,2)]) | duplicated(hp.traits[,c(1,2)],fromLast = T)),] %>% arrange(new_hname,new_pname) %>% filter(!is.na(GEO))

hp.psum.geo <- hp.traits %>% group_by(new_pname,GEO,Definitive,Intermediate,P_Taxon,Endoparasite,Complex,Trophic,Horizontal) %>% summarize(meanMaxL=mean(maxL,na.rm=T),sdMaxL=sd(maxL,na.rm=T)) # THIS IS NOT DIVIDED BY INTERMEDIATE/DEFINITIVE

hp.geo <- hp.traits %>% group_by(new_pname,GEO,P_Taxon,Endoparasite,Complex,Trophic,Horizontal) %>% summarize(meanMaxL=mean(maxL,na.rm=T),sdMaxL=sd(maxL,na.rm=T),maxMaxL=max(maxL,na.rm=T))

hp.traits2 <- hp.traits %>% select(-GEO) %>% distinct()
hp.psum.nogeo <- hp.traits2 %>% group_by(new_pname,P_Taxon,Endoparasite,Complex,Trophic,Horizontal) %>% summarize(meanMaxL=mean(maxL,na.rm=T),sdMaxL=sd(maxL,na.rm=T),maxMaxL=max(maxL,na.rm=T))


#
hp.full.simple <- left_join(hp.psum.nogeo,index,by="new_pname")
hp.full.geo.stage <-  left_join(hp.psum.geo,index,by="new_pname")
#hp.full.geo <- left_join(hp.geo,index,by="new_pname")
hp.def.geo <- left_join(index.geo, hp.geo, by=c("new_pname","GEO"))
save(hp.def.geo, file="hp.def.geo.RData")

hp.host <- left_join(hp.traits,index,by="new_pname")
hp.host.sum <- hp.host %>% group_by(new_hname) %>% summarize(maxL=max(maxL,na.rm=T),meandeg=mean(degree),meddeg=median(degree),meanSTD=mean(S_TD),meanVarS_TD=mean(VarS_TD,na.rm=T))

hp.def <- left_join(index.complex,hp.psum.nogeo,by="new_pname")
save(hp.def,file="hp.def.RData")

save(hp.full.simple,file="hp.full.simple.RData")
save(hp.full.geo.stage,file="hp.full.geo.stage.RData")
save(hp.host.sum,file="hp.host.sum.RData")
save(hp.full.geo,file="hp.full.geo.RData")


# dat.ph.drop2 <- dat.ph.drop %>% select(new_pname,P_Taxon,Endoparasite,Complex,Trophic,Horizontal,GEO,MaxL)  %>% filter(!is.na(new_pname))  %>% filter(!(new_pname=="Spiracanthus bovichthys" & P_Taxon=="M")) %>% filter(!(new_pname=="Neoechinorhynchus cylindratus" & P_Taxon=="M")) %>% arrange(new_pname,Endoparasite,Complex,Trophic,Horizontal) %>% distinct(new_pname,GEO,MaxL)
# 
# traits.noGEO <- dat.ph.drop2 %>% group_by(new_pname,Endoparasite,Complex,Trophic,Horizontal) %>% summarize(meanMaxL=mean(MaxL,na.rm=T),sdMaxL=sd(MaxL,na.rm=T))
# 
# traits.GEO <- dat.ph.drop2 %>% group_by(new_pname,Endoparasite,Complex,Trophic,Horizontal,GEO) %>% summarize(meanMaxL=mean(MaxL,na.rm=T),sdMaxL=sd(MaxL,na.rm=T))

#traits.noGEO[duplicated(traits.noGEO$new_pname) | duplicated(traits.noGEO$new_pname,fromLast = T),] %>% arrange(new_pname) %>% View()




# pairs.t1 <- dat.ph.drop %>% select(new_hname,new_pname,Host_stage) %>% distinct()
# pairs.dup <- pairs.t1[duplicated(pairs.t1[,1:2]) | duplicated(pairs.t1[,1:2], fromLast=TRUE),] %>% arrange(new_hname,new_pname) # all doubles
# 
# cleanHS <- function(d) {
#   if(any(is.na(d$Host_stage))) {
#     new <- as.character(d$Host_stage[!is.na(d$Host_stage)])
#     return(new)
#   }
#   if(any(d$Host_stage=="Both")) {return("Both")}
#   if(any(d$Host_stage=="Intermediate") & any(d$Host_stage=="Definitive")) {return("Both")}
# }
# 
# pairs.new <- pairs.dup %>% group_by(new_hname,new_pname) %>% do(newHS=cleanHS(.))
# pairs.new$newHS <- as.character(pairs.new$newHS)
# pairs.t2 <- left_join(pairs.t1,pairs.new,by=c("new_hname","new_pname"))
# pairs.t2$Host_stage[!is.na(pairs.t2$newHS)] <- pairs.t2$newHS[!is.na(pairs.t2$newHS)]
# pairs.t <- distinct(pairs.t2)
# pairs.t$newHS <- NULL
# dat.ph.drop2 <- right_join(dat.ph.drop,pairs.t,by=c("new_hname","new_pname"))
# dat.ph.drop2 <- dat.ph.drop2 %>% select(-Host_stage.x)
# dat.ph.drop3 <- unique(dat.ph.drop2)

# 
