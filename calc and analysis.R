# Generalism Score Calculation and Analysis
# 8 May 2016
###############################

# packages
library(dplyr)
library(tidyr)
library(ggplot2)


###############################
# Functions
###############################


###############################
# Read in clean data
###############################
load("~/dat/generalismtheory/hp.full.simple.RData")
load("~/dat/generalismtheory/hp.full.geo.stage.RData")
load("~/dat/generalismtheory/hp.host.sum.RData")

###############################
# Not divided by Geo & Stage
###############################

# Is maximum host body length associated with generalism?
m1 <- lm(data=hp.full.simple, log(degree) ~ meanMaxL)
summary(m1)
plot(m1)
pdf("~/dat/generalismtheory/MaxL-degree.pdf")
qplot(data=hp.full.simple,x=meanMaxL,y=degree) + geom_jitter() + theme_bw() 
dev.off()

with(hp.full.simple,cor.test(meanMaxL,degree,method="spearman"))

hp.full.simple$gen <- FALSE
hp.full.simple$gen[hp.full.simple$degree > 1] <- TRUE

m2 <- lm(data=hp.full.simple, gen ~ meanMaxL)
summary(m2)

pdf("~/dat/generalismtheory/MaxL-gen.pdf")
qplot(data=hp.full.simple,y=meanMaxL,x=gen) + geom_violin() + theme_bw() 
dev.off()

## try zip model? not really biologically relevant

# look at endos only
hp.endo <- hp.full.simple %>% filter(Endoparasite =="Endo")
qplot(data=hp.endo,y=meanMaxL,x=gen) + geom_violin() + theme_bw() 
bmp("endo-deg.bmp")
qplot(data=hp.endo,x=meanMaxL,y=degree) + geom_jitter() + theme_bw()
dev.off()
with(hp.endo,cor.test(meanMaxL,degree,method="spearman"))
with(hp.endo,cor.test(meanMaxL,as.numeric(gen),method="spearman"))
with(hp.endo,wilcox.test(meanMaxL~gen))
hp.endo %>% group_by(gen) %>% summarize(mean=mean(meanMaxL,na.rm=T),median=median(meanMaxL,na.rm=T))
# mean is higher for non-generalists 
#summary(lm(data=hp.endo,degree ~ meanMaxL))

bmp("endo-std.bmp")
qplot(data=hp.endo,x=meanMaxL,y=S_TD) + geom_point() + theme_bw() 
dev.off()
with(hp.endo,cor.test(meanMaxL,S_TD,method="spearman"))
summary(lm(data=hp.endo, S_TD ~ meanMaxL))
plot(lm(data=hp.endo, S_TD ~ meanMaxL))

# what about variance - much stronger correlation - only sd when there is more than one host...
bmp("endo-deg-var.bmp")
qplot(data=hp.endo,x=sdMaxL,y=degree) + geom_jitter() + theme_bw() 
dev.off()
with(hp.endo,cor.test(sdMaxL,degree,method="spearman"))
bmp("endo-std-var.bmp")
qplot(data=hp.endo,x=sdMaxL,y=S_TD) + geom_point() + theme_bw() 
dev.off()
with(hp.endo,cor.test(sdMaxL,S_TD,method="spearman"))


# host perspective on body size/generalism relationship
hphost2 <- hp.host.sum %>% gather("index","value",3:6)
bmp("maxL-host.bmp")
qplot(data=hphost2,x=maxL,y=value) + geom_jitter() + theme_bw() + facet_wrap(ncol=2,~index,scales="free")
dev.off()
# qplot(data=hp.host.sum,x=maxL,y=meandeg) + geom_jitter() + theme_bw() 
# qplot(data=hp.host.sum,x=maxL,y=meddeg) + geom_jitter() + theme_bw()
# qplot(data=hp.host.sum,x=maxL,y=meanSTD) + geom_jitter() + theme_bw()
# qplot(data=hp.host.sum,x=maxL,y=meanVarS_TD) + geom_jitter() + theme_bw() 

dev.off()

with(hp.host.sum,cor.test(maxL,meandeg,method="spearman"))
with(hp.host.sum,cor.test(maxL,meddeg,method="spearman"))
with(hp.host.sum,cor.test(maxL,meanSTD,method="spearman"))
with(hp.host.sum,cor.test(maxL,meanVarS_TD,method="spearman"))

######################################################

# Is host growth rate associated with generalism?
summary(lm(data=host.ind, betweenness ~ as.numeric(K)*network))
bmp("~/dat/generalismtheory/K.bmp")
qplot(data=host.ind,x=as.numeric(K),y=betweenness,color=network) + geom_point() + theme_bw() + geom_smooth(method=loess) + scale_y_log10()
dev.off()

# Is host life span associated with generalism?
summary(lm(data=host.ind, betweenness ~ as.numeric(Y)*network))
bmp("~/dat/generalismtheory/Y.bmp")
qplot(data=host.ind,x=as.numeric(Y),y=betweenness,color=network) + geom_point() + theme_bw() + geom_smooth(method=loess) + scale_y_log10()
dev.off()

# Is host age at maturity (Ym) associated with generalism?
summary(lm(data=host.ind, betweenness ~ as.numeric(Ym)*network))
bmp("~/dat/generalismtheory/Ym.bmp")
qplot(data=host.ind,x=as.numeric(Ym),y=betweenness,color=network) + geom_point() + theme_bw() + geom_smooth(method=loess) + scale_y_log10()
dev.off()

# Is host trophic level (T) associated with generalism?
summary(lm(data=host.ind, betweenness ~ as.numeric(T)*network))
bmp("~/dat/generalismtheory/T.bmp")
qplot(data=host.ind,x=as.numeric(T),y=betweenness,color=network) + geom_point() + theme_bw() + geom_smooth(method=loess) + scale_y_log10()
dev.off()

###############################
# Parasites
###############################
para.ind.long <- para.ind %>% gather(key="index",value="value",S_TD,VarS_TD,degree)
para.ind.long2 <- para.ind.long
para.ind.long2$P_Taxon <- "Total"
para.ind.long <- bind_rows(para.ind.long,para.ind.long2)

# Is there a difference in generalism scores between ecto and endo parasites?
png("~/dat/generalismtheory/endo.png")
qplot(data=para.ind.long,x=Endoparasite,y=value,colour=P_Taxon) + facet_grid(P_Taxon ~ index) + geom_boxplot() + scale_y_log10() + theme_bw() 
dev.off()

# Is there a difference in generalism scores between hermaphroditic and dioecious parasites?
png("~/dat/generalismtheory/gender.png")
qplot(data=para.ind.long,x=Gender,y=value,colour=P_Taxon) + facet_grid(P_Taxon ~ index) + geom_boxplot() + scale_y_log10() + theme_bw() 
dev.off()

# Is there a difference in generalism scores between oviparous/viviparous/ovoviviparous parasites?
png("~/dat/generalismtheory/repro.png")
qplot(data=para.ind.long,x=Reproduction,y=value,colour=P_Taxon) + facet_grid(P_Taxon ~ index) + geom_boxplot() + scale_y_log10() + theme_bw() 
dev.off()

# Is generalism correlated with the number of hosts required for the parasite to complete its life cycle (complex)?
png("~/dat/generalismtheory/host_no.png")
qplot(data=para.ind.long,x=Host_no,y=value,colour=P_Taxon) + facet_grid(P_Taxon ~ index) + geom_boxplot() + scale_y_log10() + theme_bw() 
dev.off()

# Is generalism correlated with environmental vs direct transmission?
png("~/dat/generalismtheory/env.png")
qplot(data=para.ind.long,x=Environment_LHS,y=value,colour=P_Taxon) + facet_grid(P_Taxon ~ index) + geom_boxplot() + scale_y_log10() + theme_bw() 
dev.off()

# Is generalism correlated with trophic transmission?
png("~/dat/generalismtheory/trophic.png")
qplot(data=para.ind.long,x=Trophic,y=value,colour=P_Taxon) + facet_grid(P_Taxon ~ index) + geom_boxplot() + scale_y_log10() + theme_bw() 
dev.off()


para.ind.longT <- para.ind %>% gather(key="TraitType",value="TraitValue",Endoparasite:Horizontal)
sumPara <- para.ind.longT %>% select(-pname) %>% group_by(TraitType,TraitValue,P_Taxon) %>% summarize_each(funs(mean(.,na.rm=T),q.025=quantile(.,probs=(0.025),na.rm=T),q.975=quantile(.,probs=(0.975),na.rm=T),n=length,sd(.,na.rm=T)))
sumParaTot  <- para.ind.longT %>% select(-pname,-P_Taxon) %>% group_by(TraitType,TraitValue) %>% summarize_each(funs(mean(.,na.rm=T),q.025=quantile(.,probs=(0.025),na.rm=T),q.975=quantile(.,probs=(0.975),na.rm=T),n=length,sd(.,na.rm=T))) %>% mutate(P_Taxon="Total") %>% bind_rows(sumPara)

write.csv(sumParaTot,"~/git/generalismtheory/summary_parasitetraits.csv")


