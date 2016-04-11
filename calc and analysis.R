# Generalism Score Calculation and Analysis
# 10 April 2016
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
load("~/git/generalismtheory/para.ind.RData") 
load("~/git/generalismtheory/host.ind.RData") 

###############################
# Hosts
###############################

# Is maximum host body length associated with generalism?
summary(lm(data=host.ind, betweenness ~ MaxL*network))
bmp("~/dat/generalismtheory/MaxL.bmp")
qplot(data=host.ind,x=MaxL,y=betweenness,color=network) + geom_point() + theme_bw() + geom_smooth(method=loess) + scale_y_log10()
dev.off()

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

