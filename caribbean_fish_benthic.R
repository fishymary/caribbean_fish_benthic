# -----------------------------------------------------------------------
# FISH FUNCTIONAL GROUPS & BENTHIC COVER IN THE CARIBBEAN
# Mary K. Donovan (donovan.maryk@gmail.com)
# -----------------------------------------------------------------------
rm(list=ls())

# initialization ----------------------------------------------------------
library(dplyr)
library(tidyr)
library(lme4)
library(merTools)

# data --------------------------------------------------------------------
fish <- read.csv('data/GCRMN_fish_4dryad.csv')
benthic <- read.csv('data/GCRMN_benthic_4dryad.csv')

# subset fish species
spp.incl <- read.csv("data/spp_include.csv")
fish.use <- left_join(fish, spp.incl, by="Species")
fish.use$bio_use <- ifelse(fish.use$use4bio==1,fish.use$biomass_g_m2,0)
colnames(fish.use)

# depth hist ---------------------------------------------------------------
png(file='outputs/SOM_depth.png',height=1500,width=1800,res=300)
temp <- subset(fish.use, Depth_m > 0)
hist(temp$Depth_m,main='',xlab='Depth (meters)',col='lightgrey')
dev.off()

# calculate total biomass -------------------------------------------------
tot.bio <- fish.use %>% group_by(DatasetID,Replicate,Location,Year,Habitat) %>% summarise("tot.bio"=sum(bio_use,na.rm=T)) %>% ungroup()
tot.bio.m <- tot.bio
str(tot.bio.m)

tot.bio.m$DatasetID <- as.factor(tot.bio.m$DatasetID)
tot.bio.m$Location <- as.character(tot.bio.m$Location); tot.bio.m$Location <- as.factor(tot.bio.m$Location)
hist(log(tot.bio.m$tot.bio+1)); summary(tot.bio.m$tot.bio)
plot(tot.bio.m$Location,tot.bio.m$tot.bio)

totbio.mod.l <- lmer(log(tot.bio+1) ~ Location + (1|DatasetID), data=tot.bio.m)
summary(totbio.mod.l); plot(residuals(totbio.mod.l)~fitted(totbio.mod.l))
hist(residuals(totbio.mod.l))

newdat <- data.frame(expand.grid(Location=unique(tot.bio.m$Location),DatasetID=11,Country="USVI"))
totbio.mod.pred1 <- predict(totbio.mod.l,newdata=newdat); totbio.mod.pred1 <- exp(totbio.mod.pred1)+1
totbio.mod.pred1 <- data.frame(Location=newdat$Location,tot.bio=totbio.mod.pred1)
totbio.mod.pred1 <- totbio.mod.pred1[with(totbio.mod.pred1, order(tot.bio)),]
totbio.mod.pred <- totbio.mod.pred1
totbio.mod.pred[totbio.mod.pred$tot.bio == min(totbio.mod.pred$tot.bio),]
totbio.mod.pred[totbio.mod.pred$tot.bio == max(totbio.mod.pred$tot.bio),]

totbio.mod.CI <- predictInterval(totbio.mod.l, newdata = newdat, n.sims = 10000, level=0.05)
totbio.mod.CI <- exp(totbio.mod.CI)-1
totbio.mod.CI <- data.frame(Location=newdat$Location,tot.bio=totbio.mod.CI[,1],tot.up=totbio.mod.CI[,2],tot.down=totbio.mod.CI[,3])
totbio.mod.CI[totbio.mod.CI$tot.bio == min(totbio.mod.CI$tot.bio),]
totbio.mod.CI[totbio.mod.CI$tot.bio == max(totbio.mod.CI$tot.bio),]
totbio.mod.CI[totbio.mod.CI$tot.bio == max(totbio.mod.CI$tot.bio),'tot.bio']/totbio.mod.CI[totbio.mod.CI$tot.bio == min(totbio.mod.CI$tot.bio),'tot.bio']

# calculate biomass by trophic level -------------------------------------------------

temp <- fish.use; temp$bio_use[is.na(temp$Trophic)] <- 0
trophic.bio.raw <- temp %>% 
  group_by(DatasetID,Replicate,Location,Year,Trophic,Habitat) %>% 
  summarise("sum"=sum(bio_use,na.rm=T)) %>% 
  spread(key=Trophic,value=sum) %>% 
  ungroup()
trophic.bio.raw.m <- trophic.bio.raw
trophic.bio.raw.m[is.na(trophic.bio.raw.m)] <- 0

trophic.bio.raw.m$DatasetID <- as.factor(trophic.bio.raw.m$DatasetID)
trophic.bio.raw.m$Location <- as.character(trophic.bio.raw.m$Location)

trophic.loc.mod <- data.frame(Location=unique(trophic.bio.raw.m$Location),BR=NA,GD=NA,P=NA,SC=NA,SE=NA)
trophic.loc.modCI <- data.frame(Location=unique(trophic.bio.raw.m$Location))
for(k in c(2:6)){
  temp <- trophic.bio.raw.m
  colnames(temp)[k+4] <- "resp"
  tb.mod.l <- lmer(log(resp+1) ~ Location + (1|DatasetID) , data=temp)
  plot(residuals(tb.mod.l)~fitted(tb.mod.l))
  newdat <- data.frame(expand.grid(Location=unique(temp$Location),DatasetID=11,Country="USVI"))
  # tb.mod.pred <- predict(tb.mod.l,newdata=newdat)
  # tb.mod.pred <- exp(tb.mod.pred) - 1
  # trophic.loc.mod[k] <- tb.mod.pred
  trophic.CI <- predictInterval(tb.mod.l, newdata = newdat, n.sims = 10000, level=0.05)
  trophic.CI <- exp(trophic.CI)-1
  trophic.loc.mod[k] <- trophic.CI$fit
  trophic.loc.modCI <- cbind(trophic.loc.modCI,trophic.CI)
  
}

trophic.loc.mod
colnames(trophic.loc.modCI) <- c("Location","BR.fit","BR.up","BR.down","GD.fit","GD.up","GD.down","P.fit","P.up","P.down","SC.fit","SC.up","SC.down","SE.fit","SE.up","SE.down")
str(trophic.loc.modCI)

trophic.loc.modCI[trophic.loc.modCI$P.fit == max(trophic.loc.modCI$P.fit),'P.fit']/trophic.loc.modCI[trophic.loc.modCI$P.fit == min(trophic.loc.modCI$P.fit),'P.fit']

# barplot -----------------------------------------------------------------
labels.loc <- read.csv("data/loc_labels.csv")
png(file='outputs/Fig2.png',height = 1500,width=3000,res=300)
par(mar=c(7,4,2,1),mgp=c(1.6,.7,0))
temp <- trophic.loc.mod
temp <- left_join(temp, labels.loc)
temp$sum <- rowSums(temp[2:6])
temp <- temp[with(temp, order(sum)),]
col.vec <- c("aquamarine4","aquamarine","lightgreen","#f2a54f","deepskyblue1") #red was "#e0606e"
temp$order <- seq(1:nrow(temp))
b <- barplot(t(temp[c(6,2,3,5,4)]),ylim=c(0,75),xaxt="n",col=col.vec,cex.axis=1.3
             ,ylab=expression("Biomass"~~bgroup("(","g "*m^{-2},")"))
             ,cex.lab=1.2)
text(x=b+.7, y=-3, temp$Label, xpd=TRUE, srt=35, pos=2,cex=0.8)
legend("topleft",legend=c("Predators","Sec. Consumers","Grazers","Browsers","Scrapers"),bty="n",pt.bg=rev(col.vec),pch=22,pt.cex=2)
rows <- temp[c(temp$Location=="Belize_GulfHonduras"|temp$Location=="Belize_SouthernBarrier"|temp$Location=="CaymanIslands_GrandCayman"|temp$Location=="CaymanIslands_LittleandBrac"|temp$Location== "Florida_LowerKeys"|temp$Location== "Florida_UpperKeys"|temp$Location== "Florida_MiddleKeys"|temp$Location=="Mexico_Cozumel_Leeward"),]
points(b[rows$order],temp$sum[rows$order]+2,pch=19,cex=1)
rows <- temp[c(temp$Location=="Bahamas_Remote"|temp$Location=="AandB_Antigua"|temp$Location=="AandB_Barbuda"|temp$Location=="Bahamas_Exuma"|temp$Location=="Bahamas_Other"|temp$Location=="Belize_AtollLeeward"|temp$Location=="Belize_InnerBarrier"|temp$Location=="Cuba_JardinesdelaReina"|temp$Location=="Cuba_Southwest"|temp$Location=="Curacao_NorthWest"|temp$Location=="Guadaloupe_Guadalpupe"|temp$Location=="Jamaic_MontegoBay"|temp$Location=="Jamaic_West"|temp$Location=="PuertoRico_Vieques_Culebra"|temp$Location=="TurksCaicos_TurksCaicos"|temp$Location=="USVI_StJohn"|temp$Location=="USVI_StThomas"),]
points(b[rows$order],temp$sum[rows$order]+2,pch=5,cex=1)
dev.off()

temp <- trophic.loc.modCI
temp <- left_join(temp, labels.loc)
temp$sum <- rowSums(temp[c('BR.fit','GD.fit','P.fit','SC.fit','SE.fit')])
temp <- temp[with(temp, order(sum)),]

write.csv(temp, 'outputs/SOM5.csv',row.names=F)

temp[temp$P.fit == min(temp$P.fit),]
temp[temp$P.fit == max(temp$P.fit),]

temp[temp$SC.fit == min(temp$SC.fit),]
temp[temp$SC.fit == max(temp$SC.fit),]

temp[temp$BR.fit == min(temp$BR.fit),]
temp[temp$BR.fit == max(temp$BR.fit),]

temp[temp$GD.fit == min(temp$GD.fit),]
temp[temp$GD.fit == max(temp$GD.fit),]

temp[temp$SE.fit == min(temp$SE.fit),]
temp[temp$SE.fit == max(temp$SE.fit),]
