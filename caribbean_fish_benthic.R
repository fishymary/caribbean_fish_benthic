# -----------------------------------------------------------------------
# FISH FUNCTIONAL GROUPS & BENTHIC COVER IN THE CARIBBEAN
# -----------------------------------------------------------------------
rm(list=ls())
Start <- Sys.time()

# initialization ----------------------------------------------------------
library(dplyr) # v.0.7.6
library(tidyr) # v.0.8.0
library(lme4) # v.1.1.13
library(merTools) # v.0.3.0
library(mgcv) # v.1.8-23
library(MuMIn) # v.1.15.6
library(car) # v.2.1.5
library(vegan) # v.2.4-4

# data --------------------------------------------------------------------
fish <- read.csv('data/GCRMN_fish_4dryad.csv')
benthic <- read.csv('data/GCRMN_benthic_4dryad.csv')

# subset fish species
spp.incl <- read.csv("data/spp_include.csv")
fish.use <- left_join(fish, spp.incl, by="Species")
fish.use$bio_use <- ifelse(fish.use$use4bio==1,fish.use$biomass_g_m2,0)
colnames(fish.use)

# fix a loc name
fish.use$Location[fish.use$Location=="Bahamas_Exuma"] <- "Bahamas_Other"

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
labels.loc <- data.frame(
  Location = as.factor(c("Jamaic_NorthCentral", "Mexico_SouthEastYucatan",
                         "Jamaic_MontegoBay", "Jamaic_West",
                         "Belize_GulfHonduras", "PuertoRico_Turrumote",
                         "Belize_SouthernBarrier", "Mexico_NorthEastYucatan",
                         "DominicanRepublic_North", "Belize_AtollLeeward",
                         "DominicanRepublic_South", "PuertoRico_LaPaguera",
                         "Mexico_ChinchorroBank", "USVI_StJohn", "PuertoRico_JobosBay",
                         "StVincentGrenadines_Grenadines",
                         "Belize_CentralBarrier", "Colombia_SanAndreas", "Panama_SanBlas",
                         "Guadaloupe_Guadalpupe", "TurksCaicos_TurksCaicos",
                         "StBarthelemy_StBarthelemy",
                         "StKittsNevis_StKittsNevis", "Cuba_North", "Belize_AtollWindward",
                         "Cuba_Southwest", "Honduras_BayIslands",
                         "Belize_InnerBarrier", "PuertoRico_Vieques_Culebra", "USVI_StThomas",
                         "BVI_BVI", "CaymanIslands_LittleandBrac",
                         "CaymanIslands_GrandCayman", "Martinique_Martinique",
                         "Curacao_NorthWest", "Colombia_Providencia",
                         "Bahamas_Other", "AandB_Antigua", "StEustatius_StEustatius",
                         "Bahamas_CaySalBank", "Honduras_NearShore",
                         "Cuba_JardinesdelaReina", "Mexico_Cozumel_Leeward",
                         "Bahamas_Remote", "Bahamas_Nassau", "Bahamas_Andros",
                         "AandB_Barbuda", "Florida_LowerKeys", "Bahamas_Exuma",
                         "Florida_UpperKeys", "Florida_MiddleKeys")),
  Label = as.factor(c("Jamica NC", "Mexico SE Yucatan", "Jamaica MB",
                      "Jamaica W", "Belize Gulf Honduras",
                      "Puerto Rico Turrumote", "Belize S Barrier",
                      "Mexico NE Yucatan", "DR North", "Belize Atoll Leeward", "DR South",
                      "Puerto Rico La Paguera", "Mexico Chinchorro",
                      "USVI St John", "Puerto Rico Jobos Bay", "Grenadines",
                      "Belize C Barrier", "Colombia San Andres",
                      "Panama San Blas", "Guadeloupe", "Turks & Caicos",
                      "St Barthelemy", "St Kitts & Nevis", "Cuba North",
                      "Belize Atoll Windward", "Cuba Southeast",
                      "Honduras Bay Islands", "Belize Inner Barrier",
                      "Puerto Rico Vieques", "USVI St. Thomas", "BVI", "Little Cayman",
                      "Grand Cayman", "Martinique", "Curacao Northwest",
                      "Colombia Providencia", "Bahamas other", "Antigua",
                      "St. Eustatius", "Cay Sal Bank", "Honduras Nearshore",
                      "Cuba Jardines", "Mexico Cozumel",
                      "Bahamas Remote", "Bahamas Nassau", "Bahamas Andros", "Barbuda",
                      "Florida Lower Keys", "Bahamas Exuma",
                      "Florida Upper Keys", "Florida Middle Keys"))
)

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
temp$P.fit[temp$P.fit == max(temp$P.fit)]/temp$P.fit[temp$P.fit == min(temp$P.fit)]

temp[temp$SC.fit == min(temp$SC.fit),]
temp[temp$SC.fit == max(temp$SC.fit),]

temp[temp$BR.fit == min(temp$BR.fit),]
temp[temp$BR.fit == max(temp$BR.fit),]

temp[temp$GD.fit == min(temp$GD.fit),]
temp[temp$GD.fit == max(temp$GD.fit),]

temp[temp$SE.fit == min(temp$SE.fit),]
temp[temp$SE.fit == max(temp$SE.fit),]

# correlations among fish -------------------------------------------------
panel.hist <- function(x, ...)
{ usr <- par("usr"); on.exit(par(usr))
par(usr = c(usr[1:2], 0, 1.5) )
h <- hist(x, plot = FALSE)
breaks <- h$breaks; nB <- length(breaks)
y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}
panel.blank <- function(x, y)
{ }
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex <- 0.6/strwidth(txt)
  #text(0.5, 0.5, txt, cex = cex * r)
  text(0.5, 0.5, txt, cex = cex)
}

pairs(trophic.loc.mod[c(2:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)

temp <- left_join(totbio.mod.CI,trophic.loc.modCI,by="Location")
head(temp)

png(file=file.path(getwd(),'outputs',"SOM8.png"),height=3000,width=3200,res=300)
par(mfrow=c(5,5),mar=c(1,1,1,1),oma=c(4,4,0,0),mgp=c(1.6,.7,0),xpd=T)

x <- temp$tot.bio; y <- temp$P.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,cex.axis=1.4)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

x <- temp$tot.bio; y <- temp$SC.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,cex.axis=1.4)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$P.fit; y <- temp$SC.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

x <- temp$tot.bio; y <- temp$GD.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,cex.axis=1.4)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$P.fit; y <- temp$GD.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$SC.fit; y <- temp$GD.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")
plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

x <- temp$tot.bio; y <- temp$BR.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,cex.axis=1.4)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$P.fit; y <- temp$BR.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$SC.fit; y <- temp$BR.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$GD.fit; y <- temp$BR.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

plot(c(1,1),c(1,1),type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="")

x <- temp$tot.bio; y <- temp$SE.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,cex.axis=1.4)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$P.fit; y <- temp$SE.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$SC.fit; y <- temp$SE.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$GD.fit; y <- temp$SE.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$BR.fit; y <- temp$SE.fit
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,labels=NA)
if((cor.test(x,y))$p.value < 0.05){
  text(0.1*max(x),0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

mtext("Total",outer=T,side=1,at=0.1,line=1.3,cex=1.5)
mtext("Predators",outer=T,side=1,at=0.3,line=1.3,cex=1.5)
mtext("Sec. Consumers",outer=T,side=1,at=0.5,line=1.3,cex=1.5)
mtext("Grazers",outer=T,side=1,at=0.7,line=1.3,cex=1.5)
mtext("Browsers",outer=T,side=1,at=0.9,line=1.3,cex=1.5)
# mtext(expression("Biomass"~~bgroup("(","g"*m^{-2},")")),outer=T,side=1,line=3.5,cex=1.5)

mtext("Predators",outer=T,side=2,at=0.9,line=1.3,cex=1.5)
mtext("Sec. Consumers",outer=T,side=2,at=0.7,line=1.3,cex=1.5)
mtext("Grazers",outer=T,side=2,at=0.5,line=1.3,cex=1.5)
mtext("Browsers",outer=T,side=2,at=0.3,line=1.3,cex=1.5)
mtext("Scrapers",outer=T,side=2,at=0.1,line=1.3,cex=1.5)
# mtext(expression("Biomass"~~bgroup("(","g"*m^{-2},")")),outer=T,side=2,line=3.5,cex=1.5)

dev.off()

# summarise benthic cover -------------------------------------------------

# modeled
coral.m <- benthic
coral.m <- coral.m[!is.na(coral.m$TotalCoral),]
temp <- coral.m %>% group_by(Location) %>% summarise(nrow=length(TotalCoral)) %>% ungroup(); temp <- subset(temp, nrow > 5)
coral.m <- coral.m[coral.m$Location %in% temp$Location,]
coral.m$Location <- as.character(coral.m$Location); coral.m$Location <- as.factor(coral.m$Location)
coral.m$DatasetID <- as.factor(coral.m$DatasetID); coral.m$Year <- as.factor(coral.m$Year)
str(coral.m)

coral.mod <- lmer(sqrt(sqrt(TotalCoral)) ~ Location + (1|DatasetID) , data=coral.m)
summary(coral.mod); plot(residuals(coral.mod)~fitted(coral.mod))
newdat <- data.frame(expand.grid(Location=levels(coral.m$Location),DatasetID=16))
coral.mod.pred <- predict(coral.mod,newdata=newdat); coral.mod.pred <- (coral.mod.pred)^4
coral.mod.pred <- data.frame(Location=newdat$Location,coral=coral.mod.pred)
coral.mod.pred <- coral.mod.pred[with(coral.mod.pred, order(coral)),]
coral.mod.pred
coral.mod.CI <- predictInterval(coral.mod, newdata = newdat, n.sims = 10000, level=0.05)
coral.mod.CI <- (coral.mod.CI)^4
coral.mod.CI <- data.frame(Location=newdat$Location,coral=coral.mod.CI[,1],coral.up=coral.mod.CI[,2],coral.down=coral.mod.CI[,3])
coral.mod.CI[coral.mod.CI$coral == min(coral.mod.CI$coral),]
coral.mod.CI[coral.mod.CI$coral == max(coral.mod.CI$coral),]


macro.m <- benthic
macro.m <- macro.m[!is.na(macro.m$TotalMacro),]
macro.m$DatasetID <- as.factor(macro.m$DatasetID); macro.m$Year <- as.factor(macro.m$Year)
macro.m$Location <- as.character(macro.m$Location); macro.m$Location <- as.factor(macro.m$Location)
macro.mod <- lmer(sqrt(sqrt(TotalMacro)) ~ Location + (1|DatasetID) , data=macro.m)
summary(macro.mod); plot(residuals(macro.mod)~fitted(macro.mod))
newdat <- data.frame(expand.grid(Location=unique(macro.m$Location),DatasetID=10))
macro.mod.pred <- predict(macro.mod,newdata=newdat); macro.mod.pred <- (macro.mod.pred)^4
macro.mod.pred <- data.frame(Location=newdat$Location,macro=macro.mod.pred)
macro.mod.pred <- macro.mod.pred[with(macro.mod.pred, order(macro)),]
macro.mod.pred
macro.mod.CI <- predictInterval(macro.mod, newdata = newdat, n.sims = 10000, level=0.05)
macro.mod.CI <- (macro.mod.CI)^4
macro.mod.CI <- data.frame(Location=newdat$Location,macro=macro.mod.CI[,1],macro.up=macro.mod.CI[,2],macro.down=macro.mod.CI[,3])
macro.mod.CI[macro.mod.CI$macro == min(macro.mod.CI$macro),]
macro.mod.CI[macro.mod.CI$macro == max(macro.mod.CI$macro),]

temp <- left_join(totbio.mod.CI, coral.mod.CI, by="Location")
temp <- left_join(temp, macro.mod.CI, by="Location")
temp <- left_join(temp, labels.loc)

write.csv(temp, file.path(getwd(),'outputs',"SOM4.csv"),row.names=F)

# total bio v benthic -----------------------------------------------------
totbio.pred.ben <- left_join(totbio.mod.CI, coral.mod.CI, by="Location")
totbio.pred.ben <- left_join(totbio.pred.ben, macro.mod.CI, by="Location")
plot(totbio.pred.ben$tot.bio,totbio.pred.ben$coral)
plot(totbio.pred.ben$tot.bio,totbio.pred.ben$macro)
totbio.pred.ben$cm <- log(totbio.pred.ben$coral/totbio.pred.ben$macro)

# join trophic and benthic ------------------------------------------------
trophic.benthic.m <- left_join(trophic.loc.mod, coral.mod.CI, by="Location")
trophic.benthic.m <- left_join(trophic.benthic.m, macro.mod.CI, by="Location")
str(trophic.benthic.m)
pairs(trophic.benthic.m[c(2:7,10)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
trophic.benthic.m$cm <- log(trophic.benthic.m$coral/trophic.benthic.m$macro)

# panel fig - trophic biomass -------------------------------------------------------------

plot.fish <- function(dat, x, y, i, metric){ #i is what row of mod.out to fill in
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(ncvTest(tp)$p <= 0.05){print('Non-constant Error')}

  if(metric=='biomass'){
    ylim.go <- c(0,60); if(c(i > 6 & i < 13)){ylim.go <- c(0,80)}; if(i > 12){ylim.go <- c(-2.2,2.2)}
    xlim.go <- c(0,85); if(i==2|i==8|i==14){xlim.go <- c(0,25)}; if(i==3|i==9|i==15){xlim.go <- c(0,36)}
    if(i==4|i==10|i==16){xlim.go <- c(0,8)}; if(i==5|i==11|i==17){xlim.go <- c(0,8.5)}; if(i==6|i==12|i==18){xlim.go <- c(0,20)}
  }
  if(metric=='size'){
    ylim.go <- c(0,60); if(c(i > 6 & i < 13)){ylim.go <- c(0,80)}; if(i > 12){ylim.go <- c(-2.2,2.2)}
    xlim.go <- c(5,30)
  }

  xaxt.go <- rep('n',18); if(i > 12){xaxt.go <- 's'}
  yaxt.go <- rep('n',18); if(i==1|i==7|i==13){yaxt.go <- 's'}

  if(m.sel$class[1]=='lm' & summary(tp)$coefficients[8] <= 0.05){
    tp.pred <- predict(tp,se.fit=T)
    tp.pred <- data.frame(pred= tp.pred$fit,se=tp.pred$se.fit, x=temp.sub[,x])
    tp.pred <- tp.pred[with(tp.pred, order(x)),]
    tp.pred$CIup <- tp.pred$pred + (1.96*tp.pred$se); tp.pred$CIdown <- tp.pred$pred - (1.96*tp.pred$se)

    plot(temp.sub[,y]~temp.sub[,x],pch=21,col="black",bg="grey",bty="l",xlab="",ylab=""
         ,cex.axis=1.2,cex.lab=1.3,ylim=ylim.go,xlim=xlim.go,xaxt=xaxt.go,yaxt=yaxt.go,type="n")
    # points(tp.pred$pred~tp.pred$x,type="l",lwd=2,lty=2)
    axis(1,labels=NA)
    axis(2,labels=NA)
    polygon(c(tp.pred$x,rev(tp.pred$x)),c(tp.pred$CIup,rev(tp.pred$CIdown))
              ,col=rgb(140,140,140,100,max=255),border=NA)
    points(tp.pred$pred~tp.pred$x,type="l",lwd=2)
    points(temp[,y]~temp[,x],pch=21,col="black",bg="grey")
    points(temp[cook >= 0.50,y]~temp[cook >= 0.50,x],pch=1,col="red",cex=2)
    } else {
    if(summary(tpg)$s.table[4] <= 0.05){
      tpg.pred <- predict(tpg, se.fit=T)
      tpg.pred <- data.frame(x=temp.sub[,x],y=tpg.pred$fit,up=tpg.pred$fit+1.96*tpg.pred$se.fit
                             ,down=tpg.pred$fit-1.96*tpg.pred$se.fit)
      tpg.pred <- tpg.pred[with(tpg.pred,order(x)),]
      plot(temp.sub[,y]~temp.sub[,x],pch=21,col="black",bg="grey",bty="l",xlab="",ylab=""
           ,cex.axis=1.2,cex.lab=1.3,ylim=ylim.go,xlim=xlim.go,xaxt=xaxt.go,yaxt=yaxt.go,type="n")
      axis(2,labels=NA)
      axis(1,labels=NA)
      polygon(c(tpg.pred$x,rev(tpg.pred$x)),c(tpg.pred$up,rev(tpg.pred$down)),col=rgb(169,169,169,150,max=255),border=NA)
      points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2)
      points(temp[,y]~temp[,x],pch=21,col="black",bg="grey")
      points(temp[cook >= 0.50,y]~temp[cook >= 0.50,x],pch=1,col="red",cex=2)
    } else {
      plot(temp.sub[,y]~temp.sub[,x],pch=21,col="black",bg="grey",bty="l",xlab="",ylab=""
           ,cex.axis=1.2,cex.lab=1.3,ylim=ylim.go,xlim=xlim.go,xaxt=xaxt.go,yaxt=yaxt.go,type="n")
      axis(2,labels=NA)
      axis(1,labels=NA)
      points(temp[,y]~temp[,x],pch=21,col="black",bg="grey")
    }
  }


}

png(file.path(getwd(),'outputs',"Fig3.png"),height=2000,width=3600,res=300)
par(mfrow=c(3,6),mar=c(1,1,1,1),oma=c(4,3,0,0),xpd=F,mgp=c(1.6,.7,0))
plot.fish(totbio.pred.ben,'tot.bio','coral',1,'biomass')
tp <- c('P','SC','GD','BR','SE')
for(z in 1:5){
  plot.fish(trophic.benthic.m,x=tp[z],y='coral',i=z+1,'biomass')
}
plot.fish(totbio.pred.ben,'tot.bio','macro',7,'biomass')
for(z in 1:5){
  plot.fish(trophic.benthic.m,x=tp[z],y='macro',i=z+7,'biomass')
}
plot.fish(totbio.pred.ben,'tot.bio','cm',13,'biomass')
for(z in 1:5){
  plot.fish(trophic.benthic.m,x=tp[z],y='cm',i=z+13,'biomass')
}
mtext("Coral Cover (%)", outer=T, side=2, at=0.83,line=1)
mtext("Macroalgal Cover (%)", outer=T, side=2, at=0.5,line=1)
mtext("log(Coral/Macrogalae)", outer=T, side=2, at=0.175,line=1)
mtext("Total",outer=T,side=1,at=0.09,line=1)
mtext("Predators",outer=T,side=1,at=0.248,line=1)
mtext("Sec. Consumers",outer=T,side=1,at=0.416,line=1)
mtext("Grazers",outer=T,side=1,at=0.584,line=1)
mtext("Browsers",outer=T,side=1,at=0.752,line=1)
mtext("Scrapers",outer=T,side=1,at=0.92,line=1)
mtext(expression("Biomass"~~bgroup("(","g"*m^{-2},")")),outer=T,side=1,line=2.8)

dev.off()

mod.out <- data.frame(
  y = as.character(c("coral", "coral", "coral", "coral", "coral",
                     "coral", "macro", "macro", "macro", "macro",
                     "macro", "macro", "cm", "cm", "cm", "cm", "cm", "cm")),
  x = as.character(c("total.bio", "P", "SC", "GD", "BR", "SE",
                     "total.bio", "P", "SC", "GD", "BR", "SE",
                     "total.bio", "P", "SC", "GD", "BR", "SE")),
  int = rep(NA,18),
  int.se = rep(NA,18),
  int.t = rep(NA,18),
  int.p = rep(NA,18),
  slp = rep(NA,18),
  slp.se = rep(NA,18),
  slp.t = rep(NA,18),
  slp.p = rep(NA,18),
  f = rep(NA,18),
  df = rep(NA,18),
  r2 = rep(NA,18), stringsAsFactors=F
)

for(i in c(1,7,13)){
  y <- mod.out[i,'y']
  x <- 'tot.bio'
  dat <- totbio.pred.ben
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(m.sel$class[1]=='lm'){
    mod.out[i,3] <- summary(tp)$coefficients[1]; mod.out[i,7] <- summary(tp)$coefficients[2]
    mod.out[i,4] <- summary(tp)$coefficients[3]; mod.out[i,8] <- summary(tp)$coefficients[4]
    mod.out[i,5] <- summary(tp)$coefficients[5]; mod.out[i,9] <- summary(tp)$coefficients[6]
    mod.out[i,6] <- summary(tp)$coefficients[7]; mod.out[i,10] <- summary(tp)$coefficients[8]
    mod.out[i,11] <- summary(tp)$fstatistic[1]
    mod.out[i,12] <- summary(tp)$df[2]
    mod.out[i,13] <- summary(tp)$r.squared
  }
  if(m.sel$class[1]=='gam'){
    mod.out[i,3] <- summary(tpg)$p.table[1]; mod.out[i,4] <- summary(tpg)$p.table[2]
    mod.out[i,5] <- summary(tpg)$p.table[3]; mod.out[i,6] <- summary(tpg)$p.table[4]
    mod.out[i,7] <- summary(tpg)$s.table[1]; mod.out[i,8] <- summary(tpg)$s.table[2]
    mod.out[i,10] <- summary(tpg)$s.table[4]; mod.out[i,11] <- summary(tpg)$s.table[3]
    mod.out[i,12] <- summary(tpg)$n
    mod.out[i,13] <- summary(tpg)$dev.expl
  }
}

for(i in c(2:6,8:12,14:18)){
  y <- mod.out[i,'y']
  x <- mod.out[i,'x']
  dat <- trophic.benthic.m
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(m.sel$class[1]=='lm'){
    mod.out[i,3] <- summary(tp)$coefficients[1]; mod.out[i,7] <- summary(tp)$coefficients[2]
    mod.out[i,4] <- summary(tp)$coefficients[3]; mod.out[i,8] <- summary(tp)$coefficients[4]
    mod.out[i,5] <- summary(tp)$coefficients[5]; mod.out[i,9] <- summary(tp)$coefficients[6]
    mod.out[i,6] <- summary(tp)$coefficients[7]; mod.out[i,10] <- summary(tp)$coefficients[8]
    mod.out[i,11] <- summary(tp)$fstatistic[1]
    mod.out[i,12] <- summary(tp)$df[2]
    mod.out[i,13] <- summary(tp)$r.squared
  }
  if(m.sel$class[1]=='gam'){
    mod.out[i,3] <- summary(tpg)$p.table[1]; mod.out[i,4] <- summary(tpg)$p.table[2]
    mod.out[i,5] <- summary(tpg)$p.table[3]; mod.out[i,6] <- summary(tpg)$p.table[4]
    mod.out[i,7] <- summary(tpg)$s.table[1]; mod.out[i,8] <- summary(tpg)$s.table[2]
    mod.out[i,10] <- summary(tpg)$s.table[4]; mod.out[i,11] <- summary(tpg)$s.table[3]
    mod.out[i,12] <- summary(tpg)$n
    mod.out[i,13] <- summary(tpg)$dev.expl
  }
}

write.csv(mod.out, 'outputs/Table1.csv',row.names=F)

# mean size  ------------------------------
### need to replicate the lengths based on number column
size.bio <- fish.use %>% group_by(DatasetID,Replicate,Location,Year,Length) %>% summarise("num"=sum(Number,na.rm=T)) %>% ungroup()
size.bio <- size.bio[!is.na(size.bio$Length),]
size.bio$num[size.bio$num < 1 & size.bio$num > 0] <- 1
size.bio <- subset(size.bio, Length > 0)
size.out <- data.frame(DatasetID=NA,Location=NA,size=NA)
for(i in 1:nrow(size.bio)){
  t <- rep(size.bio$Length[i],size.bio$num[i])
  t <- data.frame(DatasetID=size.bio$DatasetID[i],Location=size.bio$Location[i],size=t)
  size.out <- rbind(size.out,t)
}
str(size.out)
size.out <- size.out[2:nrow(size.out),]

size.out$DatasetID <- as.character(size.out$DatasetID); size.out$DatasetID <- as.factor(size.out$DatasetID)
size.out$Location <- as.character(size.out$Location); size.out$Location <- as.factor(size.out$Location)

hist(size.out$size)
hist(log(size.out$size+1))

size.bio.mod <- lmer(log(size+1) ~ Location + (1|DatasetID), data=size.out)
summary(size.bio.mod); plot(residuals(size.bio.mod)~fitted(size.bio.mod))
hist(residuals(size.bio.mod))
newdat <- data.frame(expand.grid(Location=unique(size.out$Location),DatasetID=48,Country="USVI"))
size.bio.mod.pred <- predict(size.bio.mod,newdata=newdat); size.bio.mod.pred <- exp(size.bio.mod.pred)+1
size.bio.mod.pred <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.pred)
size.bio.mod.pred <- size.bio.mod.pred[with(size.bio.mod.pred, order(mean.size)),]
size.bio.mod.pred <- size.bio.mod.pred
size.bio.mod.CI <- predictInterval(size.bio.mod, newdata = newdat, n.sims = 10000, level=0.05)
size.bio.mod.CI <- exp(size.bio.mod.CI)+1
size.bio.mod.CI <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.CI[,1],mean.size.up=size.bio.mod.CI[,2],mean.size.down=size.bio.mod.CI[,3])

size.bio.ben <- left_join(size.bio.mod.CI, coral.mod.pred, by="Location")
size.bio.ben <- left_join(size.bio.ben, macro.mod.pred, by="Location")
pairs(size.bio.ben[c(2,5:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)

size.bio.ben$cm <- log(size.bio.ben$coral/size.bio.ben$macro)


### by trophic level
size.bioT <- fish.use %>% group_by(DatasetID,Replicate,Location,Year,Length,Trophic) %>% summarise("num"=sum(Number,na.rm=T)) %>% ungroup()
size.bioT <- size.bioT[!is.na(size.bioT$Length),]
size.bioT$num[size.bioT$num < 1 & size.bioT$num > 0] <- 1
size.bioT <- subset(size.bioT, Length > 0)
str(size.bioT)


########### PREDATORS
size.bioP <- subset(size.bioT, Trophic=="P")
size.bioP.out <- data.frame(DatasetID=NA,Location=NA,size=NA)
for(i in 1:nrow(size.bioP)){
  t <- rep(size.bioP$Length[i],size.bioP$num[i])
  t <- data.frame(DatasetID=size.bioP$DatasetID[i],Location=size.bioP$Location[i],size=t)
  size.bioP.out <- rbind(size.bioP.out,t)
}
str(size.bioP.out)
size.bioP.out <- size.bioP.out[2:nrow(size.bioP.out),]

size.bioP.out$DatasetID <- as.character(size.bioP.out$DatasetID)
size.bioP.out$DatasetID <- as.factor(size.bioP.out$DatasetID)
size.bioP.out$Location <- as.character(size.bioP.out$Location)
size.bioP.out$Location <- as.factor(size.bioP.out$Location)

size.bio.mod.P <- lmer(log(size+1) ~ Location + (1|DatasetID), data=size.bioP.out)
summary(size.bio.mod.P)
plot(residuals(size.bio.mod.P)~fitted(size.bio.mod.P))
hist(residuals(size.bio.mod.P))
newdat <- data.frame(expand.grid(Location=unique(size.bioP.out$Location),DatasetID=48,Country="USVI"))
size.bio.mod.P.pred <- predict(size.bio.mod.P,newdata=newdat); size.bio.mod.P.pred <- exp(size.bio.mod.P.pred)+1
size.bio.mod.P.pred <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.P.pred)
size.bio.mod.P.pred <- size.bio.mod.P.pred[with(size.bio.mod.P.pred, order(mean.size)),]
size.bio.mod.P.pred <- size.bio.mod.P.pred
size.bio.P.mod.CI <- predictInterval(size.bio.mod.P, newdata = newdat, n.sims = 10000, level=0.05)
size.bio.P.mod.CI <- exp(size.bio.P.mod.CI)+1
size.bio.P.mod.CI <- data.frame(Location=newdat$Location,Pmean.size=size.bio.P.mod.CI[,1],Pmean.size.up=size.bio.P.mod.CI[,2],Pmean.size.down=size.bio.P.mod.CI[,3])

sizeP.bio.ben <- left_join(size.bio.P.mod.CI, coral.mod.pred, by="Location")
sizeP.bio.ben <- left_join(sizeP.bio.ben, macro.mod.pred, by="Location")
pairs(sizeP.bio.ben[c(2,5:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
sizeP.bio.ben$cm <- log(sizeP.bio.ben$coral/sizeP.bio.ben$macro)

########### Sec Con
size.bioSC <- subset(size.bioT, Trophic=="SC")
size.bioSC.out <- data.frame(DatasetID=NA,Location=NA,size=NA)
for(i in 1:nrow(size.bioSC)){
  t <- rep(size.bioSC$Length[i],size.bioSC$num[i])
  t <- data.frame(DatasetID=size.bioSC$DatasetID[i],Location=size.bioSC$Location[i],size=t)
  size.bioSC.out <- rbind(size.bioSC.out,t)
}
str(size.bioSC.out)
size.bioSC.out <- size.bioSC.out[2:nrow(size.bioSC.out),]

size.bioSC.out$DatasetID <- as.character(size.bioSC.out$DatasetID)
size.bioSC.out$DatasetID <- as.factor(size.bioSC.out$DatasetID)
size.bioSC.out$Location <- as.character(size.bioSC.out$Location)
size.bioSC.out$Location <- as.factor(size.bioSC.out$Location)

size.bio.mod.SC <- lmer(log(size+1) ~ Location + (1|DatasetID), data=size.bioSC.out)
summary(size.bio.mod.SC)
plot(residuals(size.bio.mod.SC)~fitted(size.bio.mod.SC))
hist(residuals(size.bio.mod.SC))
newdat <- data.frame(expand.grid(Location=unique(size.bioSC.out$Location),DatasetID=48,Country="USVI"))
size.bio.mod.SC.pred <- predict(size.bio.mod.SC,newdata=newdat); size.bio.mod.SC.pred <- exp(size.bio.mod.SC.pred)+1
size.bio.mod.SC.pred <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.SC.pred)
size.bio.mod.SC.pred <- size.bio.mod.SC.pred[with(size.bio.mod.SC.pred, order(mean.size)),]
size.bio.mod.SC.pred <- size.bio.mod.SC.pred
size.bio.SC.mod.CI <- predictInterval(size.bio.mod.SC, newdata = newdat, n.sims = 10000, level=0.05)
size.bio.SC.mod.CI <- exp(size.bio.SC.mod.CI)+1
size.bio.SC.mod.CI <- data.frame(Location=newdat$Location,SCmean.size=size.bio.SC.mod.CI[,1],SCmean.size.up=size.bio.SC.mod.CI[,2],SCmean.size.down=size.bio.SC.mod.CI[,3])

sizeSC.bio.ben <- left_join(size.bio.SC.mod.CI, coral.mod.pred, by="Location")
sizeSC.bio.ben <- left_join(sizeSC.bio.ben, macro.mod.pred, by="Location")
pairs(sizeSC.bio.ben[c(2,5:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
sizeSC.bio.ben$cm <- log(sizeSC.bio.ben$coral/sizeSC.bio.ben$macro)

########### Grazers
size.bioGD <- subset(size.bioT, Trophic=="GD")
size.bioGD.out <- data.frame(DatasetID=NA,Location=NA,size=NA)
for(i in 1:nrow(size.bioGD)){
  t <- rep(size.bioGD$Length[i],size.bioGD$num[i])
  t <- data.frame(DatasetID=size.bioGD$DatasetID[i],Location=size.bioGD$Location[i],size=t)
  size.bioGD.out <- rbind(size.bioGD.out,t)
}
str(size.bioGD.out)
size.bioGD.out <- size.bioGD.out[2:nrow(size.bioGD.out),]

size.bioGD.out$DatasetID <- as.character(size.bioGD.out$DatasetID)
size.bioGD.out$DatasetID <- as.factor(size.bioGD.out$DatasetID)
size.bioGD.out$Location <- as.character(size.bioGD.out$Location)
size.bioGD.out$Location <- as.factor(size.bioGD.out$Location)

size.bio.mod.GD <- lmer(log(size+1) ~ Location + (1|DatasetID), data=size.bioGD.out)
summary(size.bio.mod.GD)
plot(residuals(size.bio.mod.GD)~fitted(size.bio.mod.GD))
hist(residuals(size.bio.mod.GD))
newdat <- data.frame(expand.grid(Location=unique(size.bioGD.out$Location),DatasetID=48,Country="USVI"))
size.bio.mod.GD.pred <- predict(size.bio.mod.GD,newdata=newdat); size.bio.mod.GD.pred <- exp(size.bio.mod.GD.pred)+1
size.bio.mod.GD.pred <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.GD.pred)
size.bio.mod.GD.pred <- size.bio.mod.GD.pred[with(size.bio.mod.GD.pred, order(mean.size)),]
size.bio.mod.GD.pred <- size.bio.mod.GD.pred
size.bio.GD.mod.CI <- predictInterval(size.bio.mod.GD, newdata = newdat, n.sims = 10000, level=0.05)
size.bio.GD.mod.CI <- exp(size.bio.GD.mod.CI)+1
size.bio.GD.mod.CI <- data.frame(Location=newdat$Location,GDmean.size=size.bio.GD.mod.CI[,1],GDmean.size.up=size.bio.GD.mod.CI[,2],GDmean.size.down=size.bio.GD.mod.CI[,3])

sizeGD.bio.ben <- left_join(size.bio.GD.mod.CI, coral.mod.pred, by="Location")
sizeGD.bio.ben <- left_join(sizeGD.bio.ben, macro.mod.pred, by="Location")
pairs(sizeGD.bio.ben[c(2,5:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
sizeGD.bio.ben$cm <- log(sizeGD.bio.ben$coral/sizeGD.bio.ben$macro)

########### Browsers
size.bioBR <- subset(size.bioT, Trophic=="BR")
size.bioBR.out <- data.frame(DatasetID=NA,Location=NA,size=NA)
for(i in 1:nrow(size.bioBR)){
  t <- rep(size.bioBR$Length[i],size.bioBR$num[i])
  t <- data.frame(DatasetID=size.bioBR$DatasetID[i],Location=size.bioBR$Location[i],size=t)
  size.bioBR.out <- rbind(size.bioBR.out,t)
}
str(size.bioBR.out)
size.bioBR.out <- size.bioBR.out[2:nrow(size.bioBR.out),]

size.bioBR.out$DatasetID <- as.character(size.bioBR.out$DatasetID)
size.bioBR.out$DatasetID <- as.factor(size.bioBR.out$DatasetID)
size.bioBR.out$Location <- as.character(size.bioBR.out$Location)
size.bioBR.out$Location <- as.factor(size.bioBR.out$Location)

size.bio.mod.BR <- lmer(log(size+1) ~ Location + (1|DatasetID), data=size.bioBR.out)
summary(size.bio.mod.BR)
plot(residuals(size.bio.mod.BR)~fitted(size.bio.mod.BR))
hist(residuals(size.bio.mod.BR))
newdat <- data.frame(expand.grid(Location=unique(size.bioBR.out$Location),DatasetID=48,Country="USVI"))
size.bio.mod.BR.pred <- predict(size.bio.mod.BR,newdata=newdat); size.bio.mod.BR.pred <- exp(size.bio.mod.BR.pred)+1
size.bio.mod.BR.pred <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.BR.pred)
size.bio.mod.BR.pred <- size.bio.mod.BR.pred[with(size.bio.mod.BR.pred, order(mean.size)),]
size.bio.mod.BR.pred <- size.bio.mod.BR.pred
size.bio.BR.mod.CI <- predictInterval(size.bio.mod.BR, newdata = newdat, n.sims = 10000, level=0.05)
size.bio.BR.mod.CI <- exp(size.bio.BR.mod.CI)+1
size.bio.BR.mod.CI <- data.frame(Location=newdat$Location,BRmean.size=size.bio.BR.mod.CI[,1],BRmean.size.up=size.bio.BR.mod.CI[,2],BRmean.size.down=size.bio.BR.mod.CI[,3])

sizeBR.bio.ben <- left_join(size.bio.BR.mod.CI, coral.mod.pred, by="Location")
sizeBR.bio.ben <- left_join(sizeBR.bio.ben, macro.mod.pred, by="Location")
pairs(sizeBR.bio.ben[c(2,5:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
sizeBR.bio.ben$cm <- log(sizeBR.bio.ben$coral/sizeBR.bio.ben$macro)

########### Scrapers
size.bioSE <- subset(size.bioT, Trophic=="SE")
size.bioSE.out <- data.frame(DatasetID=NA,Location=NA,size=NA)
for(i in 1:nrow(size.bioSE)){
  t <- rep(size.bioSE$Length[i],size.bioSE$num[i])
  t <- data.frame(DatasetID=size.bioSE$DatasetID[i],Location=size.bioSE$Location[i],size=t)
  size.bioSE.out <- rbind(size.bioSE.out,t)
}
str(size.bioSE.out)
size.bioSE.out <- size.bioSE.out[2:nrow(size.bioSE.out),]

size.bioSE.out$DatasetID <- as.character(size.bioSE.out$DatasetID)
size.bioSE.out$DatasetID <- as.factor(size.bioSE.out$DatasetID)
size.bioSE.out$Location <- as.character(size.bioSE.out$Location)
size.bioSE.out$Location <- as.factor(size.bioSE.out$Location)

size.bio.mod.SE <- lmer(log(size+1) ~ Location + (1|DatasetID), data=size.bioSE.out)
summary(size.bio.mod.SE)
plot(residuals(size.bio.mod.SE)~fitted(size.bio.mod.SE))
hist(residuals(size.bio.mod.SE))
newdat <- data.frame(expand.grid(Location=unique(size.bioSE.out$Location),DatasetID=48,Country="USVI"))
size.bio.mod.SE.pred <- predict(size.bio.mod.SE,newdata=newdat); size.bio.mod.SE.pred <- exp(size.bio.mod.SE.pred)+1
size.bio.mod.SE.pred <- data.frame(Location=newdat$Location,mean.size=size.bio.mod.SE.pred)
size.bio.mod.SE.pred <- size.bio.mod.SE.pred[with(size.bio.mod.SE.pred, order(mean.size)),]
size.bio.mod.SE.pred <- size.bio.mod.SE.pred
size.bio.SE.mod.CI <- predictInterval(size.bio.mod.SE, newdata = newdat, n.sims = 10000, level=0.05)
size.bio.SE.mod.CI <- exp(size.bio.SE.mod.CI)+1
size.bio.SE.mod.CI <- data.frame(Location=newdat$Location,SEmean.size=size.bio.SE.mod.CI[,1],SEmean.size.up=size.bio.SE.mod.CI[,2],SEmean.size.down=size.bio.SE.mod.CI[,3])

sizeSE.bio.ben <- left_join(size.bio.SE.mod.CI, coral.mod.pred, by="Location")
sizeSE.bio.ben <- left_join(sizeSE.bio.ben, macro.mod.pred, by="Location")
pairs(sizeSE.bio.ben[c(2,5:6)], lower.panel=panel.smooth, upper.panel=panel.cor, diag.panel=panel.hist)
sizeSE.bio.ben$cm <- log(sizeSE.bio.ben$coral/sizeSE.bio.ben$macro)


#output table
size_est_tbl <- left_join(size.bio.mod.CI,size.bio.P.mod.CI,by="Location")
size_est_tbl <- left_join(size_est_tbl,size.bio.SC.mod.CI,by="Location")
size_est_tbl <- left_join(size_est_tbl,size.bio.BR.mod.CI,by="Location")
size_est_tbl <- left_join(size_est_tbl,size.bio.GD.mod.CI,by="Location")
size_est_tbl <- left_join(size_est_tbl,size.bio.SE.mod.CI,by="Location")
str(size_est_tbl)

write.csv(size_est_tbl, 'outputs/SOM6.csv',row.names=F)

# panel fig - mean size ---------------------------------------------------
png('outputs/Fig4.png',height=2000,width=3600,res=300)
par(mfrow=c(3,6),mar=c(1,1,1,1),oma=c(4,3,0,0),xpd=F,mgp=c(1.6,.7,0))
plot.fish(size.bio.ben,'mean.size','coral',1,'size')
plot.fish(sizeP.bio.ben,'Pmean.size','coral',2,'size')
plot.fish(sizeSC.bio.ben,'SCmean.size','coral',3,'size')
plot.fish(sizeGD.bio.ben,'GDmean.size','coral',4,'size')
plot.fish(sizeBR.bio.ben,'BRmean.size','coral',5,'size')
plot.fish(sizeSE.bio.ben,'SEmean.size','coral',6,'size')

plot.fish(size.bio.ben,'mean.size','macro',7,'size')
plot.fish(sizeP.bio.ben,'Pmean.size','macro',8,'size')
plot.fish(sizeSC.bio.ben,'SCmean.size','macro',9,'size')
plot.fish(sizeGD.bio.ben,'GDmean.size','macro',10,'size')
plot.fish(sizeBR.bio.ben,'BRmean.size','macro',11,'size')
plot.fish(sizeSE.bio.ben,'SEmean.size','macro',12,'size')

plot.fish(size.bio.ben,'mean.size','cm',13,'size')
plot.fish(sizeP.bio.ben,'Pmean.size','cm',14,'size')
plot.fish(sizeSC.bio.ben,'SCmean.size','cm',15,'size')
plot.fish(sizeGD.bio.ben,'GDmean.size','cm',16,'size')
plot.fish(sizeBR.bio.ben,'BRmean.size','cm',17,'size')
plot.fish(sizeSE.bio.ben,'SEmean.size','cm',18,'size')


mtext("Coral Cover (%)", outer=T, side=2, at=0.83,line=1)
mtext("Macroalgal Cover (%)", outer=T, side=2, at=0.5,line=1)
mtext("log(Coral/Macrogalae)", outer=T, side=2, at=0.175,line=1)
mtext("Total",outer=T,side=1,at=0.09,line=1)
mtext("Predators",outer=T,side=1,at=0.248,line=1)
mtext("Sec. Consumers",outer=T,side=1,at=0.416,line=1)
mtext("Grazers",outer=T,side=1,at=0.584,line=1)
mtext("Browsers",outer=T,side=1,at=0.752,line=1)
mtext("Scrapers",outer=T,side=1,at=0.92,line=1)
mtext(expression("Biomass"~~bgroup("(","g"*m^{-2},")")),outer=T,side=1,line=2.8)

dev.off()

mod.out <- data.frame(
  y = as.character(c("coral", "coral", "coral", "coral", "coral",
                     "coral", "macro", "macro", "macro", "macro",
                     "macro", "macro", "cm", "cm", "cm", "cm", "cm", "cm")),
  x = as.character(c('mean.size','Pmean.size','SCmean.size','GDmean.size','BRmean.size','SEmean.size',
                     'mean.size','Pmean.size','SCmean.size','GDmean.size','BRmean.size','SEmean.size',
                     'mean.size','Pmean.size','SCmean.size','GDmean.size','BRmean.size','SEmean.size')),
  int = rep(NA,18),
  int.se = rep(NA,18),
  int.t = rep(NA,18),
  int.p = rep(NA,18),
  slp = rep(NA,18),
  slp.se = rep(NA,18),
  slp.t = rep(NA,18),
  slp.p = rep(NA,18),
  f = rep(NA,18),
  df = rep(NA,18),
  r2 = rep(NA,18), stringsAsFactors=F
)

for(i in c(1,7,13)){
  y <- mod.out[i,'y']
  x <- mod.out[i,'x']
  dat <- size.bio.ben
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(m.sel$class[1]=='lm'){
    mod.out[i,3] <- summary(tp)$coefficients[1]; mod.out[i,7] <- summary(tp)$coefficients[2]
    mod.out[i,4] <- summary(tp)$coefficients[3]; mod.out[i,8] <- summary(tp)$coefficients[4]
    mod.out[i,5] <- summary(tp)$coefficients[5]; mod.out[i,9] <- summary(tp)$coefficients[6]
    mod.out[i,6] <- summary(tp)$coefficients[7]; mod.out[i,10] <- summary(tp)$coefficients[8]
    mod.out[i,11] <- summary(tp)$fstatistic[1]
    mod.out[i,12] <- summary(tp)$df[2]
    mod.out[i,13] <- summary(tp)$r.squared
  }
  if(m.sel$class[1]=='gam'){
    mod.out[i,3] <- summary(tpg)$p.table[1]; mod.out[i,4] <- summary(tpg)$p.table[2]
    mod.out[i,5] <- summary(tpg)$p.table[3]; mod.out[i,6] <- summary(tpg)$p.table[4]
    mod.out[i,7] <- summary(tpg)$s.table[1]; mod.out[i,8] <- summary(tpg)$s.table[2]
    mod.out[i,10] <- summary(tpg)$s.table[4]; mod.out[i,11] <- summary(tpg)$s.table[3]
    mod.out[i,12] <- summary(tpg)$n
    mod.out[i,13] <- summary(tpg)$dev.expl
  }
}
sizelist <- list(size.bio.ben,sizeP.bio.ben,sizeSC.bio.ben,sizeGD.bio.ben,sizeBR.bio.ben,sizeSE.bio.ben)
for(i in c(2:6)){
  y <- mod.out[i,'y']
  x <- mod.out[i,'x']
  dat <- sizelist[[i]]
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(m.sel$class[1]=='lm'){
    mod.out[i,3] <- summary(tp)$coefficients[1]; mod.out[i,7] <- summary(tp)$coefficients[2]
    mod.out[i,4] <- summary(tp)$coefficients[3]; mod.out[i,8] <- summary(tp)$coefficients[4]
    mod.out[i,5] <- summary(tp)$coefficients[5]; mod.out[i,9] <- summary(tp)$coefficients[6]
    mod.out[i,6] <- summary(tp)$coefficients[7]; mod.out[i,10] <- summary(tp)$coefficients[8]
    mod.out[i,11] <- summary(tp)$fstatistic[1]
    mod.out[i,12] <- summary(tp)$df[2]
    mod.out[i,13] <- summary(tp)$r.squared
  }
  if(m.sel$class[1]=='gam'){
    mod.out[i,3] <- summary(tpg)$p.table[1]; mod.out[i,4] <- summary(tpg)$p.table[2]
    mod.out[i,5] <- summary(tpg)$p.table[3]; mod.out[i,6] <- summary(tpg)$p.table[4]
    mod.out[i,7] <- summary(tpg)$s.table[1]; mod.out[i,8] <- summary(tpg)$s.table[2]
    mod.out[i,10] <- summary(tpg)$s.table[4]; mod.out[i,11] <- summary(tpg)$s.table[3]
    mod.out[i,12] <- summary(tpg)$n
    mod.out[i,13] <- summary(tpg)$dev.expl
  }
}
for(i in c(8:12)){
  y <- mod.out[i,'y']
  x <- mod.out[i,'x']
  dat <- sizelist[[i-6]]
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(m.sel$class[1]=='lm'){
    mod.out[i,3] <- summary(tp)$coefficients[1]; mod.out[i,7] <- summary(tp)$coefficients[2]
    mod.out[i,4] <- summary(tp)$coefficients[3]; mod.out[i,8] <- summary(tp)$coefficients[4]
    mod.out[i,5] <- summary(tp)$coefficients[5]; mod.out[i,9] <- summary(tp)$coefficients[6]
    mod.out[i,6] <- summary(tp)$coefficients[7]; mod.out[i,10] <- summary(tp)$coefficients[8]
    mod.out[i,11] <- summary(tp)$fstatistic[1]
    mod.out[i,12] <- summary(tp)$df[2]
    mod.out[i,13] <- summary(tp)$r.squared
  }
  if(m.sel$class[1]=='gam'){
    mod.out[i,3] <- summary(tpg)$p.table[1]; mod.out[i,4] <- summary(tpg)$p.table[2]
    mod.out[i,5] <- summary(tpg)$p.table[3]; mod.out[i,6] <- summary(tpg)$p.table[4]
    mod.out[i,7] <- summary(tpg)$s.table[1]; mod.out[i,8] <- summary(tpg)$s.table[2]
    mod.out[i,10] <- summary(tpg)$s.table[4]; mod.out[i,11] <- summary(tpg)$s.table[3]
    mod.out[i,12] <- summary(tpg)$n
    mod.out[i,13] <- summary(tpg)$dev.expl
  }
}
for(i in c(14:18)){
  y <- mod.out[i,'y']
  x <- mod.out[i,'x']
  dat <- sizelist[[i-12]]
  temp <- dat[!is.na(dat[,y]),]
  tp <- lm(temp[,y]~temp[,x])
  cook <- cooks.distance(tp)
  temp.sub <- temp[!cook >= 0.50, ]
  tp <- lm(temp.sub[, y] ~ temp.sub[, x])
  tpg <- gam(temp.sub[, y] ~ s(temp.sub[, x], k = 3), family = gaussian)
  m.sel <- model.sel(tp, tpg)
  
  if(m.sel$class[1]=='lm'){
    mod.out[i,3] <- summary(tp)$coefficients[1]; mod.out[i,7] <- summary(tp)$coefficients[2]
    mod.out[i,4] <- summary(tp)$coefficients[3]; mod.out[i,8] <- summary(tp)$coefficients[4]
    mod.out[i,5] <- summary(tp)$coefficients[5]; mod.out[i,9] <- summary(tp)$coefficients[6]
    mod.out[i,6] <- summary(tp)$coefficients[7]; mod.out[i,10] <- summary(tp)$coefficients[8]
    mod.out[i,11] <- summary(tp)$fstatistic[1]
    mod.out[i,12] <- summary(tp)$df[2]
    mod.out[i,13] <- summary(tp)$r.squared
  }
  if(m.sel$class[1]=='gam'){
    mod.out[i,3] <- summary(tpg)$p.table[1]; mod.out[i,4] <- summary(tpg)$p.table[2]
    mod.out[i,5] <- summary(tpg)$p.table[3]; mod.out[i,6] <- summary(tpg)$p.table[4]
    mod.out[i,7] <- summary(tpg)$s.table[1]; mod.out[i,8] <- summary(tpg)$s.table[2]
    mod.out[i,10] <- summary(tpg)$s.table[4]; mod.out[i,11] <- summary(tpg)$s.table[3]
    mod.out[i,12] <- summary(tpg)$n
    mod.out[i,13] <- summary(tpg)$dev.expl
  }
}

write.csv(mod.out, 'outputs/Table2.csv',row.names=F)

# multivariate trophic biomass ---------------------------------------------------

temp <- trophic.benthic.m[!is.na(trophic.benthic.m$macro),]
temp <- temp[!is.na(temp$coral),]
temp$cm <- log(temp$coral/temp$macro)
temp$coral.p <- temp$coral/100
temp$macro.p <- temp$macro/100
spp.mat <- temp[c('BR','GD','P','SC','SE')]
std <- function(x){(x-min(x))/(max(x)-min(x))}
for(i in 1:5) spp.mat[i] <- std(spp.mat[i])
summary(spp.mat)

temp$coral.ps <- std(temp$coral.p)
temp$macro.ps <- std(temp$macro.p)

trop.rda <- rda((spp.mat),(temp[c('coral.ps','macro.ps')]))
summary(trop.rda)
anova(trop.rda)
permutest(trop.rda, permutations = 1000)
RsquareAdj(trop.rda)

c <- temp['coral.ps']
m <- temp['macro.ps']
cm <- as.matrix(temp[c('coral.ps','macro.ps')])
trop.rda.f <- rda(cm~.,data=spp.mat)
anova(trop.rda.f)
permutest(trop.rda.f, permutations = 1000)
RsquareAdj(trop.rda.f)
plot(trop.rda.f,scaling=3)
summary(trop.rda.f)

anova(trop.rda.f,by="terms")

# multivariate  - size ----------------------------------------------------

all.size <- left_join(sizeP.bio.ben[c("Location","Pmean.size","coral","macro")],sizeSC.bio.ben[c("Location","SCmean.size")],by="Location")
all.size <- left_join(all.size,sizeBR.bio.ben[c("Location","BRmean.size")],by="Location")
all.size <- left_join(all.size,sizeGD.bio.ben[c("Location","GDmean.size")],by="Location")
all.size <- left_join(all.size,sizeSE.bio.ben[c("Location","SEmean.size")],by="Location")
colnames(all.size) <- c("Location","P","coral","macro","SC","BR","GD","SE")
all.size <- all.size[c(1,3,4,2,5:8)]
str(all.size)
all.size <- all.size[!is.na(all.size$coral),]
all.size <- all.size[!is.na(all.size$macro),]

spp.mat <- all.size[4:8]
for(i in 1:5) spp.mat[i] <- std(spp.mat[i])
summary(spp.mat)

all.size$coral.ps <- std(all.size$coral/100)
all.size$macro.ps <- std(all.size$macro/100)

size.rda <- rda((spp.mat),(all.size[9:10]))
plot(size.rda,scaling=3)
summary(size.rda)
anova(size.rda)
permutest(trop.rda, permutations = 10000)
RsquareAdj(size.rda)


cm <- as.matrix(all.size[c(9:10)])
size.rda.f <- rda(cm~.,data=spp.mat)
anova(size.rda.f)
permutest(size.rda.f, permutations = 1000)
RsquareAdj(size.rda.f)
plot(size.rda.f,scaling=3)
summary(size.rda.f)

anova(size.rda.f,by="terms")
anova(size.rda.f,by="axis")
anova(size.rda.f,by="margin")

# combine multivariate plots ----------------------------------------------

png(file='outputs/Fig5.png',height=1800,width=3800,res=300)
par(mfrow=c(1,2),mar=c(3,3,2,1),oma=c(2,2,0,0))

plot(trop.rda.f,scaling=3,type="n",cex.axis=1.4,xlim=c(-1,1.05))
points(scores(trop.rda.f,display="sites",scaling=3),pch=21,col="black",bg="grey",cex=1.4)
reg.arrow <- scores(trop.rda.f, display="bp", scaling=3)
arrows(0,0,reg.arrow[,1],reg.arrow[,2],length=0,lty=1,cex=3,lwd=1,col="blue")
temp <- scores(trop.rda.f, display="b", scaling=3)
text(temp[1]+.25,temp[6]-0.06,labels=c("Browsers"),col="blue",cex=1.3)
text(temp[2]+.25,temp[7]-0.05,labels=c("Grazers"),col="blue",cex=1.3)
text(temp[3]+0.20,temp[8]-0.03,labels=c("Predators"),col="blue",cex=1.3)
text(temp[4]+.32,temp[9]-0.03,labels=c("Sec. Cons."),col="blue",cex=1.3)
text(temp[5]+0.23,temp[10],labels=c("Scrapers"),col="blue",cex=1.3)
reg.arrow <- scores(trop.rda.f, display="sp", scaling=3)
arrows(0,0,reg.arrow[,1],reg.arrow[,2],length=0,lty=1,cex=3,lwd=2,col="red")
temp <- scores(trop.rda.f, display="sp", scaling=3)
text(temp[1]-0.02,temp[3]+0.07,labels=c("Coral"),col="red",cex=1.5)
text(temp[2]-0.02,temp[4]+0.05,labels=c("Macroalgae"),col="red",cex=1.5)
text(-1.15,1,"A",cex=1.5)

plot(size.rda.f,scaling=3,type="n",cex.axis=1.4,xlim=c(-0.9,1.4),ylim=c(-.7,1))
points(scores(size.rda.f,display="sites",scaling=3),pch=21,col="black",bg="grey",cex=1.4)
reg.arrow <- scores(size.rda.f, display="bp", scaling=3)
arrows(0,0,reg.arrow[,1],reg.arrow[,2],length=0,lty=1,cex=3,lwd=1,col="blue")
temp <- scores(size.rda.f, display="b", scaling=3)
text(temp[1],temp[6]-0.05,labels=c("Predators"),col="blue",cex=1.3)
text(temp[2]+0.2,temp[7]+.06,labels=c("Sec. Cons."),col="blue",cex=1.3)
text(temp[3],temp[8]-0.05,labels=c("Browsers"),col="blue",cex=1.3)
text(temp[4]+0.1,temp[9]-0.05,labels=c("Grazers"),col="blue",cex=1.3)
text(temp[5]+0.2,temp[10]+.03,labels=c("Scrapers"),col="blue",cex=1.3)

reg.arrow <- scores(size.rda.f, display="sp", scaling=3)
arrows(0,0,reg.arrow[,1],reg.arrow[,2],length=0,lty=1,cex=3,lwd=2,col="red")
temp <- scores(size.rda.f, display="sp", scaling=3)
text(temp[1]+0.18,temp[3],labels=c("Coral"),col="red",cex=1.5)
text(temp[2],temp[4]+0.05,labels=c("Macroalgae"),col="red",cex=1.5)
text(-.80,1.05,"B",cex=1.5)

mtext("RDA1",side=1,outer=T,cex=1.4)
mtext("RDA2",side=2,outer=T,cex=1.4)
dev.off()

# supplemental figure random effects --------------------------------------
t.resid <- data.frame(DatasetID=trophic.bio.raw.m$DatasetID, BR=NA, GD=NA, P=NA, SC=NA, SE=NA, Year=trophic.bio.raw.m$Year, Location=trophic.bio.raw.m$Location)
for(k in c(2:6)){
  temp <- trophic.bio.raw.m
  colnames(temp)[k+4] <- "resp"
  tb.mod.l <- lmer(log(resp+1) ~ Location + (1|DatasetID) , data=temp)
  t.resid[,k] <- residuals(tb.mod.l)
}

plot(t.resid$DatasetID,t.resid$BR)

tot.resid <- data.frame(DatasetID=tot.bio.m$DatasetID, total = residuals(totbio.mod.l))

png(file='outputs/SOM7.png',height=4000,width=4000,res=300)
par(mfcol=c(6,2),mar=c(1.5,1.5,1,1.5),oma=c(4,4,0,0),mgp=c(1.6,.7,0),xpd=F)

plot(tot.resid$DatasetID,tot.resid$total, bty="o",cex.axis=1.4,ylim=c(-4,4))
abline(h=0)
text(0.2,0.9*4,"Total Biomass",pos=4,cex=2)
plot(t.resid$DatasetID,t.resid$P, bty="o",cex.axis=1.4,ylim=c(-3,5.5))
abline(h=0)
text(0.2,0.9*5.5,"Predator Biomass",pos=4,cex=2)
plot(t.resid$DatasetID,t.resid$SC, bty="o",cex.axis=1.4,ylim=c(-3,4.5))
abline(h=0)
text(0.2,0.9*4.5,"Secondary Consumer Biomass",pos=4,cex=2)
plot(t.resid$DatasetID,t.resid$BR, bty="o",cex.axis=1.4,ylim=c(-2,4.5))
abline(h=0)
text(0.2,0.9*4.5,"Browser Biomass",pos=4,cex=2)
plot(t.resid$DatasetID,t.resid$GD, bty="o",cex.axis=1.4, ylim=c(-2,4.5))
abline(h=0)
text(0.2,0.9*4.5,"Grazer Biomass",pos=4,cex=2)
plot(t.resid$DatasetID,t.resid$SE, bty="o",cex.axis=1.4,ylim=c(-2.2,4))
abline(h=0)
text(0.2,0.9*4,"Scraper Biomass",pos=4,cex=2)

plot(size.out$DatasetID,residuals(size.bio.mod), bty="o",cex.axis=1.4,ylim=c(-1.5,3.5))
abline(h=0)
text(0.2,0.9*3.5,"Mean size overall",pos=4,cex=2)
plot(size.bioP.out$DatasetID,residuals(size.bio.mod.P), bty="o",cex.axis=1.4,ylim=c(-2,3))
abline(h=0)
text(0.2,0.9*3,"Mean Size Predators",pos=4,cex=2)
plot(size.bioSC.out$DatasetID,residuals(size.bio.mod.SC), bty="o",cex.axis=1.4,ylim=c(-2,3))
abline(h=0)
text(0.2,0.9*3,"Mean Size Secondary Consumers",pos=4,cex=2)
plot(size.bioBR.out$DatasetID,residuals(size.bio.mod.BR), bty="o",cex.axis=1.4,ylim=c(-2,2))
abline(h=0)
text(0.2,0.9*2,"Mean Size Browsers",pos=4,cex=2)
plot(size.bioGD.out$DatasetID,residuals(size.bio.mod.GD), bty="o",cex.axis=1.4,ylim=c(-2,2.5))
abline(h=0)
text(0.2,0.9*2.5,"Mean Size Grazers",pos=4,cex=2)
plot(size.bioSE.out$DatasetID,residuals(size.bio.mod.SE), bty="o",cex.axis=1.4,ylim=c(-1.5,3))
abline(h=0)
text(0.2,0.9*3,"Mean Size Scrapers",pos=4,cex=2)

mtext("Dataset",side=1,outer=T,cex=1.6,line=2)
mtext("Residuals",side=2,outer=T,cex=1.6,line=2)

dev.off()

# correlations with lat/long? ---------------------------------------------

loc.coord <- fish.use %>% group_by(Location) %>% summarise("Lat"=median(Latitude,na.rm=T),"Long"=median(Longitude,na.rm=T)) %>% ungroup()


png(file='outputs/SOM9.png',height=3000,width=2800,res=300)
par(mfcol=c(6,4),mar=c(1,1,1,1),oma=c(4,4,2,0),mgp=c(1.6,.7,0),xpd=T)

####### biomass
temp <- data.frame(total = residuals(totbio.mod.l),Location=tot.bio.m$Location)
temp <- left_join(temp,loc.coord)
str(temp)

x <- temp$Lat; y <- temp$total
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- left_join(t.resid,loc.coord)
str(temp)

x <- temp$Lat; y <- temp$P
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Lat; y <- temp$SC
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Lat; y <- temp$GD
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Lat; y <- temp$BR
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Lat; y <- temp$SE
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,cex.axis=1.4); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(total = residuals(totbio.mod.l),Location=tot.bio.m$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$total
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- left_join(t.resid,loc.coord)
str(temp)
x <- temp$Long; y <- temp$P
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Long; y <- temp$SC
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Long; y <- temp$GD
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Long; y <- temp$BR
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,labels=NA); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

x <- temp$Long; y <- temp$SE
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-89,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}


######## size
temp <- data.frame(total = residuals(size.bio.mod),Location=size.out$Location)
temp <- left_join(temp,loc.coord)
str(temp)

x <- temp$Lat; y <- temp$total
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(P = residuals(size.bio.mod.P),Location=size.bioP.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Lat; y <- temp$P
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(SC = residuals(size.bio.mod.SC),Location=size.bioSC.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Lat; y <- temp$SC
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(GD = residuals(size.bio.mod.GD),Location=size.bioGD.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Lat; y <- temp$GD
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(BR = residuals(size.bio.mod.BR),Location=size.bioBR.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Lat; y <- temp$BR
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,cex.axis=1.4); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(SE = residuals(size.bio.mod.SE),Location=size.bioSE.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Lat; y <- temp$SE
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,cex.axis=1.4); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(13,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(total = residuals(size.bio.mod),Location=size.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$total
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(P = residuals(size.bio.mod.P),Location=size.bioP.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$P
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(SC = residuals(size.bio.mod.SC),Location=size.bioSC.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$SC
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(GD = residuals(size.bio.mod.GD),Location=size.bioGD.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$GD
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,labels=NA); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(BR = residuals(size.bio.mod.BR),Location=size.bioBR.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$BR
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(2,labels=NA); axis(1,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-87,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

temp <- data.frame(SE = residuals(size.bio.mod.SE),Location=size.bioSE.out$Location)
temp <- left_join(temp,loc.coord)
x <- temp$Long; y <- temp$SE
plot(x,y,pch=21,col="black",bg="grey",cex=1.5,bty="o",xlab="",ylab="",yaxt="n",xaxt="n",cex.axis=1.4); axis(1,cex.axis=1.4); axis(2,labels=NA); abline(h=0)
if((cor.test(x,y))$p.value < 0.05){
  text(-89,0.9*max(y),round(cor(x,y),2),pos=4,cex=2,font=2)
}

mtext("Latitude",outer=T,side=1,cex=1.3,at=.125,line=1)
mtext("Longitude",outer=T,side=1,cex=1.3,at=.375,line=1)

mtext("Latitude",outer=T,side=1,cex=1.3,at=0.625,line=1)
mtext("Longitude",outer=T,side=1,cex=1.3,at=0.875,line=1)

mtext("Total",outer=T,side=2,at=0.92,line=1,cex=1.3)
mtext("Predators",outer=T,side=2,at=0.752,line=1,cex=1.3)
mtext("Sec. Consumers",outer=T,side=2,at=0.584,line=1,cex=1.3)
mtext("Grazers",outer=T,side=2,at=0.416,line=1,cex=1.3)
mtext("Browsers",outer=T,side=2,at=0.248,line=1,cex=1.3)
mtext("Scrapers",outer=T,side=2,at=0.09,line=1,cex=1.3)

mtext(expression("Biomass"~~bgroup("(","g"*m^{-2},")")),outer=T,side=3,at=.25,cex=1.3,line=-1)
mtext("Mean size (cm)",outer=T,side=3,at=.75,cex=1.3)

dev.off()

# time effect? ------------------------------------------------------------
plot(residuals(totbio.mod.l)~tot.bio.m$Year)
plot(P~Year,data=t.resid)
plot(SC~Year,data=t.resid)
plot(GD~Year,data=t.resid)
plot(BR~Year,data=t.resid)
plot(SE~Year,data=t.resid)

png('outputs/SOM10.png',height=2800,width=1900,res=300)
par(mfrow=c(6,1),mar=c(1,1,1,1),oma=c(4,3,0,0),xpd=F,mgp=c(1.6,.7,0))

plot(residuals(totbio.mod.l)~tot.bio.m$Year,bty="l",xlab="",ylab="Total",cex.axis=1.2,cex.lab=1.1,xaxt="n")
lines(lowess(residuals(totbio.mod.l)~tot.bio.m$Year),lwd=4,col="red")
axis(1,labels=NA)

plot(P~Year,data=t.resid,bty="l",xlab="",ylab="Predators",cex.axis=1.2,cex.lab=1.1,xaxt="n")
lines(lowess(t.resid$P~t.resid$Year),lwd=4,col="red")
axis(1,labels=NA)

plot(SC~Year,data=t.resid,bty="l",xlab="",ylab="Sec. Consumers",cex.axis=1.2,cex.lab=1.1,xaxt="n")
lines(lowess(t.resid$SC~t.resid$Year),lwd=4,col="red")
axis(1,labels=NA)

plot(GD~Year,data=t.resid,bty="l",xlab="",ylab="Sec. Consumers",cex.axis=1.2,cex.lab=1.1,xaxt="n")
lines(lowess(t.resid$GD~t.resid$Year),lwd=4,col="red")
axis(1,labels=NA)

plot(BR~Year,data=t.resid,bty="l",xlab="",ylab="Sec. Consumers",cex.axis=1.2,cex.lab=1.1,xaxt="n")
lines(lowess(t.resid$BR~t.resid$Year),lwd=4,col="red")
axis(1,labels=NA)

plot(SE~Year,data=t.resid,bty="l",ylab="Sec. Consumers",cex.axis=1.2,cex.lab=1.1,xlab="Year")
lines(lowess(t.resid$SE~t.resid$Year),lwd=4,col="red")

mtext("Total",outer=T,side=2,at=0.92,line=1,cex=1)
mtext("Predators",outer=T,side=2,at=0.752,line=1,cex=1)
mtext("Sec. Consumers",outer=T,side=2,at=0.584,line=1,cex=1)
mtext("Grazers",outer=T,side=2,at=0.416,line=1,cex=1)
mtext("Browsers",outer=T,side=2,at=0.248,line=1,cex=1)
mtext("Scrapers",outer=T,side=2,at=0.09,line=1,cex=1)

mtext("Year",outer=T,side=1,line=1,cex=1)

dev.off()

# total herbivore biomass -------------------------------------------------
temp <- fish.use
temp$H <- ifelse(temp$Trophic=='GD'|temp$Trophic=='BR'|temp$Trophic=='SE','H','other')
H.bio.raw <- temp %>% 
  group_by(DatasetID,Replicate,Location,Year,H) %>% 
  summarise("sum"=sum(bio_use,na.rm=T)) %>% 
  spread(key=H,value=sum) %>% 
  ungroup()
H.bio.raw$H[is.na(H.bio.raw$H)] <- 0

H.bio.raw$DatasetID <- as.factor(H.bio.raw$DatasetID)
H.bio.raw$Location <- as.character(H.bio.raw$Location); H.bio.raw$Location <- as.factor(H.bio.raw$Location)

H.mod.l <- lmer(log(H+1) ~ Location + (1|DatasetID) , data=H.bio.raw)
plot(residuals(H.mod.l)~fitted(H.mod.l))
newdat <- data.frame(expand.grid(Location=unique(H.bio.raw$Location),DatasetID=11,Country="USVI"))
# tb.mod.pred <- predict(tb.mod.l,newdata=newdat)
# tb.mod.pred <- exp(tb.mod.pred) - 1
# trophic.loc.mod[k] <- tb.mod.pred
H.CI <- predictInterval(H.mod.l, newdata = newdat, n.sims = 10000, level=0.05)
H.CI <- exp(H.CI)-1

H.loc.mod <- H.CI
H.loc.mod$Location <- newdat$Location

H.loc.mod <- left_join(H.loc.mod, coral.mod.CI, by="Location")
H.loc.mod <- left_join(H.loc.mod, macro.mod.CI, by="Location")
str(H.loc.mod)

#### plot
png(file='outputs/SOM_herb.png',height=1200,width=3200,res=300)
par(mfrow=c(1,3),mar=c(3,4.5,2,1),oma=c(2,0,0,0),mgp=c(2,.7,0))

temp <- H.loc.mod[!is.na(H.loc.mod$coral),]
tp <-lm(temp$coral~temp$fit); summary(tp)
# ncvTest(tp); op <- par(mfrow=c(2,2),mar=c(4,4,2,1)); plot(tp); par(op) # non-constant error variance test
tpg <- gam(coral~s(fit,k=3),data=temp,family=gaussian); summary(tpg)
model.sel(tp,tpg)

tpg.pred <- predict(tpg, se.fit=T)
tpg.pred <- data.frame(x=temp$fit,y=tpg.pred$fit,up=tpg.pred$fit+1.96*tpg.pred$se.fit
                       ,down=tpg.pred$fit-1.96*tpg.pred$se.fit)
tpg.pred <- tpg.pred[with(tpg.pred,order(x)),]
# points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2,lty=2)

plot(temp$fit,temp$coral,pch=21,col="black",bg="grey",bty="l",xlab="",ylab="Coral Cover (%)",cex.axis=1.2,cex.lab=1.5,ylim=c(0,60),xlim=c(0,30),type="n")
# axis(1,labels=NA)

polygon(c(tpg.pred$x,rev(tpg.pred$x)),c(tpg.pred$up,rev(tpg.pred$down)),col=rgb(169,169,169,150,max=255),border=NA)
points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2)
# }
points(temp$fit,temp$coral,pch=21,col="black",bg="grey")


temp <- H.loc.mod[!is.na(H.loc.mod$macro),]
tp <-lm(temp$macro~temp$fit); summary(tp)
# ncvTest(tp); op <- par(mfrow=c(2,2),mar=c(4,4,2,1)); plot(tp); par(op) # non-constant error variance test
tpg <- gam(macro~s(fit,k=3),data=temp,family=gaussian); summary(tpg)
model.sel(tp,tpg)

tpg.pred <- predict(tpg, se.fit=T)
tpg.pred <- data.frame(x=temp$fit,y=tpg.pred$fit,up=tpg.pred$fit+1.96*tpg.pred$se.fit
                       ,down=tpg.pred$fit-1.96*tpg.pred$se.fit)
tpg.pred <- tpg.pred[with(tpg.pred,order(x)),]
# points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2,lty=2)

plot(temp$fit,temp$macro,pch=21,col="black",bg="grey",bty="l",xlab="",ylab="Macroalgal Cover (%)",cex.axis=1.2,cex.lab=1.5,ylim=c(0,60),xlim=c(0,30),type="n")
# axis(1,labels=NA)

polygon(c(tpg.pred$x,rev(tpg.pred$x)),c(tpg.pred$up,rev(tpg.pred$down)),col=rgb(169,169,169,150,max=255),border=NA)
points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2)
# }
points(temp$fit,temp$macro,pch=21,col="black",bg="grey")

temp <- H.loc.mod[!is.na(H.loc.mod$macro),]
temp <- temp[!is.na(temp$coral),]
temp$cm <- log(temp$coral/temp$macro)
tp <-lm(temp$cm~temp$fit); summary(tp)
# ncvTest(tp); op <- par(mfrow=c(2,2),mar=c(4,4,2,1)); plot(tp); par(op) # non-constant error variance test
tpg <- gam(cm~s(fit,k=3),data=temp,family=gaussian); summary(tpg)
model.sel(tp,tpg)

tpg.pred <- predict(tpg, se.fit=T)
tpg.pred <- data.frame(x=temp$fit,y=tpg.pred$fit,up=tpg.pred$fit+1.96*tpg.pred$se.fit
                       ,down=tpg.pred$fit-1.96*tpg.pred$se.fit)
tpg.pred <- tpg.pred[with(tpg.pred,order(x)),]
# points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2,lty=2)

plot(temp$fit,temp$cm,pch=21,col="black",bg="grey",bty="l",xlab="",ylab="log(Coral/Macroalgae)",cex.axis=1.2,cex.lab=1.5,xlim=c(0,30),type="n")
# axis(1,labels=NA)

polygon(c(tpg.pred$x,rev(tpg.pred$x)),c(tpg.pred$up,rev(tpg.pred$down)),col=rgb(169,169,169,150,max=255),border=NA)
points(tpg.pred$x,tpg.pred$y,type="l",col="black",lwd=2)
# }
points(temp$fit,temp$cm,pch=21,col="black",bg="grey")

mtext(expression("Total Herbivore Biomass"~~bgroup("(","g "*m^{-2},")")),side=1,outer=T,cex=1)
dev.off()

# SOM 1 & 3 -------------------------------------------------------------------
str(tot.bio.m)
temp <- tot.bio.m
temp$DatasetID <- as.character(temp$DatasetID); temp$DatasetID <- as.integer(temp$DatasetID)
meta.grp <- data.frame(
  DatasetID = c(11L, 33L, 40L, 47L, 48L, 56L, 59L, 60L, 200L, 201L, 202L,
                 203L, 205L, 206L, 207L, 209L, 210L, 211L, 212L, 216L, 218L,
                 223L, 262L, 581L, 699L, 700L, 700L),
  grp = as.character(c("1 Alan & Jim", "2 Alan CO", "3 NOAA LP", "5 UVI ",
                    "4 NOAA USVI", "6 Pete", "7 Marah", "7 Marah",
                    "8 AGRRA", "8 AGRRA", "8 AGRRA", "8 AGRRA", "8 AGRRA",
                    "8 AGRRA", "8 AGRRA", "8 AGRRA", "8 AGRRA", "8 AGRRA",
                    "8 AGRRA", "8 AGRRA", "8 AGRRA", "8 AGRRA",
                    "8 AGRRA", "7 Marah", "9 FL", "8 Waitt", "8 Waitt"))
)

temp <- left_join(temp,meta.grp,by='DatasetID')
str(temp)
temp %>% distinct(DatasetID,grp)
SOM1 <- temp %>% group_by(grp) %>% summarise('n'=length(tot.bio),'minYr'=min(Year),'maxYr'=max(Year))
sum(SOM1$n)
write.csv(SOM1, 'outputs/SOM1.csv',row.names=F)

som3 <- temp %>% group_by(Location) %>% summarise("n"=length((Replicate))) %>% ungroup()
temp <- temp %>% group_by(Location) %>% summarise("nDat"=length(unique(DatasetID))) %>% ungroup()
som3 <- full_join(som3,temp,by='Location')
write.csv(som3,'outputs/SOM3.csv',row.names = FALSE)
sum(som3$n)

# end ---------------------------------------------------------------------

Sys.time()-Start
