############       Fishing artificial reefs Lauwersoog

## the following two commands remove remnants from previous analyses and cleans the workspace
rm(list=ls())
graphics.off()


## you need to install and open the following libraries
library(FSAdata)  # for data
library(data.table)
library(tidyr)
library(vegan)
library(plyr)
library(dplyr)
library(lattice)
library(car)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggpubr)
library(reshape2)
library(lme4)
library(Rmisc)


#loading data

setwd("C:/RUG/Arbetsyta/1 Desktop/2022/REEFS-Lauwersoog - NEW/REEFS_sampling/Fishing/Analyses 2022_data")
setwd("D:/Arbetsyta/1 Desktop/2023/REEFS-Lauwersoog/REEFS_sampling/Fishing/Analyses 2023_data")
setwd("C:/Users/P310512/Documents/Groningen/Sound library")

######################## SHIETFUIKEN ############################

data <- read.table("data_table2.txt", header=T, sep="\t")
SH=data[which(data$gear=="schietfuik"),] # schietfuiken

###################### MULTIVARIATE ANALYSES

# 2023

SH=SH[which(SH$year==2023),]
SHK<- data[which(data$year==2023),]
# remove invertebrates and empty

fishSH=SH[, -which(names(SH)=="gewone.garnaal")]
fishSH=fishSH[,-which(names(fishSH)=="steurgarnaal")]
fishSH=fishSH[,-which(names(fishSH)=="strandkrab"),]
fishSH=fishSH[,-which(names(fishSH)=="gewone.zwemkrab")]
fishSH=fishSH[,-which(names(fishSH)=="wolhandkrab")]
fishSH=fishSH[,-which(names(fishSH)=="heremietkreeft")]
fishSH=fishSH[,-which(names(fishSH)=="nordzeekrab")]

fishSH=fishSH[,-which(names(fishSH)=="empty")]
fishSH=fishSH[,-which(names(fishSH)=="na")]

fishSH=fishSH[,-which(names(fishSH)=="grondel")] # one unidentified individual
fishSH=fishSH[,-which(names(fishSH)=="snoekbaars")]
fishSH=fishSH[,-which(names(fishSH)=="sprat")]
fishSH=fishSH[,-which(names(fishSH)=="koornaarvis")]
fishSH=fishSH[,-which(names(fishSH)=="zeedonderpad")]

fishSHK = SHK[,which(names(SHK)!="steurgarnaal" & names(SHK)!="strandkrab" & names(SHK)!="gewone zwemkrab" & names(SHK)!="gewone.zwemkrab" & 
                       names(SHK)!="wolhandkrab" &
                      names(SHK)!="heremietkreeft" & names(SHK)!="noordzeekrab" & names(SHK)!="empty" & names(SHK)!="na" & 
                       names(SHK)!="snoekbaars" & names(SHK)!="sprat")]
# checking for empty rows and deleting them
fishSH$S=specnumber(fishSH[,12:32])
fishSH=fishSH[-which(fishSH$haul_ID=="SF_124"), ]
fishSH = subset(fishSH, month2!="August")

fishSHK$S<-specnumber(fishSHK[,12:37])
fishSHKsub<-subset(fishSHK,S>0)
fishSHKsub<-subset(fishSHKsub,month2!="August")

#Matching fykes to sound locations
#fishSHsub<-subset(fishSH,treatment_ID=="SF_1" | treatment_ID == "SF_5" | treatment_ID == "SF_3" | treatment_ID == "SF_4" | treatment_ID == "SF_9")
fishSHsub<-subset(fishSH,block!="BE")
fishSHsub[6,12:32]<-fishSHsub[6,12:32]+fishSHsub[7,12:32]
fishSHsub[10,12:32]<-fishSHsub[10,12:32]+fishSHsub[11,12:32]
fishSHsub[13,12:32]<-fishSHsub[13,12:32]+fishSHsub[14,12:32]
fishSHsub[15,12:32]<-fishSHsub[15,12:32]+fishSHsub[16,12:32]
fishSHsub[17,12:32]<-fishSHsub[17,12:32]+fishSHsub[18,12:32]
fishSHsub2<-fishSHsub[c(1:6,8:10,12:13,15,17),]
fishSHsub2$treatment_ID[fishSHsub2$treatment_ID=="SF_1"]<-"W_on"
fishSHsub2$treatment_ID[fishSHsub2$treatment_ID=="SF_2"]<-"W_on"
fishSHsub2$treatment_ID[fishSHsub2$treatment_ID=="SF_5"]<-"E_on"
fishSHsub2$treatment_ID[fishSHsub2$treatment_ID=="SF_3"]<-"E_off"
fishSHsub2$treatment_ID[fishSHsub2$treatment_ID=="SF_4"]<-"E_off"
fishSHsub2$treatment_ID[fishSHsub2$treatment_ID=="SF_9"]<-"W_off"
fishSHsub2$treatment_ID2<-as.factor(fishSHsub2$treatment_ID)
#Create sound file
soundtypes<-unique(nonoise$Comments)
soundtypes
allsounds<-data.frame(location = c("E_off","E_on","W_off","W_on","E_off","E_on","W_off","W_on","E_off","E_on","W_off","W_on","E_off","E_on","W_off","W_on"),
                         month = c("February","February","February","February","April","April","April","April","May","May","May","May","June","June","June","June"))

for(n in 1:length(soundtypes)){
  spec<-soundtypes[n]
  sub<-subset(nonoise,Comments==spec)
  tb<-table(sub$location,sub$month2)
  soundcatch<-as.data.frame(tb)
  allsounds[,n+2]<-soundcatch[1:16,3]
}

names(allsounds)<-c("location","month","drum_vibrations","bubble_growl","clicks","striped_grunt","fast_clicking","fish_croacking","grunt_2","higher_grunt",
                    "grunt_4","hollow_drumming","hollow_growl","low_drumming","bursting","water_croacking","wiping","knocking","hard_knocking","low_creaking",
                    "lower_knocking","bowl_growl","low_growl","vibrating_sound","seal_sound","chirping","baby_seal")
#soundtable<-spread(nonoise,Comments,Start)
#allsounds2<-cbind(allsounds,allsoundcatch[,1:2])

#Remove rare sounds
allsounds3<-allsounds[,c(1:3,5,6,10,12,14,15,18,19,21:24)]
allsounds4<-allsounds3[order(allsounds3$location),]
allsoundssub<-allsounds4[-c(6,11,15),]
fishSHsub3<-fishSHsub2[order(fishSHsub2$treatment_ID2),]
#allsounds3<-rbind(allsounds2[3,],allsounds2[1,],allsounds2[4,],allsounds2[2,],allsounds2[7,],allsounds2[8,],allsounds2[5],allsounds2[9,],allsounds2[10,],allsounds2[15,],allsounds2[16,],allsounds2[13,],allsounds2[14,])
##allsounds2$Var1[allsounds2$Var1=="W_on"]<-"reef"
#allsounds2$Var1[allsounds2$Var1=="E_on"]<-"reef"
#allsounds2$Var1[allsounds2$Var1=="W_off"]<-"no_reef"

# multivariate datafiles

species=fishSHKsub[,12:37] # just a way to get an easy command to the species data, you can also just use sounds[,6:16]
env=fishSHKsub[,c(3:11,39:43)] #  just a way to get an easy command to the factors data, you can also just use sounds[,1:4]
#env = envfit(ord, allsounds2, permutations = 999, na.rm = TRUE)
species = fishSHsub3[,12:32]
names(species)<-c("flounder","rock gunnel","common goby","sand goby","three-spined stickleback","greater pipefish","herring","cod","lesser pipefish",
                  "five-bearded rockling","European eel","eelpout","dab","plaice","smelt","pouting","sole","whiting","seabass","turbot","lemon sole")
env = allsoundssub

env$month[env$month=="April"]<-"May"
env$month[env$month=="May"]<-"May"
##### MULTVARIATES ######

ord <- metaMDS(species)
ordiplot(ord, type = "t")
ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))

### misc figs

# Fig 1a - reef area or not
plot(ord, display = "species", type = "n", ylim=c(-1.5,1))
points(ord, display = "sites", pch = c(2, 1) [as.factor(env$reef_type_2)])
legend("topleft", c("on reef", "off reef"), pch = c(1,2), title = "on/off reef")


plot(ord, display = "species", type = "n", ylim=c(-1.5,1))
points(ord, display = "sites", pch = c(2, 1,3,4) [as.factor(env$Var1)])
legend("topleft", c("E_off", "E_on","W_off","W_on"), pch = c(2,1,3,4), title = "Location")

# Fig 1b - reef area or not
ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))
# with(env, ordiellipse(ord, Site, kind = "sd", col=2:3, draw = "polygon", label = TRUE))
with(env, ordiellipse(ord, reef_type_2, kind="sd",col="green3", draw = "polygon",  label = TRUE, show.groups = "reef"))
with(env, ordiellipse(ord, reef_type_2, kind="sd",col="darkviolet", draw = "polygon", label = TRUE, show.groups = "no_reef"))


# Fig 1b - reef area or not
ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))
# with(env, ordiellipse(ord, Site, kind = "sd", col=2:3, draw = "polygon", label = TRUE))
with(env, ordiellipse(ord, Var1, kind="sd", draw = "polygon", col='blue', label = TRUE,show.groups = 'E_off'))
with(env, ordiellipse(ord, Var1, kind="sd",col="darkviolet", draw = "polygon", label = TRUE, show.groups = "E_on"))
with(env, ordiellipse(ord, Var1, kind="sd",col="darkgreen", draw = "polygon", label = TRUE, show.groups = "W_off"))
with(env, ordiellipse(ord, Var1, kind="sd",col="red", draw = "polygon", label = TRUE, show.groups = "W_on"))

# Fig 1c - season

ordiplot(ord, display = "species", type = "none")
points(ord, display = "sites", pch = c(3, 2, 1) [as.factor(env$season)])
legend("topleft", c("winter", "spring", "summer"), pch = c(1,3,2), title = "season")

ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))
with(env, ordiellipse(ord, season, kind="sd",col="white", draw = "polygon",  label = TRUE, show.groups = "winter"))
with(env, ordiellipse(ord, season, kind="sd",col="green3", draw = "polygon",  label = TRUE, show.groups = "spring"))
with(env, ordiellipse(ord, season, kind="sd",col="red", draw = "polygon",  label = TRUE, show.groups = "summer"))


#Fig 1d - months
ordiplot(ord, display = "species", type = "none")
points(ord, display = "sites", pch = c(4, 3, 2, 1) [as.factor(env$Var2)])
legend("topleft", c("February", "April", "May","June"), pch = c(4,3,2,1), title = "season")

ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))
with(env, ordiellipse(ord, month2, kind="sd",col="white", draw = "polygon",  label = TRUE, show.groups = "February"))
with(env, ordiellipse(ord, month2, kind="sd",col="green3", draw = "polygon",  label = TRUE, show.groups = "April"))
with(env, ordiellipse(ord, month2, kind="sd",col="orange", draw = "polygon",  label = TRUE, show.groups = "May"))
with(env, ordiellipse(ord, month2, kind="sd",col="red", draw = "polygon",  label = TRUE, show.groups = "June"))

env2<-allsoundssub[,3:15]
sigsounds<-allsoundssub[,c(3:5,11)]
en = envfit(ord, env2, permutations = 999, na.rm = TRUE)

ordiplot(ord,display = "species", type='none')
with(species,text(ord,disp="sp"))
plot(en)
with(env, ordiellipse(ord, month, kind="sd",col="white", draw = "polygon",  label = TRUE, show.groups = "February"))
with(env, ordiellipse(ord, month, kind="sd",col="green3", draw = "polygon",  label = TRUE, show.groups = "May"))
with(env, ordiellipse(ord, month, kind="sd",col="orange", draw = "polygon",  label = TRUE, show.groups = "May2"))
with(env, ordiellipse(ord, month, kind="sd",col="red", draw = "polygon",  label = TRUE, show.groups = "June"))


#Making an ordination plot for sounds
sounds=allsounds3[,3:15] # just a way to get an easy command to the species data, you can also just use sounds[,6:16]
env=allsounds[,1:2] #  just a way to get an easy command to the factors data, you can also just use sounds[,1:4]

ord <- metaMDS(sounds)
ordiplot(ord, type = "t")
ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))
with(env, ordiellipse(ord, month, kind="sd",col="white", draw = "polygon",  label = TRUE, show.groups = "February"))
with(env, ordiellipse(ord, month, kind="sd",col="green3", draw = "polygon",  label = TRUE, show.groups = "April"))
with(env, ordiellipse(ord, month, kind="sd",col="orange", draw = "polygon",  label = TRUE, show.groups = "May"))
with(env, ordiellipse(ord, month, kind="sd",col="red", draw = "polygon",  label = TRUE, show.groups = "June"))


#Trying a new way of blocking the data to match the sound recordings
require(tidyverse)
fishSHKsub$month2[fishSHKsub$month2=="April"]<-"May"
fishSHKsub$block[fishSHKsub$block=="B1"]<-"W_on"
fishSHKsub$block[fishSHKsub$block=="B2"]<-"E_off"
fishSHKsub$block[fishSHKsub$block=="B3"]<-"E_on"
fishSHKsub$block[fishSHKsub$block=="BW"]<-"W_off"

blockspec<- fishSHKsub %>%
  group_by(month2,block) %>%
  summarize(flounder = sum(bot),rock.gunnel = sum(botervis),common.goby = sum(brakwatergrondel),sand.goby = sum(dikkopje),
            threespined.stickleback = sum(driedoornige.stekelbaars),greater.pipefish = sum(grote.zeenaald),herring = sum(haring),cod = sum(kabeljauw),
  lesser.pipefish = sum(kleine.zeenaald), fivebearded.rockling = sum(meun.vijfdradige),European.eel = sum(paling),eelpout = sum(puitaal),dab = sum(schar),plaice = sum(schol),smelt = sum(spiering),pouting = sum(steenbolk),sole = sum(tong),whiting = sum(wijting),
  seabass = sum(zeebaars),northorn.sculpin = sum(zeedonderpad), sand.smelt = sum(koornaarvis),turbot = sum(Tarbot),lemon.sole = sum(tongschar))

blockspec_sub<-subset(blockspec,block!="BE")

allsounds$month[allsounds$month=="April"]<-"May"

blocksounds<- allsounds %>%
  group_by(month,location) %>%
  summarize(drum_vibrations = sum(drum_vibrations),bubble_growl = sum(bubble_growl),clicks = sum(clicks),striped_grunt = sum(striped_grunt),
            fast_clicking = sum(fast_clicking), fish_croacking = sum(fish_croacking), higher_grunt = sum(higher_grunt),hollow_drumming = sum(hollow_drumming),
            hollow_growl = sum(hollow_growl),low_drumming = sum(low_drumming),bursting = sum(bursting),water_croacking = sum(water_croacking),
            wiping = sum(wiping),knocking = sum(knocking),hard_knocking = sum(hard_knocking),low_creaking = sum(low_creaking),
            lower_knocking = sum(lower_knocking),bowl_growl = sum(bowl_growl),low_growl = sum(low_growl), vibrating_sound = sum(vibrating_sound),
            chirping = sum(chirping))

#New NMDS plot with correct sample dates
env = blocksounds[,1:2]
env2=blocksounds[,c(3,4:23)] # just a way to get an easy command to the species data, you can also just use sounds[,6:16]
species=blockspec_sub[,3:25] #  just a way to get an easy command to the factors data, you can also just use sounds[,1:4]

ord <- metaMDS(species)
en = envfit(ord, env2, permutations = 999, na.rm = TRUE)

ordiplot(ord, display = "species", type = "none")
with(species, text(ord, disp = "sp"))
plot(en)
with(env, ordiellipse(ord, month, kind="sd",col="white", draw = "polygon",  label = TRUE, show.groups = "February"))
with(env, ordiellipse(ord, month, kind="sd",col="orange", draw = "polygon",  label = TRUE, show.groups = "May"))
with(env, ordiellipse(ord, month, kind="sd",col="red", draw = "polygon",  label = TRUE, show.groups = "June"))


### Fish catches vs sounds for specific sound types ###

#Lower knocking - European seabass
require(patchwork)
mytheme<-theme(axis.text.x = element_text(colour="black",size=14), 
               axis.text.y=element_text(size=14,colour="black"), axis.title.x=element_text(size=16,colour="black"), 
               axis.title.y=element_text(size=16,colour="black"), axis.title.y.right=element_text(size=16,colour="#999900"), panel.grid = element_line(colour=colors()[355]), panel.background=element_rect(fill= "transparent",colour =NA),
               panel.border=element_rect(colour="black",fill="transparent"),
               legend.background = element_rect(fill="transparent",colour=NA), axis.line.x = element_line(colour="black"),legend.title = element_text(size=16),legend.text = element_text(size=14))

fishSHKsub$month2<-factor(fishSHKsub$month2,levels=c("February","May","June"))

g1<-ggplot(data=fishSHKsub,aes(x=month2,y=zeebaars))+geom_bar(stat="identity",fill="#336666")+labs(title="Fish catches",x="Month",y="Sea bass")+mytheme
g2<-ggplot(data=allsounds,aes(x=month,y=lower_knocking))+geom_bar(stat="identity",fill="#66CCCC")+labs(title="Sounds",x="Month",y="Detections")+mytheme
g1/g2
