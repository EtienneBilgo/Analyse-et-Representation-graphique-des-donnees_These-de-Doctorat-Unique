# Analyse-et-Representation-graphique-des-donnees_These-de-Doctorat-Unique
Annexe 1 : R Codes détaillés du traitement, des représentations graphiques et de l’analyse des  données de la virulence et du test tunnel avec Chromobacterium violaceum contre les moustiques Anopheles coluzzii
#Codes des analyses virulence Chromobacterium # violaceum et test tunnel
#Etienne Bilgo,
#IRSS/Centre Muraz,
# bilgo02@yahoo.fr
# Data importation
dat <- read.csv("~/Desktop/Thesis data_PhD_Def/virulence_CSP_def.csv", sep=";")
View(dat)
colnames(dat)
# Data manupilations
colnames(dat)
library(reshape)
#melt(Data set, column names to keep, column names to restructure)
colnames(dat)
df=melt(dat, colnames(dat)[c(1,5)], colnames(dat)[2:4])
View(df)
library(plyr)
df2=ddply(df, .(Replicate,variable), transform, nval=sum(value, na.rm=TRUE))
View(df2)
df3=ddply(df2, .(Replicate,variable), transform, Dead=cumsum(value))
View(df3)
df3$Survival=(1-df3$Dead/df3$nval)*100
View(df3)
colnames(df3)[3]="Traitements"
View(df3)
df5=subset(df3, Day!="Alives")
df5$Day=as.numeric(as.character(df5$Day))
View(df5)
df6=ddply(df5, .(Day,Traitements), summarize, mean=mean(Survival), replicates=length(Survival), se=sd(Survival)/sqrt(length(Survival)))
View(df6)
library(ggplot2)
limits=aes(ymax=mean+se, ymin=mean-se)
theme = theme_bw()+theme(text = element_text(size=25), axis.title.x = element_text(size=25), axis.title.x = element_text(size=25), title = element_text(size=25), legend.title = element_text(size=25), legend.text = element_text(size=20))
cbPalette <-c("#105576","#9ed900","#ff7373")
Plt=ggplot(df6, aes(Day, mean, color=Traitements))+geom_line(size=2)+geom_errorbar(limits, width=.1, size=1)+theme+scale_colour_manual(values=cbPalette)+ylab("Moyenne survie (%)")+xlab("Jours")
#Displays the plot
Plt

# Data analysis tunnel test_chromobacterium 
data_tunnel <- read.csv("~/Desktop/data_tunnel.csv", sep=";")
View(data_tunnel)
colnames(data_tunnel)
dat=data_tunnel[c(2,4:9)]
View(dat)
colnames(dat)
str(dat)
colnames(dat)[4]="Empty"
colnames(dat)[5]="Status"
colnames(dat)[6]="GP"
colnames(dat)[7]="Vitality"
colnames(dat)
levels(dat$Status)
levels(dat$Status)=c("Sugarfed", "Sugarfed", "Unfed", "Unfed")
levels(dat$Status)
levels(dat$Vitality)=c("Alive", "Dead", "Alive", "Dead")
View(dat)
colnames(dat)
library(reshape2)
dat=melt(dat, colnames(dat)[c(1:3, 5, 7)], colnames(dat)[c(4,6)])
View(dat)
library(plyr)
colnames(dat)[2]="DaysPostInfection"
dat=ddply(subset(dat, Vitality=="Alive"), .(Replicate, DaysPostInfection, Treatment, variable), summarize, value=sum(value))
View(dat)
dat2=ddply(dat, .(Replicate, DaysPostInfection, Treatment), transform, Percent=value/sum(value))
View(dat2)
dat3=ddply(dat2, .(DaysPostInfection, Treatment, variable), summarize, Mosquitoes=mean(Percent), se=(sd(Percent)/sqrt(length(Percent))), Replicates=length(value))
View(dat3)
# We can start the t.test
subsetters=dat$Treatment=="control"&dat$variable=="GP"
testdat=dat[subsetters,]
pairwise.t.test(testdat$value, testdat$DaysPostInfection)
subsetters=dat$Treatment=="CSP"&dat$variable=="GP"
testdat=dat[subsetters,]
pairwise.t.test(testdat$value, testdat$DaysPostInfection)
## Results Pairwise t.test
#Pairwise comparisons using t tests with pooled SD 
# data:  testdat$value and testdat$DaysPostInfection 
#  6    7    8   
#7 0.93 -    -   
#8 1.00 1.00 -   
#9 0.54 1.00 0.93
#P value adjustment method: holm
subsetters=dat$DaysPostInfection==6 & dat$variable=="GP"
testdat=dat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment)
# Results below
#airwise comparisons using t tests with pooled SD 
# data:  testdat$value and testdat$Treatment 
# CSP 
# control 0.14
#P value adjustment method: holm 
subsetters=dat$DaysPostInfection==7 & dat$variable=="GP"
testdat=dat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment)
#Results below
#Pairwise comparisons using t tests with pooled SD 
#data:  testdat$value and testdat$Treatment 
#       CSP  
#control 0.063
#P value adjustment method: holm 
subsetters=dat$DaysPostInfection==8 & dat$variable=="GP"
testdat=dat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment)
#Pairwise comparisons using t tests with pooled SD 
#data:  testdat$value and testdat$Treatment 
# CSP  
# control 0.052
#P value adjustment method: holm 
subsetters=dat$DaysPostInfection==9 & dat$variable=="GP"
testdat=dat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment)
#Pairwise comparisons using t tests with pooled SD 
#data:  testdat$value and testdat$Treatment 
#CSP 
#control 0.02 # P<0.05 here
#P value adjustment method: holm 
subsetters=dat$variable=="GP"
testdat=dat[subsetters,]
summary(aov(testdat$value~testdat$Treatment*testdat$DaysPostInfection))
#                                           Df Sum Sq Mean Sq F value   Pr(>F)    
#testdat$Treatment                            1 1152.0  1152.0   24.90 2.85e-05 ***
#testdat$DaysPostInfection                    1   34.2    34.2    0.74    0.397    
#testdat$Treatment:testdat$DaysPostInfection  1   34.2    34.2    0.74    0.397    
#Residuals                                   28 1295.5    46.3                     
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#Plot according days post infection
library(ggplot2)
library(scales)
plt=ggplot(dat, aes(Treatment, value, fill=variable))+geom_bar(stat="identity", position="dodge")+facet_wrap(~DaysPostInfection)
plt
limits=aes(ymax=Mosquitoes+se, ymin=Mosquitoes-se)
theme = theme_bw()+theme(text = element_text(size=25),axis.title.x = element_text(size=30), title = element_text(size=35))
cbPalette <- c("#cc0000","#007acc")
plt=ggplot(subset(dat3, variable=="GP"), aes(DaysPostInfection, Mosquitoes, color=Treatment))+geom_line(size=3)+geom_errorbar(limits, width=.2, size=2)+theme+scale_colour_manual(values=cbPalette)+xlab("Days post bacteria infection")+ylab("Blood-feeding interest")+ggtitle("Impact of chromobacterium infection on\nblood-feeding interest")+scale_y_continuous(labels=percent)
plt
library(plotly)
ggplotly(plt)

Annexe 2 : R Codes du traitement, des représentations graphiques et de l’analyse des  données du chapitre 7 

# Etienne Bilgo,
#IRSS/Centre Muraz,
# bilgo02@yahoo.fr
####Load data and install packages####
#Mortality data
Mdat <- read.csv("Mortality.csv", header=TRUE)
#Tunnel data
Tdat <- read.csv("Tunnel.csv", header=TRUE)
#PCR data
#Toxin arsenal data
TxArdat <- read.csv("Toxin Arsenal.csv", header=TRUE)
#Spore dose data
#Load packages
#install.packages(c("reshape","plyr","ggplot2", "scales", "MASS"))
#uncomment to load packages
library(reshape)
library(plyr)
library(ggplot2)
library(scales)
library(MASS)
####Mortality Analysis####
#Standardize species names
levels(Mdat$Species)
levels(Mdat$Species)[1]="An.coluzzii"
levels(Mdat$Species)[2]="An.gambiae ss"
#Reshape and transpose data
Mdat2=melt(Mdat, colnames(Mdat)[c(1,11)], colnames(Mdat)[c(2:10)])
colnames(Mdat2)[3]="Treatment"
#Split combined variable
Mdat2$Replicate=Mdat2$Treatment
levels(Mdat2$Treatment)[grep("RFP", levels(Mdat2$Treatment))]="RFP"
levels(Mdat2$Treatment)[grep("Hybrid", levels(Mdat2$Treatment))]="Hybrid"
levels(Mdat2$Treatment)[grep("Controle",
levels(Mdat2$Treatment))]="Control"
levels(Mdat2$Replicate)[grep("1", levels(Mdat2$Replicate))]="1"
levels(Mdat2$Replicate)[grep("2", levels(Mdat2$Replicate))]="2"
levels(Mdat2$Replicate)[grep("3", levels(Mdat2$Replicate))]="3"
#Analyze survival from mortality data
Mdat3=ddply(Mdat2, .(Replicate, Treatment, Species), transform,
Percent=1-cumsum(value)/sum(value), n=sum(value))
Mdat3=subset(Mdat3, day!="Alive")
Mdat3$day=as.numeric(as.character(Mdat3$day))
Mdat4=ddply(Mdat3, .(Species, Treatment, day), summarize,
mean=mean(Percent), se=sd(Percent)/sqrt(length(Percent)))
colnames(Mdat4)[3]="Day"
#Calculate LT50 for each replicate with an LT50
LTdat=Mdat3[c(1:7)]
colnames(LTdat)[1]="Day"
attach(LTdat)
LTdat$Alive=LTdat$Percent*LTdat$n
LTdat$Dead=LTdat$n-LTdat$Alive
surv.per=0.50
LTdat2=ddply(LTdat,.(Species, Treatment, Replicate), summarize,
LT=as.numeric(dose.p(glm(cbind(Alive,Dead)~Day,binomial),p=surv.per)))
LTdat2[LTdat2$LT<0,]$LT=NA
LT50.Error=ddply(LTdat2, .(Treatment, Species), summarize, "LT50
Mean"=mean(LT,na.rm=T), se=sd(LT, na.rm=T)/sqrt(length(LT[!is.na(LT)])),
Replicates=length(LT[!is.na(LT)]))
#Run t tests for LT50s
testdat=LTdat2
subsetters=testdat$Species=="An.coluzzii"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$LT, testdat$Treatment,  p.adj="none")
testdat=LTdat2
subsetters=testdat$Species=="An.gambiae ss"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$LT, testdat$Treatment,  p.adj="none")
testdat=LTdat2
subsetters=testdat$Species=="An.kisumu"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$LT, testdat$Treatment,  p.adj="none")
testdat=LTdat2
summary(aov(testdat$LT~testdat$Treatment*testdat$Species))
#Calculate LT80 for each replicate with an LT80
surv.per=0.20
LTdat2=ddply(LTdat,.(Species, Treatment, Replicate), summarize,
LT=as.numeric(dose.p(glm(cbind(Alive,Dead)~Day,binomial),p=surv.per)))
LTdat2[LTdat2$LT<0,]$LT=NA
LT80.Error=ddply(LTdat2, .(Treatment, Species), summarize, "LT50
Mean"=mean(LT,na.rm=T), se=sd(LT, na.rm=T)/sqrt(length(LT[!is.na(LT)])),
Replicates=length(LT[!is.na(LT)]))
#Plot survival data
limits=aes(ymax=mean+se, ymin=mean-se)
theme = theme_bw()+theme(text = element_text(size=20), axis.title.x =
element_text(size=30), axis.text.x = element_text(angle=90, hjust=1,
vjust=.5, size=20), axis.text.y = element_text(size=25), title =
element_text(size=35), legend.title = element_text(size=25), legend.text
= element_text(size=20))
cbPalette <- c("#cc0000", "#29a329", "#007acc")
m.plt1=ggplot(Mdat4, aes(Day, mean, color=Treatment))+geom_line(size=3)+
  geom_errorbar(limits, width=.4,
size=2)+theme+scale_colour_manual(values=cbPalette)+
  xlab("Days after exposure")+ylab("Percent survival")+ggtitle("Survival
following WHO tube exposure")+
scale_y_continuous(labels=percent)+scale_x_continuous(breaks=0:max(Mdat4$
Day))+facet_wrap(~Species)
m.plt1
####Tunnel Analysis####
#Manipulate and rename data
Tdat=Tdat[c(2,4:8)]
colnames(Tdat)=c("Replicate", "DaysPostInfection", "Treatment", "Empty",
"Status", "GuineaPig")
#Split combined variable
Tdat2=melt(Tdat, colnames(Tdat)[c(1:3, 5)], colnames(Tdat)[c(4,6)])
Tdat2$Vitality=Tdat2$Status
levels(Tdat2$Vitality)[grep("dea", levels(Tdat2$Vitality))]="Dead"
levels(Tdat2$Vitality)[grep("aliv", levels(Tdat2$Vitality))]="Alive"
levels(Tdat2$Status)[grep("sug", levels(Tdat2$Status))]="Sugarfed"
levels(Tdat2$Status)[grep("un", levels(Tdat2$Status))]="Unfed"
Tdat2$Replicate=as.factor(Tdat2$Replicate)
levels(Tdat2$Treatment)=c("Control", "Hybrid", "No Host", "RFP")
#Calculate host preferences
Tdat3=ddply(subset(Tdat2, Vitality=="Alive"), .(Replicate,
DaysPostInfection, Treatment, variable), summarize, value=sum(value))
Tdat4=ddply(Tdat3, .(Replicate, DaysPostInfection, Treatment), transform,
Percent=value/sum(value))
Tdat5=ddply(Tdat4, .(DaysPostInfection, Treatment, variable), summarize,
Mosquitoes=mean(Percent), se=(sd(Percent)/sqrt(length(Percent))),
Replicates=length(value))
#Run t test on preferences
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed")
subsetters=testdat$Treatment=="Control"&testdat$variable=="GuineaPig"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$DaysPostInfection,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed")
subsetters=testdat$DaysPostInfection==5&testdat$variable=="GuineaPig"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed")
subsetters=testdat$Treatment=="Hybrid"&testdat$variable=="GuineaPig"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$DaysPostInfection,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed")
subsetters=testdat$variable=="GuineaPig"
testdat=testdat[subsetters,]
summary(aov(testdat$value~testdat$Treatment*testdat$DaysPostInfection))
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed"&
variable=="GuineaPig")
subsetters=testdat$DaysPostInfection==0|testdat$DaysPostInfection==1
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed"&
variable=="GuineaPig")
subsetters=testdat$DaysPostInfection==0|testdat$DaysPostInfection==2
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed"&
variable=="GuineaPig")
subsetters=testdat$DaysPostInfection==0|testdat$DaysPostInfection==3
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed"&
variable=="GuineaPig")
subsetters=testdat$DaysPostInfection==0|testdat$DaysPostInfection==4
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed"&
variable=="GuineaPig")
subsetters=testdat$DaysPostInfection==0|testdat$DaysPostInfection==5
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
testdat=subset(Tdat2, Vitality=="Alive" & Status=="Sugarfed"&
variable=="GuineaPig")
subsetters=testdat$DaysPostInfection==0|testdat$Treatment=="Hybrid"
testdat=testdat[subsetters,]
pairwise.t.test(testdat$value, testdat$Treatment,  p.adj="none")
#Plot host preference
theme = theme_bw()+theme(text = element_text(size=20), axis.title.x =
element_text(size=30), axis.text.x = element_text(hjust=1, vjust=.5,
size=25), axis.text.y = element_text(size=25), title =
element_text(size=35), legend.title = element_text(size=25), legend.text
= element_text(size=20),  legend.key = element_blank())
limits=aes(ymax=Mosquitoes+se, ymin=Mosquitoes-se)
attach(Tdat5)
cbPalette <- c("#007acc", "#cc0000", "#29a329")
t.plt1=ggplot(subset(Tdat5, variable=="GuineaPig"& Treatment!="No Host"),
aes(DaysPostInfection, Mosquitoes, color=Treatment))+
  geom_line(size=3)+geom_errorbar(limits, width=.2, size=2,
linetype="solid")+
  theme+scale_colour_manual(values=cbPalette)+
  xlab("Days post fungal infection")+
  ylab("% mosquitoes in chamber\nnear host")+ggtitle("Impact of fungal
infection on\nblood-feeding")+
  scale_y_continuous(labels=percent)
t.plt1 
#Combine host preference and mortality for transmission
Mort=rbind(subset(Mdat4, Day==1), subset(Mdat4, Day==2), subset(Mdat4,
Day==3), subset(Mdat4, Day==4), subset(Mdat4, Day==5))
Mort=Mort[1:4]
Trans=Tdat5
Trans$Day=as.numeric(Trans$DaysPostInfection)
Trans=subset(Trans, variable=="GuineaPig")
Trans=Trans[Trans$Treatment!="No Host", c(2,4,7)]
colnames(Trans)[2]="Interest"
dat.m=merge(Mort, Trans)
dat.m$Transmission=dat.m$mean*dat.m$Interest
#Plot transmission data
cbPalette <- c("#cc0000", "#29a329", "#007acc")
theme = theme_bw()+theme(text = element_text(size=20), axis.title.x =
element_text(size=30), axis.text.x = element_text(angle=90, hjust=1,
vjust=.5, size=25), axis.text.y = element_text(size=25), title =
element_text(size=35), legend.title = element_text(size=25), legend.text
= element_text(size=20))
t.plt2=ggplot(subset(dat.m, Species=="An.coluzzii"))+geom_area(aes(Day,
mean, color=Treatment,
      fill=Treatment, position="identity"),
alpha=0.4)+scale_fill_manual(values=cbPalette)+
      geom_area(aes(Day, Transmission, color=Treatment, fill=Treatment,
position="identity"),alpha=0.5)+
      geom_hline(yintercept=.2, size=1, color="white",
linetype="dashed")+
      geom_hline(yintercept=.5, size=1, color="white",
linetype="dashed")+
      theme+scale_colour_manual(values=cbPalette)+xlab("Days after
exposure")+ylab("Percent survival")+
      scale_y_continuous(breaks=pretty_breaks(n = 10), labels=percent)+
      scale_x_continuous(breaks=0:max(dat.m$Day))+facet_wrap(~Treatment)
t.plt2
####Two-Toxin Analysis####
#Reshape and reformat dataset
Txdat2=melt(Txdat, colnames(Txdat)[c(1:3)],
colnames(Txdat)[c(4:length(colnames(Txdat)))])
colnames(Txdat2)[4]="Day"
levels(Txdat2$Day)=c("1", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5",
"6", "7", "8")
Txdat2$Day=as.numeric(as.character(Txdat2$Day))
levels(Txdat2$Treatment)[c(1,4)]=c("AaIT", "Hybrid/AaIT")
colnames(Txdat2)[5]="Dead"
#Calculate number of alive mosquitoes
Txdat2$Alive=30-Txdat2$Dead
#Calculate LT50 for each replicate across treatments and concentration
attach(Txdat2)
surv.per=.5
Txdat2.LT=ddply(Txdat2,.(Treatment, Replicate, Concentration), summarize,
LT=as.numeric(dose.p(glm(cbind(Dead,Alive)~Day,binomial),p=surv.per)))
#Calculate SE and mean LT50 across treatments and concentration
Txdat2.LT50.Error=ddply(Txdat2.LT, .(Treatment, Concentration),
summarize, "LT50 Mean"=mean(LT,na.rm=T), se=sd(LT,
na.rm=T)/sqrt(length(LT)), Replicates=length(LT[!is.na(LT)]))
#Create csv with results
#write.csv(Txdat2.LT50.Error, file = "Toxin LT50.csv")
#Reformat dataset for LC50 calculations
Txdat2$Day=as.factor(Txdat2$Day)
Txdat2$Concentration=as.numeric(as.character(Txdat2$Concentration))
Txdat2$Concentration=10^Txdat2$Concentration
#Calculate LC50 for each replicate across Treatments for 5 days post
infection
attach(subset(Txdat2, Treatment!="Untreated"& Day %in%
c("4.5","5","5.5","6")))
surv.per=.5
Txdat2.LC=ddply(subset(Txdat2, Treatment!="Untreated" & Day == "5"),
.(Treatment, Replicate, Day), summarize,
LC=as.numeric(dose.p(glm(cbind(Dead,Alive)~Concentration,binomial),p=surv
.per)))
Txdat2.LC[Txdat2.LC$LC<0,]$LC=NA
#Calculate SE and mean LC50 across treatments for 5 days post infection
Txdat2.LC50.Error=ddply(Txdat2.LC, .(Treatment, Day), summarize, "LC50
Mean"=mean(LC,na.rm=T), se=sd(LC, na.rm=T)/sqrt(length(LC)),
Replicates=length(LC[!is.na(LC)]))
#Create csv with results
#write.csv(Txdat2.LC50.Error, file = "Toxin LC50.csv")
#Pull LT50 means for visualization
point.dat=Txdat2.LT50.Error
point.dat$Concentration=as.factor(point.dat$Concentration)
#Create dataset representing Fang et. al 2011 spore count data
Load.dat=data.frame(Concentration=as.factor(c(7,6,5)),value=c(10000000,
1000000,100000), Spores=c(201,41, 4))
#Represent control data in each concentration
point.dat.c=data.frame(Treatment=rep("Untreated", 3),
Concentration=c(5:7))
point.dat.c=merge(point.dat.c, point.dat[c(1,3:5)])
#Add control treatements for each concentration to original dataset
point.dat=rbind(subset(point.dat, Treatment!="Untreated"), point.dat.c)
#Merge spore count data  and LC50 data to LT50 means
point.dat=merge(point.dat, Load.dat[c(1,3)], all.x=T)
point.dat=merge(point.dat, Txdat2.LC50.Error[c(1,3)], all.x=T)
#Generate linear model to calculate inoculum based on concentration
spore.model=glm(Spores~value, dat=Load.dat)
#Estimate inoculum relative to LC50
point.dat$Est.LC50.Spores=predict(spore.model,newdata=data.frame(value=po
int.dat$`LC50 Mean`), type="response")
#Calculate dose (Inoculum/LC50)
point.dat$dose=point.dat$Spores/point.dat$Est.LC50.Spores
#Fill in zero for untreated mosquitoes
point.dat[point.dat$Treatment=="Untreated",]$`LC50 Mean`=0
point.dat[point.dat$Treatment=="Untreated",]$Est.LC50.Spores=0
point.dat[point.dat$Treatment=="Untreated",]$dose=0
#Reformat data for visualization
point.dat$Spores=as.factor(point.dat$Spores)
point.dat$Treatment <- factor(point.dat$Treatment, levels =
levels(point.dat$Treatment)[c(8,4,5,6,7,1,2,3)])
levels(point.dat$Spores)=c("~4 Spores", "41Â±6 Spores", "200Â±9 Spores")
#Plot LT50s by dose
limits=aes(ymax=`LT50 Mean`+se, ymin=`LT50 Mean`-se)
theme = theme_bw()+theme(text = element_text(size=25), axis.title.x =
element_text(size=30), axis.title.x = element_text(size=30), title =
element_text(size=30))
cbPalette <- c("#a812ce", "#ce8f12", "#54ce12", "#ce2e12", "#129cce")
ar.plt1=ggplot(point.dat, aes(dose, `LT50 Mean`, color=Treatment))+
  geom_point(stat="identity", position="dodge",
size=5)+geom_errorbar(limits, position=position_dodge(.9), width=.2,
size=2)+
  theme+xlab("Inoculum/LC50")+scale_colour_manual(values=cbPalette)+
  ylab("LT50 (Days)")+ggtitle("LT50 of Single, Dual Toxin Strains by
Inoculum/LC50")+facet_wrap(~Spores)
ar.plt1
#Test for synergism
#Calculate SE and mean LC50 across treatments for 5 days post infection
for toxin expressing strains
syn.tab=ddply(subset(ddply(subset(Txdat2, Treatment!="Untreated" & Day
%in% c(4.5,5)),
    .(Treatment, Replicate, Day), summarize,
LC=as.numeric(dose.p(glm(cbind(Dead,Alive)~Concentration,binomial),p=surv
.per))),
    Treatment!="WT"), .(Treatment, Day), summarize,
Mean=mean(LC,na.rm=T),
    se=sd(LC, na.rm=T)/sqrt(length(LC)))
#Reformat dataset for synergy calculations
syn.tab2=melt(syn.tab, colnames(syn.tab)[c(1:2)],
colnames(syn.tab)[c(3,4)])
syn.tab2=cast(syn.tab2, Day~Treatment+variable)
#Calculate expected value for combined toxins and ratio with observed
value
syn.tab2$Expected=(0.5/syn.tab2$AaIT_Mean+0.5/syn.tab2$Hybrid_Mean)^-1
syn.tab2$Ratio=syn.tab2$`Hybrid/AaIT_Mean`/syn.tab2$Expected
#Reformat table
syn.tab3=data.frame(Day=syn.tab$Day, Treatment=syn.tab$Treatment,
LC50=paste(format(syn.tab$Mean, scientific=T, digits=3),
format(syn.tab$se, scientific=T, digits=3), sep="Â±"))
syn.tab3=subset(syn.tab3, Treatment=="Hybrid/AaIT")
colnames(syn.tab3)[3]="LC50 (MeanÂ±SE)"
syn.tab3$Expected=NA
syn.tab3$`Ratio (Obs/Exp)`=NA
syn.tab3$Expected=format(syn.tab2$Expected, scientific=T, digits=3)
syn.tab3$`Ratio (Obs/Exp)`=syn.tab2$Ratio
#Create csv with results
#write.csv(syn.tab3, file = "Synergism table.csv")
####Toxin Arsenal Analysis####
#Reformat dataset
TxArdat2=melt(TxArdat, colnames(TxArdat)[1:2],
colnames(TxArdat)[3:length(colnames(TxArdat))])
levels(TxArdat2$variable)[1:11]=1:11
#Pull n values and add as column
TxArdat2.n=subset(TxArdat2, variable=="n")
TxArdat2$variable=as.numeric(as.character(TxArdat2$variable))
colnames(TxArdat2.n)[4]="n"
TxArdat2=subset(TxArdat2, variable!="n")
TxArdat2=merge(TxArdat2, TxArdat2.n[c(1:2, 4)])
colnames(TxArdat2)[3:4]=c("Day", "Dead")
#Calculate alive mosquitoes for LT50
TxArdat2$Alive=TxArdat2$n-TxArdat2$Dead
#Calculate LT50 for each replicate of each toxin
attach(TxArdat2)
surv.per=.5
TxArdat2.LT=ddply(TxArdat2,.(Treatment, Replicate), summarize,
LT=as.numeric(dose.p(glm(cbind(Dead,Alive)~Day,binomial),p=surv.per)))
#Calculate SE and mean LT50 for each toxin
TxArdat2.LT50.Error=ddply(TxArdat2.LT, .(Treatment), summarize, "LT50
Mean"=mean(LT,na.rm=T), se=sd(LT, na.rm=T)/sqrt(length(LT)),
Replicates=length(LT[!is.na(LT)]))
#Pull 1x10^6 concentrations from two-toxin bioassays for comparison
Master.Tx.LT50.Error=rbind(TxArdat2.LT50.Error, subset(Txdat2.LT50.Error,
Concentration==6)[c(1,3:5)])
#Pull LT values for each replicate and statistics
Master.Tx.LT=TxArdat2.LT[1:2]
#Standardizing dataset format
Master.Tx.LT$Concentration="6"
Master.Tx.LT$LT=TxArdat2.LT$LT
#Combine LT50s of dual toxin and toxin arsenal analyses
Master.Tx.LT=rbind(Master.Tx.LT, Txdat2.LT)
#Run t test on LT50s
subsetter=Master.Tx.LT$Concentration=="5"|Master.Tx.LT$Concentration==0
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Treatment, p.adj="none")
subsetter=Master.Tx.LT$Concentration=="6"|Master.Tx.LT$Concentration==0
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Treatment, p.adj="none")
subsetter=Master.Tx.LT$Concentration=="7"|Master.Tx.LT$Concentration==0
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Treatment, p.adj="none")
subsetter=Master.Tx.LT$Treatment=="WT"
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Concentration, p.adj="none")
subsetter=Master.Tx.LT$Treatment=="AaIT"
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Concentration, p.adj="none")
subsetter=Master.Tx.LT$Treatment=="Hybrid"
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Concentration, p.adj="none")
subsetter=Master.Tx.LT$Treatment=="Hybrid/AaIT"
testdat=Master.Tx.LT[subsetter,]
pairwise.t.test(testdat$LT,testdat$Concentration, p.adj="none")
#Create csv with results
#write.csv(Master.Tx.LT50.Error, file = "Toxin Arsenal LT50s.csv")
####Spore Count Analysis####
#Reformat spore count data
Sdat=Sdat[c(2,4,5)]
levels(Sdat$Replicate)=1:4
#Calculate the mean across replicates
Sdat2=ddply(Sdat, .(Replicate), summarize, mn=mean(spores.mosq),
count=length(spores.mosq))
#Calculate the mean among replicates
Sdat3=ddply(subset(Sdat2, Replicate!=3), .(), summarize, mn.all=mean(mn),
sd=sd(mn))
####Log-Rank Testing####
#install.packages("survival")
library(survival)
#Individualize dual toxin data for Log-Rank Tests
Txdat3=Txdat
Txdat3$Alive=30
for (i in 16:5){
  Txdat3[i]=Txdat3[i]-Txdat3[i-1]
}
Txdat4=melt(Txdat3, colnames(Txdat3)[1:3], colnames(Txdat3)[4:16])
colnames(Txdat4)[4]="Day"
levels(Txdat4$Day)=c("1", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5",
"6", "7", "8", "Alive")
levels(Txdat4$Treatment)[c(1,4)]=c("AaIT", "Hybrid/AaIT")
Treatment=as.vector(rep(Txdat4$Treatment, Txdat4$value))
Concentration=as.vector(rep(Txdat4$Concentration, Txdat4$value))
Replicate=as.vector(rep(Txdat4$Replicate, Txdat4$value))
Day=as.vector(rep(Txdat4$Day, Txdat4$value))
LogR.dat=data.frame(Treatment, Concentration, Replicate, Day)
LogR.dat$status=0
LogR.dat[Day!="Alive",]$status=1
LogR.dat[Day=="Alive",]$Day="8"
LogR.dat$Day=as.numeric(as.character(LogR.dat$Day)



