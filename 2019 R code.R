#bring in datasets

bowls <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/Insect%20ID%202019%20-%20Bowl_natural.csv?token=GHSAT0AAAAAABSFZTB2FZOSER2YXMUQZEDAYS5UCHQ",na.strings = NULL)
summary(bowls)
str(bowls) 
ramps <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/Insect%20ID%202019%20-%20Ramp_natural.csv?token=GHSAT0AAAAAABSFZTB2DXWH7WA4B2PQ5HFOYS5UC6Q",na.strings = NULL)
summary(ramps)
str(ramps)
sticky <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/Insect%20ID%202019%20-%20Sticky%20card_natural.csv?token=GHSAT0AAAAAABSFZTB32TJ5QMTROVWHS32SYS5UDJQ",na.strings = NULL)
summary(sticky)
str(sticky)

#add trap type as a column on each data file
bowls$Trap="bowl"
ramps$Trap="ramp"
sticky$Trap="sticky"

#combine data tables 
library (plyr)
bowlramp <- rbind.fill (bowls, ramps)
allbugs <-rbind.fill (bowlramp, sticky)

#add column for region (south central north)
allbugs$region<-ifelse(allbugs$Site=="WLR", "South",
                       ifelse(allbugs$Site=="WPR", "South",
                              ifelse(allbugs$Site=="SNY", "South",
                                     ifelse(allbugs$Site=="DAL", "North", 
                                            ifelse(allbugs$Site=="BAL", "North",
                                                   ifelse(allbugs$Site == "CHA", "North", "Central")
                                            )))))
str(allbugs)

#To obtain richness counts
allbugs.rowsums <- rowSums(allbugs[,4:48]>0)
allbugs$richness <- allbugs.rowsums

#To obtain abundance counts
allbugs.abun <- rowSums(allbugs[,4:48])
allbugs$abundance <- allbugs.abun

#load vegan
library(vegan)

#calculate Shannon diversity
diversity <-diversity(allbugs[,4:48])
allbugs$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(allbugs[,4:48]))
allbugs$evenness <- evenness

#look at data set
summary(allbugs)
str(allbugs)


##richness linear mixed effects model
library (emmeans) #for pairwise comparisons
library(lme4)
library(lmerTest) #to obtain p values

richmodel <- lmer(richness~Site+region+(1|Date)+(1|Trap),data=allbugs)  #AIC = 1057
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(richmodel)
AIC(richmodel)
anova(richmodel)

rich.emm<-emmeans(richmodel,pairwise~region) #comparing region richness
rich.emm
#results: same for all
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld 

rich.emm.s<-emmeans(richmodel,pairwise~Site) #comparing site richness
rich.emm.s
#results: differences between BFB-DGM, DGM-SSH, WLR-WPR
rich.cld.s<-multcomp::cld(rich.emm.s, alpha = 0.05, Letters = LETTERS)
rich.cld.s

#

##abundance linear mixed effects model
abunmodel <- lmer(abundance~Site+region+(1|Date)+(1|Trap),data=allbugs)  #AIC = 3152
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(abunmodel)
AIC(abunmodel)
anova(abunmodel)

abun.emm<-emmeans(abunmodel,pairwise~region) #comparing region abundance
abun.emm
#results: same for all
abun.cld<-multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS)
abun.cld 

abun.emm.s<-emmeans(abunmodel,pairwise~Site) #comparing site abundance
abun.emm.s
#results: differences btw BAL-DAL, CHA-DAL
abun.cld.s<-multcomp::cld(abun.emm.s, alpha = 0.05, Letters = LETTERS)
abun.cld.s

#

##diversity linear mixed effects model
divmodel <- lmer(diversity~Site+region+(1|Date)+(1|Trap),data=allbugs)  #AIC = 291
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(divmodel)
AIC(divmodel)
anova(divmodel)

div.emm<-emmeans(divmodel,pairwise~region) #comparing region diversity
div.emm
#results: same for all
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld 

div.emm.s<-emmeans(divmodel,pairwise~Site) #comparing site diversity
div.emm.s
#results: no differences between sites
div.cld.s<-multcomp::cld(div.emm.s, alpha = 0.05, Letters = LETTERS)
div.cld.s

#

##evenness linear mixed effects model
evenmodel <- lmer(evenness~Site+region+(1|Date)+(1|Trap),data=allbugs)  #AIC = -92
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(evenmodel)
AIC(evenmodel)
anova(evenmodel)

even.emm<-emmeans(evenmodel,pairwise~region) #comparing region evenness
even.emm
#results: same for all
even.cld<-multcomp::cld(even.emm, alpha = 0.05, Letters = LETTERS)
even.cld 

even.emm.s<-emmeans(evenmodel,pairwise~Site) #comparing site richness
even.emm.s
#results: no differences between sites
even.cld.s<-multcomp::cld(even.emm.s, alpha = 0.05, Letters = LETTERS)
even.cld.s

#######
#ggplot box plots
library (ggplot2)

#site richness by region
richness.plot<-ggplot(allbugs, aes(x = factor(region,level = c("South","Central","North")), y = richness, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position="bottom")+
  labs(title="", x="", y="Richness")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
richness.plot

#site abundance by region
abundance.plot<-ggplot(allbugs, aes(x = factor(region, level = c("South","Central","North")), y = abundance, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="bottom")+
  labs(title="", x="", y="Abundance (log10)")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_y_continuous(trans="log10")+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
abundance.plot

#site diversity by region
diversity.plot<-ggplot(allbugs, aes(x = factor(region,level = c("South","Central","North")), y = diversity, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position="bottom")+
  labs(title="", x="", y="Shannon diversity")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
diversity.plot

#site evenness by region
evenness.plot<-ggplot(allbugs, aes(x = factor(region,level = c("South","Central","North")), y = evenness, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position="bottom")+
  labs(title="", x="Sites by Region", y="Evenness")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
evenness.plot

###
#mush together plots
library(ggpubr) 
boxplot <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                     #labels = c("A", "B", "C", "D"),
                     ncol = 1, nrow = 4,
                     common.legend = TRUE, legend = "bottom")
boxplot

pdf("boxplot.pdf", height=8, width=8) #height and width in inches
boxplot
dev.off()
###

#find colors in "Paired" color palette
library(RColorBrewer)
display.brewer.all(colorblindFriendly = TRUE)
display.brewer.pal(n = 9, name = 'Paired')
brewer.pal(n = 9, name = "Paired") #only has 8 colors, but we have 9 sites
#"#A6CEE3" "#1F78B4" "#B2DF8A" "#33A02C" "#FB9A99" "#E31A1C" "#FDBF6F" "#FF7F00" "#CAB2D6"
#######

#NMDS of insect community 
library (vegan)

#Create matrix of environmental variables
env.matrix<-allbugs[c(1:3,51:52)]
#create matrix of community variables
com.matrix<-allbugs[c(4:50)]

#ordination by NMDS
NMDS<-metaMDS(com.matrix, distance="bray", k=2, autotransform=FALSE, trymax=300)
###stress = no convergence (but once it said 0.15)
stressplot(NMDS)

#plot NMDS for region
#might need to change colors
#8 x 10
plot(NMDS, disp='sites', type="n")
#title(main="Arthropod community composition by region", cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS, env.matrix$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
ordiellipse(NMDS, env.matrix$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South") 
#add data points
points(NMDS, display="sites", select=which(env.matrix$region=="North"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(env.matrix$region=="Central"), pch=17, col="#009E73")
points(NMDS, display="sites", select=which(env.matrix$region=="South"), pch=15, col="#CC79A7")
#add legend
legend(2.62,2.95, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=1.5, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (regions)
fit<-adonis(com.matrix ~ region, data = env.matrix, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$region))
#P-value = 0.6367 -- assumes homogeneity of multivariate dispersion

###

#plot NMDS for Northern sites 
#8 x 10
plot(NMDS, disp='sites', type="n")
title(main="North", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
#add data points
points(NMDS, display="sites", select=which(env.matrix$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS, display="sites", select=which(env.matrix$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS, display="sites", select=which(env.matrix$Site=="DAL"), pch=15, col="#33A02C")
#add legend
legend(1.43,2.95, title=NULL, pch=c(19,17,15), col=c("#A6CEE3","#B2DF8A","#33A02C"), cex=1.55, legend=c("Beaton alvar", "Cape Hurd alvar", "Davis alvar"))

#plot NMDS for Central sites 
#8 x 10
plot(NMDS, disp='sites', type="n")
title(main="Central", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
#add data points
points(NMDS, display="sites", select=which(env.matrix$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS, display="sites", select=which(env.matrix$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS, display="sites", select=which(env.matrix$Site=="BFB"), pch=15, col="#1F78B4")
#add legend
legend(0.15,2.95, title=NULL, pch=c(19,17,15), col=c("#FDBF6F","#FB9A99","#1F78B4"), cex=1.55, legend=c("Slate shale hill", "Dusty Goldenrod meadow", "Bedford barren"))

#plot NMDS for Southern sites 
#8 x 10
plot(NMDS, disp='sites', type="n")
title(main="South", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#CAB2D6",kind="sd", conf=0.95, label=FALSE, show.groups = "WPR")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#E31A1C",kind="sd", conf=0.95, label=FALSE, show.groups = "SNY")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FF7F00",kind="sd", conf=0.95, label=FALSE, show.groups = "WLR")
#add data points
points(NMDS, display="sites", select=which(env.matrix$Site=="WPR"), pch=17, col="#CAB2D6")
points(NMDS, display="sites", select=which(env.matrix$Site=="SNY"), pch=15, col="#E31A1C")
points(NMDS, display="sites", select=which(env.matrix$Site=="WLR"),pch=19, col="#FF7F00")
#add legend
legend(1.67,2.95, title=NULL, pch=c(19,17,15), col=c("#FF7F00","#CAB2D6","#E31A1C"), cex=1.55, legend=c("W ladder", "W picnic rock", "Synder hollow"))

#bootstrapping and testing for differences between the groups (sites)
fit<-adonis(com.matrix ~ Site, data = env.matrix, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Site))
#P-value = 0.1694 -- assumes homogeneity of multivariate dispersion

####

#do all of above again, but with data organized by ORDER

bowls_order <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/Insect%20ID%202019%20-%20Bowl_natural_order.csv?token=GHSAT0AAAAAABSFZTB2I25KABWPGWQIGUAAYS5UPCA",na.strings = NULL)
 
ramps_order <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/Insect%20ID%202019%20-%20Ramp_natural_order.csv?token=GHSAT0AAAAAABSFZTB2DWOHVZEZBSMF3FYKYS5UPOQ",na.strings = NULL)

sticky_order <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/Insect%20ID%202019%20-%20Sticky%20card_natural_order.csv?token=GHSAT0AAAAAABSFZTB3QJ7FKN3AJGOCONLAYS5UPYQ",na.strings = NULL)


#add trap type as a column on each data file
bowls_order$Trap="bowl"
ramps_order$Trap="ramp"
sticky_order$Trap="sticky"

#combine data tables 
library (plyr)
bowlramp_order <- rbind.fill (bowls_order, ramps_order)
allbugs_order <-rbind.fill (bowlramp_order, sticky_order)

#add column for region (south central north)
allbugs_order$region<-ifelse(allbugs_order$Site=="WLR", "South",
                       ifelse(allbugs_order$Site=="WPR", "South",
                              ifelse(allbugs_order$Site=="SNY", "South",
                                     ifelse(allbugs_order$Site=="DAL", "North", 
                                            ifelse(allbugs_order$Site=="BAL", "North",
                                                   ifelse(allbugs_order$Site == "CHA", "North", "Central")
                                            )))))
str(allbugs_order)

#To obtain richness counts
allbugs_order.rowsums <- rowSums(allbugs_order[,4:17]>0)
allbugs_order$richness <- allbugs_order.rowsums

#To obtain abundance counts
allbugs_order.abun <- rowSums(allbugs_order[,4:17])
allbugs_order$abundance <- allbugs_order.abun

#load vegan
library(vegan)

#calculate Shannon diversity
diversity_order <-diversity(allbugs_order[,4:17])
allbugs_order$diversity <-diversity_order

#calculate Evenness
evenness_order <-diversity_order/log(specnumber(allbugs_order[,4:17]))
allbugs_order$evenness <- evenness_order

#look at data set
summary(allbugs_order)
str(allbugs_order)


##richness linear mixed effects model
library (emmeans) #for pairwise comparisons
library(lme4)
library(lmerTest) #to obtain p values

richmodel_order <- lmer(richness~Site+region+(1|Date)+(1|Trap),data=allbugs_order)  #AIC = 766
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(richmodel_order)
AIC(richmodel_order)
anova(richmodel_order)

rich.emm_order<-emmeans(richmodel_order,pairwise~region) #comparing region richness
rich.emm_order
#results: same for all
rich.cld_order<-multcomp::cld(rich.emm_order, alpha = 0.05, Letters = LETTERS)
rich.cld_order 

rich.emm.s_order<-emmeans(richmodel_order,pairwise~Site) #comparing site richness
rich.emm.s_order
#results: differences between WLR-WPR
rich.cld.s_order<-multcomp::cld(rich.emm.s_order, alpha = 0.05, Letters = LETTERS)
rich.cld.s_order

#

##abundance linear mixed effects model
abunmodel_order <- lmer(abundance~Site+region+(1|Date)+(1|Trap),data=allbugs_order)  #AIC = 3153
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(abunmodel_order)
AIC(abunmodel_order)
anova(abunmodel_order)

abun.emm_order<-emmeans(abunmodel_order,pairwise~region) #comparing region abundance
abun.emm_order
#results: same for all
abun.cld_order<-multcomp::cld(abun.emm_order, alpha = 0.05, Letters = LETTERS)
abun.cld_order

abun.emm.s_order<-emmeans(abunmodel_order,pairwise~Site) #comparing site abundance
abun.emm.s_order
#results: differences btw BAL-DAL, CHA-DAL
abun.cld.s_order<-multcomp::cld(abun.emm.s_order, alpha = 0.05, Letters = LETTERS)
abun.cld.s_order

#

##diversity linear mixed effects model
divmodel_order <- lmer(diversity~Site+region+(1|Date)+(1|Trap),data=allbugs_order)  #AIC = 183
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(divmodel_order)
AIC(divmodel_order)
anova(divmodel_order)

div.emm_order<-emmeans(divmodel_order,pairwise~region) #comparing region diversity
div.emm_order
#results: same for all
div.cld_order<-multcomp::cld(div.emm_order, alpha = 0.05, Letters = LETTERS)
div.cld_order 

div.emm.s_order<-emmeans(divmodel_order,pairwise~Site) #comparing site diversity
div.emm.s_order
#results: differences between CHA-DAL
div.cld.s_order<-multcomp::cld(div.emm.s_order, alpha = 0.05, Letters = LETTERS)
div.cld.s_order

#

##evenness linear mixed effects model
evenmodel_order <- lmer(evenness~Site+region+(1|Date)+(1|Trap),data=allbugs_order)  #AIC = -39
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(evenmodel_order)
AIC(evenmodel_order)
anova(evenmodel_order)

even.emm_order<-emmeans(evenmodel_order,pairwise~region) #comparing region evenness
even.emm_order
#results: same for all
even.cld_order<-multcomp::cld(even.emm_order, alpha = 0.05, Letters = LETTERS)
even.cld_order 

even.emm.s_order<-emmeans(evenmodel_order,pairwise~Site) #comparing site richness
even.emm.s_order
#results: no differences between sites
even.cld.s_order<-multcomp::cld(even.emm.s_order, alpha = 0.05, Letters = LETTERS)
even.cld.s_order

#######
#ggplot box plots
library (ggplot2)

#site richness by region
richness.plot_order<-ggplot(allbugs_order, aes(x = factor(region,level = c("South","Central","North")), y = richness, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position="bottom")+
  labs(title="", x="", y="Richness")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
richness.plot

#site abundance by region
abundance.plot_order<-ggplot(allbugs_order, aes(x = factor(region, level = c("South","Central","North")), y = abundance, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="bottom")+
  labs(title="", x="", y="Abundance (log10)")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_y_continuous(trans="log10")+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
abundance.plot

#site diversity by region
diversity.plot_order<-ggplot(allbugs_order, aes(x = factor(region,level = c("South","Central","North")), y = diversity, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position="bottom")+
  labs(title="", x="", y="Shannon diversity")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
diversity.plot

#site evenness by region
evenness.plot_order<-ggplot(allbugs_order, aes(x = factor(region,level = c("South","Central","North")), y = evenness, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position="bottom")+
  labs(title="", x="Sites by Region", y="Evenness")+
  #theme (plot.title = element_text(hjust=0.5))+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
evenness.plot

###
#mush together plots
library(ggpubr) 
boxplot_order <- ggarrange(richness.plot_order, abundance.plot_order, diversity.plot_order, evenness.plot_order,
                     #labels = c("A", "B", "C", "D"),
                     ncol = 1, nrow = 4,
                     common.legend = TRUE, legend = "bottom")
boxplot_order

pdf("boxplot_order.pdf", height=8, width=8) #height and width in inches
boxplot_order
dev.off()

#NMDS of insect community 
library (vegan)

#Create matrix of environmental variables
env.matrix_order<-allbugs_order[c(1:3,18:19)]
#create matrix of community variables
com.matrix_order<-allbugs_order[c(4:17)]

#ordination by NMDS
NMDS_order<-metaMDS(com.matrix_order, distance="bray", k=2, autotransform=FALSE, trymax=300)
###stress = 0.12
stressplot(NMDS_order)

#plot NMDS for region
#8 x 10
plot(NMDS_order, disp='sites', type="n")
#title(main="Arthropod community composition by region", cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS_order, env.matrix_order$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
ordiellipse(NMDS_order, env.matrix_order$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South") 
#add data points
points(NMDS_order, display="sites", select=which(env.matrix_order$region=="North"),pch=19, col="#E69F00")
points(NMDS_order, display="sites", select=which(env.matrix_order$region=="Central"), pch=17, col="#009E73")
points(NMDS_order, display="sites", select=which(env.matrix_order$region=="South"), pch=15, col="#CC79A7")
#add legend
legend(1.62,2.95, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=1.5, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (regions)
fit<-adonis(com.matrix_order ~ region, data = env.matrix_order, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_order)
anova(betadisper(distances_data, env.matrix_order$region))
#P-value = 0.3606 -- assumes homogeneity of multivariate dispersion

###

#plot NMDS for Northern sites 
#8 x 10
plot(NMDS_order, disp='sites', type="n")
title(main="North", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
#add data points
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="DAL"), pch=15, col="#33A02C")
#add legend
legend(0.55,2.95, title=NULL, pch=c(19,17,15), col=c("#A6CEE3","#B2DF8A","#33A02C"), cex=1.55, legend=c("Beaton alvar", "Cape Hurd alvar", "Davis alvar"))

#plot NMDS for Central sites 
#8 x 10
plot(NMDS_order, disp='sites', type="n")
title(main="Central", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
#add data points
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="BFB"), pch=15, col="#1F78B4")
#add legend
legend(-0.59,2.95, title=NULL, pch=c(19,17,15), col=c("#FDBF6F","#FB9A99","#1F78B4"), cex=1.55, legend=c("Slate shale hill", "Dusty Goldenrod meadow", "Bedford barren"))

#plot NMDS for Southern sites 
#8 x 10
plot(NMDS_order, disp='sites', type="n")
title(main="South", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#CAB2D6",kind="sd", conf=0.95, label=FALSE, show.groups = "WPR")
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#E31A1C",kind="sd", conf=0.95, label=FALSE, show.groups = "SNY")
ordiellipse(NMDS_order, env.matrix_order$Site, draw="polygon", col="#FF7F00",kind="sd", conf=0.95, label=FALSE, show.groups = "WLR")
#add data points
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="WPR"), pch=17, col="#CAB2D6")
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="SNY"), pch=15, col="#E31A1C")
points(NMDS_order, display="sites", select=which(env.matrix_order$Site=="WLR"),pch=19, col="#FF7F00")
#add legend
legend(0.77,2.95, title=NULL, pch=c(19,17,15), col=c("#FF7F00","#CAB2D6","#E31A1C"), cex=1.55, legend=c("W ladder", "W picnic rock", "Synder hollow"))

#bootstrapping and testing for differences between the groups (sites)
fit<-adonis(com.matrix ~ Site, data = env.matrix, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Site))
#P-value = 0.1694 -- assumes homogeneity of multivariate dispersion


#########

#Plant analysis

#bring in plant data
plants <- read.csv("https://raw.githubusercontent.com/BahlaiLab/Manning_K/master/2019%20Built%20by%20Nature/2019%20insect%20ID%20analysis/2019%20Insect%20ID/2019%20plants_p.a..csv?token=GHSAT0AAAAAABSFZTB2OJC75GU3NE3AGI2WYS7CP3A",na.strings = NULL)
summary(plants)









###################################################################################

#NOT USING ANYTHING AFTER HERE

#plot NMDS for traps
plot(NMDS, disp='sites', type="n")
#title(main="Arthropod community composition by trap type", cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "bowl")
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "ramp") 
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "sticky") 
#add data points
points(NMDS, display="sites", select=which(env.matrix$Trap=="bowl"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(env.matrix$Trap=="ramp"), pch=17, col="#009E73")
points(NMDS, display="sites", select=which(env.matrix$Trap=="sticky"), pch=15, col="#CC79A7")
#add legend
legend(3.5,2.6, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=.8, legend=c("Bowl", "Ramp", "Sticky card"))

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix ~ Trap, data = env.matrix, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Trap))
#P-value = 0.005 -- CANNOT assume homogeneity of multivariate dispersion
##############

#NMDS for trap type data

#first let's prep the data

#add column for region (south central north)
bowls$region<-ifelse(bowls$Site=="WLR", "South",
                       ifelse(bowls$Site=="WPR", "South",
                              ifelse(bowls$Site=="SNY", "South",
                                     ifelse(bowls$Site=="DAL", "North", 
                                            ifelse(bowls$Site=="BAL", "North",
                                                   ifelse(bowls$Site == "CHA", "North", "Central")
                                            )))))
str(bowls)
ramps$region<-ifelse(ramps$Site=="WLR", "South",
                     ifelse(ramps$Site=="WPR", "South",
                            ifelse(ramps$Site=="SNY", "South",
                                   ifelse(ramps$Site=="DAL", "North", 
                                          ifelse(ramps$Site=="BAL", "North",
                                                 ifelse(ramps$Site == "CHA", "North", "Central")
                                          )))))
str(ramps)
sticky$region<-ifelse(sticky$Site=="WLR", "South",
                     ifelse(sticky$Site=="WPR", "South",
                            ifelse(sticky$Site=="SNY", "South",
                                   ifelse(sticky$Site=="DAL", "North", 
                                          ifelse(sticky$Site=="BAL", "North",
                                                 ifelse(sticky$Site == "CHA", "North", "Central")
                                          )))))
str(sticky)

#Create matrix of environmental variables
env.matrix_b<-bowls[c(1:3,64:65)]
env.matrix_r<-ramps[c(1:3,64:65)]
env.matrix_s<-sticky[c(1:3,64:65)]
#create matrix of community variables
com.matrix_b<-bowls[c(4:63)]
com.matrix_r<-ramps[c(4:63)]
com.matrix_s<-sticky[c(4:63)]

#ordination by NMDS for bowls
NMDS_b<-metaMDS(com.matrix_b, distance="bray", k=2, autotransform=FALSE, trymax=100)
###stress = .18
stressplot(NMDS_b)

#ordination by NMDS for ramps
NMDS_r<-metaMDS(com.matrix_r, distance="bray", k=2, autotransform=FALSE, trymax=100)
###stress = .17
stressplot(NMDS_r)

#ordination by NMDS for sticky cards
NMDS_s<-metaMDS(com.matrix_s, distance="bray", k=2, autotransform=FALSE, trymax=100)
###stress = .13
stressplot(NMDS_s)

#Now plot

#plot bowl NMDS for region
plot(NMDS_b, disp='sites', type="n")
title(main="Bowl", adj = 0.05, line = -2, cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_b, env.matrix_b$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS_b, env.matrix_b$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
ordiellipse(NMDS_b, env.matrix_b$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South") 
#add data points
points(NMDS_b, display="sites", select=which(env.matrix_b$region=="North"),pch=19, col="#E69F00")
points(NMDS_b, display="sites", select=which(env.matrix_b$region=="Central"), pch=17, col="#009E73")
points(NMDS_b, display="sites", select=which(env.matrix_b$region=="South"), pch=15, col="#CC79A7")
#add legend
legend(1.5,1.5, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=.8, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (sites - bowl)
fit<-adonis(com.matrix_b ~ region, data = env.matrix_b, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_b)
anova(betadisper(distances_data, env.matrix_b$region))
#P-value = 0.028 -- CANNOT assume homogeneity of multivariate dispersion

#

#plot ramp NMDS for region
plot(NMDS_r, disp='sites', type="n")
title(main="Ramp", adj = 0.05, line = -2, cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_r, env.matrix_r$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS_r, env.matrix_r$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
ordiellipse(NMDS_r, env.matrix_r$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South") 
#add data points
points(NMDS_r, display="sites", select=which(env.matrix_r$region=="North"),pch=19, col="#E69F00")
points(NMDS_r, display="sites", select=which(env.matrix_r$region=="Central"), pch=17, col="#009E73")
points(NMDS_r, display="sites", select=which(env.matrix_r$region=="South"), pch=15, col="#CC79A7")
#add legend
legend(2,1, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=.8, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (sites - ramp)
fit<-adonis(com.matrix_r ~ region, data = env.matrix_r, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_r)
anova(betadisper(distances_data, env.matrix_r$region))
#P-value = 0.0017 -- CANNOT assume homogeneity of multivariate dispersion

#

#plot sticky card NMDS for region
plot(NMDS_s, disp='sites', type="n")
title(main="Sticky card", adj = 0.05, line = -2, cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_s, env.matrix_s$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS_s, env.matrix_s$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
ordiellipse(NMDS_s, env.matrix_s$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South") 
#add data points
points(NMDS_s, display="sites", select=which(env.matrix_s$region=="North"),pch=19, col="#E69F00")
points(NMDS_s, display="sites", select=which(env.matrix_s$region=="Central"), pch=17, col="#009E73")
points(NMDS_s, display="sites", select=which(env.matrix_s$region=="South"), pch=15, col="#CC79A7")
#add legend
legend(2,1, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=.8, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (sites - ramp)
fit<-adonis(com.matrix_s ~ region, data = env.matrix_s, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_s)
anova(betadisper(distances_data, env.matrix_s$region))
#P-value = 0.32 -- assumes homogeneity of multivariate dispersion

###

#plot bowl NMDS by sites 
plot(NMDS_b, disp='sites', type="n")
title(main="", adj = 0.05, line = -2, cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#CAB2D6",kind="sd", conf=0.95, label=FALSE, show.groups = "WPR")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#E31A1C",kind="sd", conf=0.95, label=FALSE, show.groups = "SNY")
ordiellipse(NMDS_b, env.matrix_b$Site, draw="polygon", col="#FF7F00",kind="sd", conf=0.95, label=FALSE, show.groups = "WLR")
#add data points
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="DAL"), pch=15, col="#33A02C")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="BFB"), pch=15, col="#1F78B4")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="WPR"), pch=17, col="#CAB2D6")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="SNY"), pch=15, col="#E31A1C")
points(NMDS_b, display="sites", select=which(env.matrix_b$Site=="WLR"),pch=19, col="#FF7F00")

#bootstrapping and testing for differences between the groups (sites - bowl)
fit<-adonis(com.matrix_b ~ Site, data = env.matrix_b, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_b)
anova(betadisper(distances_data, env.matrix_b$Site))
#P-value = 0.3171 -- assumes homogeneity of multivariate dispersion

#

#plot ramp NMDS by sites 
plot(NMDS_r, disp='sites', type="n")
title(main="", adj = 0.05, line = -2, cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#CAB2D6",kind="sd", conf=0.95, label=FALSE, show.groups = "WPR")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#E31A1C",kind="sd", conf=0.95, label=FALSE, show.groups = "SNY")
ordiellipse(NMDS_r, env.matrix_r$Site, draw="polygon", col="#FF7F00",kind="sd", conf=0.95, label=FALSE, show.groups = "WLR")
#add data points
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="DAL"), pch=15, col="#33A02C")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="BFB"), pch=15, col="#1F78B4")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="WPR"), pch=17, col="#CAB2D6")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="SNY"), pch=15, col="#E31A1C")
points(NMDS_r, display="sites", select=which(env.matrix_r$Site=="WLR"),pch=19, col="#FF7F00")

#bootstrapping and testing for differences between the groups (sites - ramp)
fit<-adonis(com.matrix_r ~ Site, data = env.matrix_r, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_r)
anova(betadisper(distances_data, env.matrix_r$Site))
#P-value = 0.0086 -- CANNOT assume homogeneity of multivariate dispersion

#

#plot sticky card NMDS by sites 
plot(NMDS_s, disp='sites', type="n")
title(main="", adj = 0.05, line = -2, cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#CAB2D6",kind="sd", conf=0.95, label=FALSE, show.groups = "WPR")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#E31A1C",kind="sd", conf=0.95, label=FALSE, show.groups = "SNY")
ordiellipse(NMDS_s, env.matrix_s$Site, draw="polygon", col="#FF7F00",kind="sd", conf=0.95, label=FALSE, show.groups = "WLR")
#add data points
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="DAL"), pch=15, col="#33A02C")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="BFB"), pch=15, col="#1F78B4")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="WPR"), pch=17, col="#CAB2D6")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="SNY"), pch=15, col="#E31A1C")
points(NMDS_s, display="sites", select=which(env.matrix_s$Site=="WLR"),pch=19, col="#FF7F00")

#bootstrapping and testing for differences between the groups (sites - sticky cards)
fit<-adonis(com.matrix_s ~ Site, data = env.matrix_s, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_s)
anova(betadisper(distances_data, env.matrix_s$Site))
#P-value = 0.04 -- CANNOT assume homogeneity of multivariate dispersion
#########

#LMs by trap type

#bowl
#To obtain richness counts
bowls.rowsums <- rowSums(bowls[,4:63]>0)
bowls$richness <- bowls.rowsums
summary(bowls)
str(bowls)

richmodel_b <- lm(richness~Site+Data+region,data=bowls)  #AIC = 262.09
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(richmodel_b)
anova(richmodel_b)
AIC(richmodel_b)
emmeans(richmodel_b,pairwise~region) #comparing region richness
#results: sig diff btw Central-South
emmeans(richmodel_b,pairwise~Site) #comparing site richness
#results: sig differences btw BFB-DAL, BFB-WLR

#ramp
#To obtain richness counts
ramps.rowsums <- rowSums(ramps[,4:63]>0)
ramps$richness <- ramps.rowsums
summary(ramps)
str(ramps)

richmodel_r <- lm(richness~Site+Date+region,data=ramps)  #AIC = 311.33
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(richmodel_r)
anova(richmodel_r)
AIC(richmodel_r)
emmeans(richmodel_r,pairwise~region) #comparing region richness
#results: no sig differences
emmeans(richmodel_r,pairwise~Site) #comparing site richness
#results: sig differences - BFB-DGM, BFB-BAL, BFB-SNY, BFB-WLR

#sticky cards
#To obtain richness counts
sticky.rowsums <- rowSums(sticky[,4:63]>0)
sticky$richness <- sticky.rowsums
summary(sticky)
str(sticky)

richmodel_s <- lm(richness~Site+Date+region,data=sticky)  #AIC = 517.82
#region doesn't do anything in GLM, but you need it in to get values for site comparisons (and then can also get region comparisons)
summary(richmodel_s)
anova(richmodel_s)
AIC(richmodel_s)
emmeans(richmodel_s,pairwise~region) #comparing region richness
#results: sig diff btw Central-North & Central-South
emmeans(richmodel_s,pairwise~Site) #comparing site richness
#results: sig differences - BFB-BAL, BFB-CHA, BFB-SNY, BFB-WLR, DGM-WLR, SSH-BAL, SSH-CHA, SSH-SNY, SSH-WLR, BAL-DAL, CHA-DAL, DAL-SNY, DAL-WLR, WLR-WPR

