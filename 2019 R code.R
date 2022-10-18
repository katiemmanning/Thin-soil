#insect communities

#bring in datasets
bowls <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/Insect%20ID%202019%20-%20Bowl_natural.csv",na.strings = NULL)
summary(bowls)
str(bowls) 
ramps <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/Insect%20ID%202019%20-%20Ramp_natural.csv",na.strings = NULL)
summary(ramps)
str(ramps)
sticky <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/Insect%20ID%202019%20-%20Sticky%20card_natural.csv",na.strings = NULL)
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
allbugs.rowsums <- rowSums(allbugs[,5:50]>0)
allbugs$richness <- allbugs.rowsums

#To obtain abundance counts
allbugs.abun <- rowSums(allbugs[,5:50])
allbugs$abundance <- allbugs.abun

#load vegan
library(vegan)

#calculate Shannon diversity
diversity <-diversity(allbugs[,5:50])
allbugs$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(allbugs[,5:50]))
allbugs$evenness <- evenness

#look at data set
summary(allbugs)
str(allbugs) #trap and region are listed as character 
allbugs$Trap <- as.factor(allbugs$Trap)
allbugs$region <- as.factor(allbugs$region)
str(allbugs) #now trap and region are listed as a factor

#models and checking assumptions
library (emmeans) #for pairwise comparisons
library(lme4)
library(lmerTest) #to obtain p values
library (multcompView) #to view letters

if (!suppressWarnings(require(nortest))) install.packages("nortest")
citation("nortest")
if (!suppressWarnings(require(car))) install.packages("car")
citation("car")
if (!suppressWarnings(require(bbmle))) install.packages("bbmle")
citation("bbmle")
if (!suppressWarnings(require(DHARMa))) install.packages("DHARMa")
citation("DHARMa")
if (!suppressWarnings(require(ggplot2))) install.packages("ggplot2")
citation("ggplot2")
if (!suppressWarnings(require(sjPlot))) install.packages("sjPlot")
citation("sjPlot")
if (!suppressWarnings(require(jtools))) install.packages("jtools")
citation("jtools")
if (!suppressWarnings(require(interactions))) install.packages("interactions")
citation("interactions")

##richness linear mixed effects model
richmodel <- lmer(richness~region + Date + Trap + (1|Site:Replicate), data=allbugs)  #AIC = 1062
summary(richmodel)
AIC(richmodel)
anova(richmodel) #region not sig 

rich.emm<-emmeans(richmodel,pairwise~region) #comparing region richness
rich.emm
#results: South-Central are different, no difference between C-N or N-S
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld 

#check assumptions
dotchart(allbugs$richness, main = "richness", group = allbugs$region) # way to visualize outliers

with(allbugs, ad.test(richness)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = 1.793e-06

with(allbugs, bartlett.test(richness ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.3217

plot(richmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(richmodel))
qqline(resid(richmodel))

plot(simulateResiduals(richmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.58265
#dispersion test: p = 0.584
#outlier test: p = 0.05265
#no significant problems detected 

densityPlot(rstudent(richmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(richmodel)
influenceIndexPlot(richmodel, vars = c("Cook"), id = list(n = 3))

#

##abundance linear mixed effects model
abunmodel <- glmer(abundance~region + Date + Trap + (1|Site:Replicate),data=allbugs, family = negative.binomial(2)) #AIC 2501

summary(abunmodel)
AIC(abunmodel)
anova(abunmodel) #region sig

abun.emm<-emmeans(abunmodel,pairwise~region) #comparing region abundance
abun.emm
#results: same btw central-north, different btw central-south and north-south
abun.cld<-multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS)
abun.cld 

#check assumptions
dotchart(allbugs$abundance, main = "abundance", group = allbugs$region) # way to visualize outliers

with(allbugs, ad.test(abundance)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = < 2.2e-16

with(allbugs, bartlett.test(abundance ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = < 2.2e-16

plot(abunmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(abunmodel))
qqline(resid(abunmodel))

plot(simulateResiduals(abunmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.37758
#dispersion test: p = 0.024 *sig deviation
#outlier test: p = 0.56
#no significant problems detected 

densityPlot(rstudent(abunmodel)) # check density estimate of the distribution of residuals

#Can't use with glmer
# check for outliers influencing the data
outlierTest(abunmodel)
influenceIndexPlot(abunmodel, vars = c("Cook"), id = list(n = 3))

#

##diversity linear mixed effects model
divmodel <- lmer(diversity~region + Date + Trap + (1|Site:Replicate), data=allbugs)  #AIC = 293
summary(divmodel)
AIC(divmodel)
anova(divmodel) #region sig

div.emm<-emmeans(divmodel,pairwise~region) #comparing region diversity
div.emm
#results: same for all
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld 

#check assumptions
dotchart(allbugs$diversity, main = "diversity", group = allbugs$region) # way to visualize outliers

with(allbugs, ad.test(diversity)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = 0.02046

with(allbugs, bartlett.test(diversity ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.9444

plot(divmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(divmodel))
qqline(resid(divmodel))

plot(simulateResiduals(divmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.69357
#dispersion test: p = 0.704
#outlier test: p = 0.72825
#no significant problems detected  

densityPlot(rstudent(divmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(divmodel)
influenceIndexPlot(divmodel, vars = c("Cook"), id = list(n = 3))

#

##evenness linear mixed effects model
evenmodel <- lmer(evenness~region + Date + Trap + (1|Site:Replicate), data=allbugs)  #AIC = -83
summary(evenmodel)
AIC(evenmodel)
anova(evenmodel) #region sig

even.emm<-emmeans(evenmodel,pairwise~region) #comparing region evenness
even.emm
#results: same btw central-north, different btw C-N and N-S
even.cld<-multcomp::cld(even.emm, alpha = 0.05, Letters = LETTERS)
even.cld 

#check assumptions
dotchart(allbugs$evenness, main = "evenness", group = allbugs$region) # way to visualize outliers

with(allbugs, ad.test(evenness)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = 4.478e-05

with(allbugs, bartlett.test(evenness ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.09482

plot(evenmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(evenmodel))
qqline(resid(evenmodel))

plot(simulateResiduals(evenmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.41025
#dispersion test: p = 0.704
#outlier test: p = 1
#no significant problems detected 

densityPlot(rstudent(evenmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(evenmodel)
influenceIndexPlot(evenmodel, vars = c("Cook"), id = list(n = 3))

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
allbugs_boxplot <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                     #labels = c("A", "B", "C", "D"),
                     ncol = 1, nrow = 4,
                     common.legend = TRUE, legend = "bottom")
allbugs_boxplot

pdf("allbugs_boxplot.pdf", height=8, width=8) #height and width in inches
allbugs_boxplot
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

#bring in data pooled by site
bowls_pooled <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/Insect%20ID%202019%20-%20Bowl_natural_pooled.csv",na.strings = NULL)
ramps_pooled <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/Insect%20ID%202019%20-%20Ramp_natural_pooled.csv",na.strings = NULL)
sticky_pooled <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/Insect%20ID%202019%20-%20Sticky%20card_natural_pooled.csv",na.strings = NULL)

#add trap type as a column on each data file
bowls_pooled$Trap="bowl"
ramps_pooled$Trap="ramp"
sticky_pooled$Trap="sticky"

#combine data tables 
library(plyr)
bowlramp_pooled <- rbind.fill (bowls_pooled, ramps_pooled)
allbugs_pooled <-rbind.fill (bowlramp_pooled, sticky_pooled)

#add column for region (south central north)
allbugs_pooled$region<-ifelse(allbugs_pooled$Site=="WLR", "South",
                       ifelse(allbugs_pooled$Site=="WPR", "South",
                              ifelse(allbugs_pooled$Site=="SNY", "South",
                                     ifelse(allbugs_pooled$Site=="DAL", "North", 
                                            ifelse(allbugs_pooled$Site=="BAL", "North",
                                                   ifelse(allbugs_pooled$Site == "CHA", "North", "Central")
                                            )))))
str(allbugs_pooled)

#Create matrix of environmental variables    
env.matrix<-allbugs_pooled[c(1:2,49:50)]

#create matrix of community variables
com.matrix<-allbugs_pooled[c(3:48)]

#change to presence/absence
com.matrix[com.matrix > 0] <- 1
str(com.matrix)
rowSums(com.matrix)

#ordination by NMDS
NMDS<-metaMDS(com.matrix, distance="bray", k=2, autotransform=TRUE, trymax=300)
NMDS
###stress = 0.22
stressplot(NMDS)

#plot NMDS for region
#might need to change colors
#8 x 13
plot(NMDS, disp='sites', type="n")
#title(main="Arthropod community composition by region", cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South")
ordiellipse(NMDS, env.matrix$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS, env.matrix$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
#add data points
points(NMDS, display="sites", select=which(env.matrix$region=="North"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(env.matrix$region=="Central"), pch=17, col="#009E73")
points(NMDS, display="sites", select=which(env.matrix$region=="South"), pch=15, col="#CC79A7")
#add legend
legend(0.888,0.83, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=1.5, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (regions)
fit<-adonis(com.matrix ~ region, data = env.matrix, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$region))
#P-value = 0.001 -- cannot assume homogeneity of multivariate dispersion

install.packages("devtools")
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
citation("pairwiseAdonis")

pairwise.adonis(com.matrix, env.matrix$region) #south-central sig

###

#plot NMDS for Northern sites 
#8 x 13
plot(NMDS, disp='sites', type="n")
title(main="North", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
#add data points
points(NMDS, display="sites", select=which(env.matrix$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS, display="sites", select=which(env.matrix$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS, display="sites", select=which(env.matrix$Site=="DAL"), pch=15, col="#33A02C")
#add legend
legend(0.7,1.12, title=NULL, pch=c(19,17,15), col=c("#A6CEE3","#B2DF8A","#33A02C"), cex=1.55, legend=c("Beaton alvar", "Cape Hurd alvar", "Davis alvar"))

#plot NMDS for Central sites 
#8 x 13
plot(NMDS, disp='sites', type="n")
title(main="Central", adj = 0.01, line = -2, cex.main=2.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
#add data points
points(NMDS, display="sites", select=which(env.matrix$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS, display="sites", select=which(env.matrix$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS, display="sites", select=which(env.matrix$Site=="BFB"), pch=15, col="#1F78B4")
#add legend
legend(0.40,1.12, title=NULL, pch=c(19,17,15), col=c("#FDBF6F","#FB9A99","#1F78B4"), cex=1.55, legend=c("Slate shale hill", "Dusty Goldenrod meadow", "Bedford barren"))

#plot NMDS for Southern sites 
#8 x 13
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
legend(0.85,1.12, title=NULL, pch=c(19,17,15), col=c("#FF7F00","#CAB2D6","#E31A1C"), cex=1.55, legend=c("W ladder", "W picnic rock", "Synder hollow"))

#bootstrapping and testing for differences between the groups (sites)
fit<-adonis(com.matrix ~ Site, data = env.matrix, permutations = 999, method="bray")
fit
#P= 0.01

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Site))
#P-value = 0.1672 -- assumes homogeneity of multivariate dispersion

pairwise.adonis(com.matrix, env.matrix$Site) #none sig


#multi-panel NMDS with just sites
library(ggvegan)

pdf("site NMDS.pdf", height=6.5, width=15)
par(mfrow=c(1,3), mar=c(4.1, 4.8, 1.5, 8.1),xpd=TRUE) 

plot(NMDS, disp='sites', type="n")
title(main="South", adj = 0.01, line = -2, cex.main=2.5)
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#CAB2D6",kind="sd", conf=0.95, label=FALSE, show.groups = "WPR")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#E31A1C",kind="sd", conf=0.95, label=FALSE, show.groups = "SNY")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FF7F00",kind="sd", conf=0.95, label=FALSE, show.groups = "WLR")
points(NMDS, display="sites", select=which(env.matrix$Site=="WPR"), pch=17, col="#CAB2D6")
points(NMDS, display="sites", select=which(env.matrix$Site=="SNY"), pch=15, col="#E31A1C")
points(NMDS, display="sites", select=which(env.matrix$Site=="WLR"),pch=19, col="#FF7F00")
legend(-.45,-1, title=NULL, pch=c(19,17,15), col=c("#FF7F00","#CAB2D6","#E31A1C"), cex=1.55, legend=c("W ladder", "W picnic rock", "Synder hollow"))

plot(NMDS, disp='sites', type="n")
title(main="Central", adj = 0.01, line = -2, cex.main=2.5)
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FB9A99",kind="sd", conf=0.95, label=FALSE, show.groups = "DGM")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#1F78B4",kind="sd", conf=0.95, label=FALSE, show.groups = "BFB")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#FDBF6F",kind="sd", conf=0.95, label=FALSE, show.groups = "SSH")
points(NMDS, display="sites", select=which(env.matrix$Site=="SSH"),pch=19, col="#FDBF6F")
points(NMDS, display="sites", select=which(env.matrix$Site=="DGM"), pch=17, col="#FB9A99")
points(NMDS, display="sites", select=which(env.matrix$Site=="BFB"), pch=15, col="#1F78B4")
legend(-.65,-1, title=NULL, pch=c(19,17,15), col=c("#FDBF6F","#FB9A99","#1F78B4"), cex=1.55, legend=c("Slate shale hill", "Dusty Goldenrod meadow", "Bedford barren"))

plot(NMDS, disp='sites', type="n")
title(main="North", adj = 0.01, line = -2, cex.main=2.5)
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#B2DF8A",kind="sd", conf=0.95, label=FALSE, show.groups = "CHA") 
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#33A02C",kind="sd", conf=0.95, label=FALSE, show.groups = "DAL")
ordiellipse(NMDS, env.matrix$Site, draw="polygon", col="#A6CEE3",kind="sd", conf=0.95, label=FALSE, show.groups = "BAL")
points(NMDS, display="sites", select=which(env.matrix$Site=="BAL"),pch=19, col="#A6CEE3")
points(NMDS, display="sites", select=which(env.matrix$Site=="CHA"), pch=17, col="#B2DF8A")
points(NMDS, display="sites", select=which(env.matrix$Site=="DAL"), pch=15, col="#33A02C")
legend(-.5,-1, title=NULL, pch=c(19,17,15), col=c("#A6CEE3","#B2DF8A","#33A02C"), cex=1.55, legend=c("Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
dev.off()

######

#Bee analysis (species/genus)
#bring in data 
bee_bowls <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/2019%20bees%20-%20Bowl_species.csv",na.strings = NULL)
summary(bee_bowls)
str(bee_bowls) 
bee_ramps <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/2019%20bees%20-%20Ramp_species.csv",na.strings = NULL)
summary(bee_ramps)
str(bee_ramps)

#add trap type as a column on each data file
bee_bowls$Trap="bowl"
bee_ramps$Trap="ramp"

#combine data tables 
library (plyr)
bees <- rbind.fill (bee_bowls, bee_ramps)

#add column for region (south central north)
bees$region<-ifelse(bees$Site=="WLR", "South",
                       ifelse(bees$Site=="WPR", "South",
                              ifelse(bees$Site=="SNY", "South",
                                     ifelse(bees$Site=="DAL", "North", 
                                            ifelse(bees$Site=="BAL", "North",
                                                   ifelse(bees$Site == "CHA", "North", "Central")
                                            )))))
str(bees)

#To obtain richness counts
bees.rowsums <- rowSums(bees[,4:29]>0)
bees$richness <- bees.rowsums

#To obtain abundance counts
bees.abun <- rowSums(bees[,4:29])
bees$abundance <- bees.abun

#load vegan
library(vegan)

#calculate Shannon diversity
diversity <-diversity(bees[,4:29])
bees$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(bees[,4:29]))
bees$evenness <- evenness

#look at data set
summary(bees)
str(bees)

##richness linear mixed effects model
library (emmeans) #for pairwise comparisons
library(lme4)
library(lmerTest) #to obtain p values

richmodel <- lm(richness~region+Date+Trap+Site, data=bees)  #AIC = 171
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
#results: same for all
rich.cld.s<-multcomp::cld(rich.emm.s, alpha = 0.05, Letters = LETTERS)
rich.cld.s

#check assumptions
dotchart(bees$richness, main = "richness", group = bees$region) # way to visualize outliers

with(bees, ad.test(richness)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = 2.2e-16

with(bees, bartlett.test(richness ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.0001

plot(richmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(richmodel))
qqline(resid(richmodel))

plot(simulateResiduals(richmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.46063
#dispersion test: p = 0.168
#outlier test: p = 1
#no significant problems detected 

densityPlot(rstudent(richmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(richmodel)
influenceIndexPlot(richmodel, vars = c("Cook"), id = list(n = 3))

#

##abundance linear mixed effects model
abunmodel <- lm(abundance~region+Date+Trap+Site, data=bees)  #AIC = 252
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
#results: same for all
abun.cld.s<-multcomp::cld(abun.emm.s, alpha = 0.05, Letters = LETTERS)
abun.cld.s

#check assumptions
dotchart(bees$abundance, main = "abundance", group = bees$region) # way to visualize outliers

with(bees, ad.test(abundance)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = 5.013e-13

with(bees, bartlett.test(abundance ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.00099 

plot(abunmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(abunmodel))
qqline(resid(abunmodel))

plot(simulateResiduals(abunmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.06964
#dispersion test: p = 0.168
#outlier test: p = 0.3611
#no significant problems detected 

densityPlot(rstudent(abunmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(abunmodel)
influenceIndexPlot(abunmodel, vars = c("Cook"), id = list(n = 3))

#

##diversity linear mixed effects model
divmodel <- lm(diversity~region+Date+Trap+Site,data=bees)  #AIC = 67
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
#results: same for all
div.cld.s<-multcomp::cld(div.emm.s, alpha = 0.05, Letters = LETTERS)
div.cld.s

#check assumptions
dotchart(bees$diversity, main = "diversity", group = bees$region) # way to visualize outliers

with(bees, ad.test(diversity)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = < 2.2e-16

with(bees, bartlett.test(diversity ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.036

plot(divmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(divmodel))
qqline(resid(divmodel))

plot(simulateResiduals(divmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.56221
#dispersion test: p = 0.168
#outlier test: p = 1
#no significant problems detected  

densityPlot(rstudent(divmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(divmodel)
influenceIndexPlot(divmodel, vars = c("Cook"), id = list(n = 3))

#

##evenness linear mixed effects model
evenmodel <- lm(evenness~region+Date+Trap+Site,data=bees)  #AIC = -34
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
#results: same for all
even.cld.s<-multcomp::cld(even.emm.s, alpha = 0.05, Letters = LETTERS)
even.cld.s

#check assumptions
dotchart(bees$evenness, main = "evenness", group = bees$region) # way to visualize outliers

with(bees, ad.test(evenness)) #Anderson-darling test for normality (good for small sample sizes), low p-value means assumption is violated
#p-value = 0.00028

with(bees, bartlett.test(evenness ~ region)) #Bartlett test for homogeneity of variance, low p-value means assumption is violated
#p-value = 0.1592

plot(evenmodel) # check distribution of residuals

# check normality with these figures, are there outliers at either end
qqnorm(resid(evenmodel))
qqline(resid(evenmodel))

plot(simulateResiduals(evenmodel)) # another way to check for normailty and homogeneity of variance
#KS test: p = 0.44695
#dispersion test: p = 0.088
#outlier test: p = 1
#no significant problems detected 

densityPlot(rstudent(evenmodel)) # check density estimate of the distribution of residuals

# check for outliers influencing the data
outlierTest(evenmodel)
influenceIndexPlot(evenmodel, vars = c("Cook"), id = list(n = 3))

##
library (ggplot2)

#site richness by region
richness.plot<-ggplot(bees, aes(x = factor(region,level = c("South","Central","North")), y = richness, fill=Site))+
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
abundance.plot<-ggplot(bees, aes(x = factor(region, level = c("South","Central","North")), y = abundance, fill=Site))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="bottom")+
  labs(title="", x="", y="Abundance")+
  #theme (plot.title = element_text(hjust=0.5))+
  #scale_y_continuous(trans="log10")+
  scale_fill_brewer(palette="Paired",name="Sites:",
                    breaks=c("SNY", "WLR", "WPR", "BFB", "DGM", "SSH", "BAL", "CHA", "DAL"),
                    labels=c("Snyder hollow","W ladder", "W picnic rock", "Bedford barren","Dusty goldenrod meadow", "Slate shale hill", "Beaton alvar", "Cape Hurd alvar", "Davis alvar"))
abundance.plot

#site diversity by region
diversity.plot<-ggplot(bees, aes(x = factor(region,level = c("South","Central","North")), y = diversity, fill=Site))+
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
evenness.plot<-ggplot(bees, aes(x = factor(region,level = c("South","Central","North")), y = evenness, fill=Site))+
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
bees_boxplot <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                             #labels = c("A", "B", "C", "D"),
                             ncol = 1, nrow = 4,
                             common.legend = TRUE, legend = "bottom")
bees_boxplot

pdf("bees_boxplot.pdf", height=8, width=8) #height and width in inches
bees_boxplot
dev.off()
###

#NMDS of bee community 
library (vegan)

#Create matrix of environmental variables
env.matrix_bees<-bees[c(1:3,30:35)]
#create matrix of community variables
com.matrix_bees<-bees[c(4:29)]

#ordination by NMDS
NMDS_bees<-metaMDS(com.matrix_bees, distance="bray", k=2, autotransform=TRUE, trymax=300)
NMDS_bees
###stress =  8.59e-05 -- insufficient data
stressplot(NMDS_bees)

#plot bee NMDS for region
#8 x 10
plot(NMDS_bees, disp='sites', type="n")
#title(main="Bee community composition by region", cex.main=1.5)
#add ellipsoids with ordiellipse
ordiellipse(NMDS_bees, env.matrix_bees$region, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "South")
ordiellipse(NMDS_bees, env.matrix_bees$region, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "North")
ordiellipse(NMDS_bees, env.matrix_bees$region, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "Central") 
#add data points
points(NMDS_bees, display="sites", select=which(env.matrix_bees$region=="North"),pch=19, col="#E69F00")
points(NMDS_bees, display="sites", select=which(env.matrix_bees$region=="Central"), pch=17, col="#009E73")
points(NMDS_bees, display="sites", select=which(env.matrix_bees$region=="South"), pch=15, col="#CC79A7")
#add legend
#legend(8.9,1.12, title=NULL, pch=c(19,17,15), col=c("#E69F00","#009E73","#CC79A7"), cex=1.5, legend=c("North", "Central", "South"))

#bootstrapping and testing for differences between the groups (regions)
fit<-adonis(com.matrix_bees ~ region, data = env.matrix_bees, permutations = 999, method="bray")
fit
#P=0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_bees)
anova(betadisper(distances_data, env.matrix_bees$region))
#P-value = 0.3521 -- assumes homogeneity of multivariate dispersion

pairwise.adonis(com.matrix_bees, env.matrix_bees$region) #none sig
pairwise.adonis(com.matrix_bees, env.matrix_bees$Site) #none sig

###

#species accumulation
library (BiodiversityR)
library(ggplot2)

#separate bees by region
#bring in data

north <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/2019%20bees_north.csv",na.strings = NULL)
central <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/2019%20bees_central.csv",na.strings = NULL)
south <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/2019%20bees_south.csv",na.strings = NULL)

#individual curves for each region
north.com.matrix<-north[c(4:29)]
north_curve<-accumresult(north.com.matrix, method = "exact", permutations = 1000)

central.com.matrix<-central[c(4:29)]
central_curve<-accumresult(central.com.matrix, method = "exact", permutations = 1000)

south.com.matrix<-south[c(4:29)]
south_curve<-accumresult(south.com.matrix, method = "exact", permutations = 1000)


#first-order jackknife estimates are based on the number of singletons
#second-order jackknife estimates are based on the number of singletons and doubletons

#calculates species richness for each sample
specnumber(com.matrix_bees) #ranges from 1 to 6

#calculates species richness by treatment (region)
specnumber(com.matrix_bees, groups = bees$region) #north=5; central=19; south=16

#total richness and jackknife
rich <- diversityresult(com.matrix_bees, y=NULL, index = "richness")
rich # 26
j1 <- diversityresult(com.matrix_bees, y=NULL, index = "jack1")
j1 # 37.785714
#69%
j2 <- diversityresult(com.matrix_bees, y=NULL, index = "jack2")
j2 # 45.57013
#57%

#north jackknife; richness = 5
j1.n <- diversityresult(north.com.matrix, y=NULL, index = "jack1")
j1.n # 6.8
#74%
j2.n <- diversityresult(north.com.matrix, y=NULL, index = "jack2")
j2.n # 7.6888889
#65%

#central jackknife; richness = 19
j1.c <- diversityresult(central.com.matrix, y=NULL, index = "jack1")
j1.c # 26.733333
#71%
j2.c <- diversityresult(central.com.matrix, y=NULL, index = "jack2")
j2.c # 32.397701
#59%

#south jackknife; richness = 16
j1.s <- diversityresult(south.com.matrix, y=NULL, index = "jack1")
j1.s # 25.375
#63%
j2.s <- diversityresult(south.com.matrix, y=NULL, index = "jack2")
j2.s # 32.491667
#49%

#BiodiversityR::accumcomp
Accum.1_functional <- accumcomp(com.matrix_bees, y=env.matrix_bees, factor='region', 
                                method='random', conditioned=FALSE, plotit=FALSE)
Accum.1_functional

#BiodiversityR::accumcomp.long
accum.long1_functional <- accumcomp.long(Accum.1_functional, ci=NA, label.freq=5)
head(accum.long1_functional)

#plot
#empty canvas
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

accum <- ggplot(data=accum.long1_functional, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_color_manual(values=c("#009E73","#E69F00","#CC79A7"))+
  scale_shape_manual(values=c(19,17,15,25))+
  geom_line(aes(colour=Grouping), size=0.1) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE, linetype = 0) + 
  geom_point(data=subset(accum.long1_functional, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) +
  BioR.theme +
  labs(x = "Number of Samples", y = "Bee Richness", colour = "region", shape = "region")
accum

pdf("accum curve.pdf", height=6, width=8) #height and width in inches
accum
dev.off()


######

#Plant analysis

#bring in plant data (presence/absence)
plants <- read.csv("https://raw.githubusercontent.com/katiemmanning/Thin-soil/main/Data/2019%20plants_new.csv",na.strings = NULL)

#To obtain richness counts -- number of plants at each site
plants.rowsums <- rowSums(plants[,4:42]>0)
plants$richness <- plants.rowsums

#look at data set
summary(plants)
str(plants)

#indicator species analysis
library (indicspecies)

abund = plants[,4:ncol(plants)]
region = plants$Region
site = plants$Site

inv.region = multipatt(abund, region, func = "r.g", control = how(nperm=9999))
summary (inv.region)
#Aster for central
#Gaultheria procumbens and Pinus virginiana for south
#none for north

inv.site = multipatt(abund, site, func = "r.g", control = how(nperm=9999))
summary (inv.site)
#none