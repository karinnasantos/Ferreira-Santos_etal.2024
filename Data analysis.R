setwd (choose.dir())
getwd()

require(vegan)
require(ggplot2)
library(viridis)
require(ggarrange)
library("FactoMineR")
library("factoextra")

#install.packages("ggarrange")
windowsFonts(Times=windowsFont("Times New Roman"))
loadfonts(device="win")

# Functional and taxomic compositios


pca <- rda(ab.hel_samp, scale = T)
summary(pca)$sites
summary(pca)
scores(pca)$species
summary(pca)
plot(pca)

a<-summary(pca)$sites



res.pca <- PCA(ab.hel_samp, graph = FALSE)

res.desc <- dimdesc(res.pca, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc$Dim.1

res.desc$Dim.2


################################################################################
resampling <-read.table ("resampling.txt", h=T)

ab.hel_resamp <- decostand(resampling , "hel")

pca2 <- rda(ab.hel_resamp, scale = T)
summary(pca2)$sites
summary(pca2)
scores(pca2)$species
summary(pca2)
plot(pca2)


res.pca2 <- PCA(ab.hel_resamp, graph = FALSE)

res.desc2 <- dimdesc(res.pca2, axes = c(1,2), proba = 0.05)
# Description of dimension 1
res.desc2$Dim.1

res.desc2$Dim.2

b<-summary(pca2)$sites


#######################
pro_all <- procrustes(pca, pca2)

summary(pro_all)
plot(pro_all)
plot(pro_all, kind=2)
residuals(pro_all)


ctest_all <- data.frame(pca1=pro_all$Yrot[,1],
                        pca2=pro_all$Yrot[,2],xpca1=pro_all$X[,1],
                        xpca2=pro_all$X[,2])



#rowScores<-as.data.frame(pca.summary$sites)
rowScores.bai <- ctest_all[c(1:16, 49, 62:70, 73:80),]
rowScores.inc<-ctest_all[c(17:41, 60, 71, 72, 81:86),]
rowScores.top<-ctest_all[c(42:48, 50:59, 61),]


p <- ggplot()+
 # geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
  #geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_point(data = rowScores.bai, aes(x=pca1, y=pca2),size = 3,
             color = "#111111", fill="#111111", shape=21) +
  geom_point(data = rowScores.bai,aes(x=xpca1, y=xpca2),
             size=3,
             color="#111111", fill="#111111", shape=1) +
  geom_segment(data = rowScores.bai,aes(x=pca1,y=pca2,xend=xpca1,yend=xpca2),arrow=arrow(length=unit(0.2,"cm","inches")))+
  
  geom_point(data = rowScores.inc, aes(x=pca1, y=pca2), size=3,
             color="#a9a9a9", fill="#a9a9a9", shape=17) +
  geom_point(data = rowScores.inc,aes(x=xpca1, y=xpca2), size=3,
             color="#a9a9a9", fill="#a9a9a9", shape=2) +
  geom_segment(data = rowScores.inc,aes(x=pca1,y=pca2,xend=xpca1,yend=xpca2),arrow=arrow(length=unit(0.2,"cm", "inches")))+
  
  geom_point(data = rowScores.top, aes(x=pca1, y=pca2), size=3,
             color="#191970", fill="#191970", shape=22) +
  geom_point(data = rowScores.top, aes(x=xpca1, y=xpca2), size=3,
             color="#191970", fill="#191970", shape=0) +
  geom_segment( data = rowScores.top,aes(x=pca1,y=pca2,xend=xpca1,yend=xpca2),arrow=arrow(length=unit(0.2,"cm","inches")))+
  
  #ggtitle("All habitats")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 15, angle = 0),
        axis.text.y = element_text(color ="black", size = 15, angle = 0))+
  theme(text=element_text(family="Times", face="plain", size=17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p



p <- p +
  xlab("Dim 1") +
  theme(axis.title.x=element_text(angle = 0, size = 17)) + #, face = "bold"
  ylab("Dim 2") +
  theme(axis.title.y=element_text(angle = 90, size = 17))
p


################################################################################
#rowScores<-as.data.frame(pca.summary$sites)
rowScores.bai <- ctest_all[c(1:16, 49, 62:70, 73:80),]

pca1 <- ggplot() +
  #geom_point(data=rowScores, aes(x= CCA1, y= CCA2),
  #          size=4,
  #         fill = "white",
  #        shape=22) +
  #geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
  #geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_point(data = rowScores.bai, aes(x=pca1, y=pca2),size = 3,
             color = "#111111", fill="#111111", shape=21) +
  geom_point(data = rowScores.bai,aes(x=xpca1, y=xpca2),
             size=3,
             color="#111111", fill="#111111", shape=1) +
  geom_segment(data = rowScores.bai,aes(x=pca1,y=pca2,xend=xpca1,yend=xpca2),
               arrow=arrow(length=unit(0.2,"cm","inches")))+
  #ggtitle("Valley")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 15, angle = 0),
        axis.text.y = element_text(color ="black", size = 15, angle = 0))+
  theme(text=element_text(family="Times", face="plain", size=17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pca1

pca1 <- pca1 +
  xlab("Dim 1") +
  theme(axis.title.x=element_text(angle = 0, size = 17)) + #, face = "bold"
  ylab("Dim 2") +
  theme(axis.title.y=element_text(angle = 90, size = 17))
pca1


#####################################################################

rowScores.inc<-ctest_all[c(17:41, 60, 71, 72, 81:86),]

pca2 <- ggplot() +
  #geom_point(data=rowScores, aes(x= CCA1, y= CCA2),
  #          size=4,
  #         fill = "white",
  #        shape=22) +
  #geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
  #geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_point(data = rowScores.inc, aes(x=pca1, y=pca2), size=3,
             color="#a9a9a9", fill="#a9a9a9", shape=17) +
  geom_point(data = rowScores.inc,aes(x=xpca1, y=xpca2), size=3,
             color="#a9a9a9", fill="#a9a9a9", shape=2) +
  geom_segment(data = rowScores.inc,aes(x=pca1,y=pca2,xend=xpca1,yend=xpca2),
               arrow=arrow(length=unit(0.2,"cm")))+
  #ggtitle("slope")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 15, angle = 0),
        axis.text.y = element_text(color ="black", size = 15, angle = 0))+
  theme(text=element_text(family="Times", face="plain", size=17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pca2

pca2 <- pca2 +
  xlab("Dim 1") +
  theme(axis.title.x=element_text(angle = 0, size = 17)) + #, face = "bold"
  ylab("Dim 2") +
  theme(axis.title.y=element_text(angle = 90, size = 17))
pca2

#####################################################################

rowScores.top<-ctest_all[c(42:48, 50:59, 61),]

pca3 <- ggplot() +
  #geom_point(data=rowScores, aes(x= CCA1, y= CCA2),
  #          size=4,
  #         fill = "white",
  #        shape=22) +
 # geom_vline(xintercept = 0, alpha=0.9, linetype = "dashed") +
 # geom_hline(yintercept = 0, alpha=0.9, linetype = "dashed") +
  geom_point(data = rowScores.top, aes(x=pca1, y=pca2), size=3,
             color="#191970", fill="#191970", shape=22) +
  geom_point(data = rowScores.top, aes(x=xpca1, y=xpca2), size=3,
             color="#191970", fill="#191970", shape=0) +
  geom_segment( data = rowScores.top,aes(x=pca1,y=pca2,xend=xpca1,yend=xpca2),
                arrow=arrow(length=unit(0.2,"cm")))+
  # ggtitle("Ridge")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 15, angle = 0),
        axis.text.y = element_text(color ="black", size = 15, angle = 0))+
  theme(text=element_text(family="Times", face="plain", size=17))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


pca3

pca3 <- pca3 +
  xlab("Dim 1") +
  theme(axis.title.x=element_text(angle = 0, size = 17)) + #, face = "bold"
  ylab("Dim 2") +
  theme(axis.title.y=element_text(angle = 90, size = 17))
pca3

############################

library("cowplot")


taxo <- ggarrange(p, pca1, pca2, pca3, 
                  labels = c("A", "B", "C", "D"),
                  ncol = 2, nrow = 2)
taxo

#####################
require(cowplot)
tiff(filename="tax_2_PCAs.tiff",res=600, 
     height=600/72*500, width=450/72*800, compression= "lzw")
plot_grid(
  taxo,
  ncol = 1,
  nrow = 1
)
dev.off()
###########################################################################################
# GLMM
require(lme4)
require(glmm)
library(car)
library(psycho)
library(tidyverse)
library(lmerTest)
require(MASS)
require(postHoc)
library(multcomp)
require(emmeans)


# define variables
PC1 <- as.numeric(tax$PC1)
habitat <- as.factor(tax$Habitat)
date <- as.factor(tax$date)
plot <- as.factor(tax$plot)

hist(PC1)

recog <- tax$PC1 + 1
qqp(recog, "norm")

str(tax)
head(tax)

model1<- lmer(PC1 ~ date + Habitat  +(1 | plot), data = tax,REML=T)
summary(model1)
anova(model1)
#Anova(model1, type=3)
model1.1 <- anova(model1, Type="II")
model1.1

cc <- coef(summary(model1))
str(cc)

fixef(model1)

summary(model1)
Anova(model1, Type="III")

emmeans(model1, pairwise ~ date)

marginal = emmeans(model1,
                   ~ date)

pairs(marginal,
      adjust="tukey")

#install.packages("multcompView")
library(multcompView)
library(multcomp)

cld(marginal,
    alpha=0.05,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey") 



emmeans(model1, pairwise ~ Habitat)

marginal = emmeans(model3,
                   ~ Habitat)

pairs(marginal,
      adjust="tukey")

library(multcomp)

cld(marginal,
    alpha=0.05,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey") 


# Type III anova table with p-values for F-tests based on Satterthwaite's
(aov <- anova(model1))

Anova(model1, Type="II")




## backward elimination of non-significant effects:
step_result <- step(model1)

# Extract the model that step found:
final_model <- get_model(step_result)
summary(final_model)

#####################################
model2<- lmer(PC2 ~ date + Habitat  +(1 | plot), data = tax,REML=FALSE)

summary(model2)
anova(model2)
Anova(model2, type=3)
Anova(model2, Type="II")

vc <- vcov(model2, useScale = TRUE)
b <- fixef(model2)
se <- sqrt(diag(vc))
z <- b / sqrt(diag(vc))
P <- 2 * (1 - pnorm(abs(z)))

cbind(b, se, z, P)


emmeans(model2, pairwise ~ date)

emmeans(model2, pairwise ~ Habitat)

#emmeans(model2, pairwise ~ date*Habitat)

###############################################################################
#Functional

func <-read.table ("func_glmm.txt", h=T)

hist(func$PC1)

func1 <- func$PC1 + 1
qqp(func1, "norm")

# lnorm means lognormal
qqp(func1, "lnorm")
gamma <- fitdistr(func1, "gamma")
qqp(func1, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])

model3 <- lmer(log(PC1) ~ date + Habitat + (1 | plot), data = func)

summary(model3)
#anova(model1)
Anova(model3, type=3)
Anova(model3, Type="II")

anova(model3)

emmeans(model3, pairwise ~ date)

marginal = emmeans(model3,
                   ~ date)

pairs(marginal,
      adjust="tukey")

library(multcomp)

cld(marginal,
    alpha=0.05,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey") 



emmeans(model3, pairwise ~ Habitat)

marginal = emmeans(model3,
                   ~ Habitat)

pairs(marginal,
      adjust="tukey")

library(multcomp)

cld(marginal,
    alpha=0.05,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey") 


## backward elimination of non-significant effects:
step_result2 <- step(model2)

# Extract the model that step found:
final_model2 <- get_model(step_result2)
summary(final_model2)
#####################################
model4<- lmer(PC2 ~ date + Habitat  +(1 | plot), data = func,REML=FALSE)

summary(model4)
anova(model4)
Anova(model4, type=3)
Anova(model4, Type="II")



emmeans(model3, pairwise ~ date)

marginal = emmeans(model3,
                   ~ date)

pairs(marginal,
      adjust="tukey")

library(multcomp)

cld(marginal,
    alpha=0.05,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey") 



emmeans(model4, pairwise ~ Habitat)

marginal = emmeans(model4,
                   ~ Habitat)

pairs(marginal,
      adjust="tukey")

library(multcomp)

cld(marginal,
    alpha=0.05,
    Letters=letters,      ### Use lower-case letters for .group
    adjust="tukey") 


###############################################################################
# define variables
mass_loss <- as.numeric(decomp_glmm$perda_massa)
habitat <- as.factor(decomp_glmm$Habitat)
date <- as.factor(decomp_glmm$date)
plot <- as.factor(decomp_glmm$plot)

hist(mass_loss)

recog <- decomp_glmm$perda_massa + 1
qqp(recog, "norm")

# lnorm means lognormal
qqp(recog, "lnorm")

gamma <- fitdistr(recog, "gamma")
qqp(recog, "gamma", shape = gamma$estimate[[1]], rate = gamma$estimate[[2]])


# log normal
decomp_log <-log1p(decomp_glmm$perda_massa)
str(decomp_log)
decomp_glmm2<-cbind(decomp_glmm, decomp_log)
str(decomp_glmm2)
head(decomp_glmm2)
tail(decomp_glmm2)
decomp_glmm2

model5 <- lmer(decomp_log ~ date + Habitat + (1| plot), data = decomp_glmm2)

summary(model5)
anova(model5)
Anova(model5, type=3)
Anova(model5, Type="II")


vc <- vcov(model5, useScale = TRUE)
b <- fixef(model5)
se <- sqrt(diag(vc))
z <- b / sqrt(diag(vc))
P <- 2 * (1 - pnorm(abs(z)))

cbind(b, se, z, P)

emmeans(model5, pairwise ~ date)

emmeans(model5, pairwise ~ Habitat)
###########################################################################################
# Decompositios experiments
require(vegan)
require(ggplot2)
require(ggpubr)

# Valley
# Gr?fico

jitter_mag <- 0
gg2 <- ggplot(data=baixada, aes(x=year, y=perda_massa, group=parcela, label = parcela )) +
  geom_line(aes(group=parcela), size=.6,
            alpha=1, 
            position=position_jitter(w=0, h=jitter_mag)) +
  # geom_point( )+
  ylab("CWM of mass loss") +
  xlab("") +
  ylim(0, 0.1)+
  theme_pubr() +
  NULL+
  #ggtitle("Valley")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 20, angle = 0),
        axis.text.y = element_text(color ="black", size = 20, angle = 0))+
  theme(text=element_text(face="plain", size=25))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
scale_x_discrete(limits=c("2016(year)", "2020(year)"),
                 labels=c(" ", " "))

gg2

#gg2 <- gg2 + theme(legend.position = c(.15, 0.75))

#gg2

#gg2 <- gg2 + annotate ("text", x=0.7, y=0.1, label="Valley", size=7)
#gg2

gg2 <- gg2 + ylim(-0, 0.1)
gg2

################################################################################
# Gr?fico
jitter_mag <- 0
gg3 <- ggplot(data=slope, aes(x=year, y=perda_massa, group=parcela, label = parcela )) +
  geom_line(aes(group=parcela), size=.6,
            alpha=0.8, 
            position=position_jitter(w=0, h=jitter_mag)) +
  # geom_point( )+
  ylab(" ") +
  xlab("") +
  theme_pubr() +
  NULL+
  #ggtitle("Valley")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 20, angle = 0),
        axis.text.y = element_text(color ="black", size = 20, angle = 0))+
  theme(text=element_text(family="Times", face="plain", size=25))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
#scale_x_discrete(limits=c("2016(year)", "2020(year)"),
#           labels=c(" ", " "))

gg3

gg3 <- gg3 + theme(legend.position="none")

gg3

#gg3 <- gg3 + annotate ("text", x=0.7, y=0.1, label="Slope", size=7)
#gg3

gg3 <- gg3 + ylim(-0, 0.1)
gg3



################################################################################
# Gr?fico
jitter_mag <- 0
gg4 <- ggplot(data=ridge, aes(x=year, y=perda_massa, group=parcela, label = parcela )) +
  geom_line(aes(group=parcela), size=.6,
            alpha=0.8, 
            position=position_jitter(w=0, h=jitter_mag)) +
  # geom_point( )+
  ylab(" ") +
  xlab("") +
  theme_pubr() +
  NULL+
  #ggtitle("Valley")+
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.grid.minor = element_line(colour = NULL),
        panel.grid.major = element_line(colour = NULL),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "black"),
        axis.text.x = element_text(color ="black", size = 20, angle = 0),
        axis.text.y = element_text(color ="black", size = 20, angle = 0))+
  theme(text=element_text(family="Times", face="plain", size=25))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# scale_x_discrete(limits=c("2016(year)", "2020(year)"),
#                  labels=c("During the drought", "Post-drought"))

gg4

gg4 <- gg4 + theme(legend.position="none")

gg4

#gg4 <- gg4 + annotate ("text", x=0.7, y=0.1, label="Ridge", size=7)
#gg4

gg4 <- gg4 + ylim(-0, 0.1)
gg4

#################################################
require(cowplot)

tiff(filename="decomposition_all2.tiff",res=600, 
     height=500/72*500, width=700/72*800, compression= "lzw")
plot_grid(
  gg2,gg3,gg4,
  labels = c('A', 'B', 'C'),
  align="hv",
  ncol = 3,
  nrow = 1
)
dev.off()