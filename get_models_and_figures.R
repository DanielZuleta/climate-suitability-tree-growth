#### This scripts generate manin results for the paper
# After estimating growth.l5 and suitability (i.e., code: ), this code tests their relationships 
# Daniel Zuleta, Manuel Bernal, Ken Feeley

rm(list = ls()) # removes all objects in ram

######################################
################ PATHS ###############
######################################

if(grepl("/Users/danielzuleta", getwd())) { 
  path.to.data <- "//Users/danielzuleta/Google Drive/Coauthoring papers/1. ITRDB Project-Metanalisis/growth-suitability-DZ/data/"
  path.to.results <- "/Users/danielzuleta/Google Drive/Coauthoring papers/1. ITRDB Project-Metanalisis/growth-suitability-DZ/results/"
}

########################################################
############ LOAD PACKAGES AND FUNCTIONALITY ###########
########################################################

library(lme4)
library(lmerTest)
library(cAIC4)
library(nlme)
library(ggplot2)
library(ggeffects)
library(MuMIn)
library(lattice)
library(sjPlot)

####Control to converge in lmer
controlConver_daniel=lmerControl(optimizer = "Nelder_Mead",
                                 check.conv.grad = "ignore",
                                 check.conv.singular = "ignore",
                                 check.conv.hess = "ignore",
                                 check.scaleX = "ignore",
                                 optCtrl = list(maxfun = 1000000),
                                 check.nobs.vs.nlev="ignore"
)

tr <- function(color, alpha=0.5)
{
  # add transparency to a given color
  rgb <- as.vector(col2rgb(color))
  return(rgb(rgb[1], rgb[2], rgb[3], alpha=alpha*255, maxColorValue=255))
}

######################################
############### LOAD DATA ############
######################################

crecimiento.aic <- read.csv(paste0(path.to.data, "crecimiento-aic.csv"))
crecimiento.aic$sp.plot <- paste(crecimiento.aic$species, crecimiento.aic$plot, sep = "_")

#####################################
####### PHYLUM: GYMNO & ANGIO #######
#####################################
taxo <- read.csv(paste0(path.to.data, "sp_list.csv")) 
crec.old <- read.csv(paste0(path.to.data, "crecimiento_varianza.csv")) # contains the taxo
# add taxo
taxo$genus <- crec.old$Genus[match(taxo$Taxon.name, crec.old$Species)]
taxo$family <- crec.old$Family[match(taxo$Taxon.name, crec.old$Species)]
taxo$phylum <- crec.old$Phylum[match(taxo$Taxon.name, crec.old$Species)]
names(taxo)[which(names(taxo)=='Taxon.name')] <- 'species'

# adding the phylum, genus, family to the crec table
crecimiento.aic$phylum <- crec.old$Phylum[match(crecimiento.aic$species, crec.old$Species)]
crecimiento.aic$family <- crec.old$Family[match(crecimiento.aic$species, crec.old$Species)]
crecimiento.aic$genus <- crec.old$Genus[match(crecimiento.aic$species, crec.old$Species)]

######################################
######## MODELS ########
######################################

# 1. G~S across space
# -----
mod.gs.across.space <- lmer(growth.l5 ~ suitability.m1915 +
                       (1 + suitability.m1915 | species / year), 
                     data = crecimiento.aic, 
                     control=controlConver_daniel,
                     REML = F)
summary(mod.gs.across.space)
r.squaredGLMM(mod.gs.across.space) ####
dotplot(ranef(mod.gs.across.space, condVar = T), strip = T, scales = list(relation='free'))$species
ranef(mod.gs.across.space, condVar = T)$species
tab_model(mod.gs.across.space, digits=4, digits.p=4, digits.re=4)
anova(mod.gs.across.space)

# 2. G~S across time 
# -----
mod.gs.across.time <- lmer(growth.l5 ~ suitability.m1915 +
                              (1 + suitability.m1915 | species / plot), 
                            data = crecimiento.aic, 
                            control=controlConver_daniel,
                            REML = F)

# 3. G~S across time 
# -----
mod.gs.across.space.and.time <- lmer(growth.l5 ~ suitability.m1915 +
                             (1 + suitability.m1915 | species), 
                           data = crecimiento.aic, 
                           control=controlConver_daniel,
                           REML = F)

# 4. G~T (1 + T|species) 
# --------
mod.gt <- lmer(growth.l5 ~ year + # fixed
                   (1 + year | species), 
                 data = crecimiento.aic, 
                 control = controlConver_daniel,
                 REML = F)

# 5. G~S*T + (1 + S | species/plot)  
mod.gst <- lmer(growth.l5 ~ suitability.m1915*year + # fixed
                       (1 + suitability.m1915 | species/plot),
                     data = crecimiento.aic,
                     control = controlConver_daniel,
                     REML = F)

#### An additional model of interest: Suitabiliy ~ time
mod.st <- lmer(suitability.m1915 ~ year + # fixed
                  (1 + year | species),
                data = crecimiento.aic,
                control = controlConver_daniel,
                REML = F)

#####  Models phylum: based on variance decomposition

##### mod.gstp
mod.gstp <- lmer(growth.l5 ~ suitability.m1915*year*phylum + # fixed
                    (1 + suitability.m1915 | species/plot),
                   data = crecimiento.aic,
                   control = controlConver_daniel,
                   REML = F)

######################################
######## FIGURES ########
######################################

##### Fig 2 ###########
pdf(paste0(path.to.results, "Fig2.pdf"), width = 10, height = 4)
par(mfrow = c(1,2), oma = c(0,0,0,0))

# Panel a: M5 over time
par(mar = c(4.2,5.1,2.5, 1), xpd = T)
pred <- as.data.frame(ggpredict(mod.gst, terms = c('year', 'suitability.m1915'), type = "fe", fullrange = TRUE))
col.gs <- c('darkgoldenrod3','brown4','darkolivegreen4')
plot(crecimiento.aic$year, crecimiento.aic$growth.l5, pch = 16, col = tr("grey30"), cex = .5,
     ylim = c(-0.14,0.2), xlim = c(1915,1995),
     type = "n",
     las = 1, ylab = "", xlab = "", xaxt='n', yaxt='n', ann=FALSE)
axis(1, at = seq(1915, 1995, 20),  line = 0, cex = 1.1)
axis(2, at = seq(-0.10, 0.2, 0.05),  line = 0, las  = 2, cex = 1.1)
mtext("Standardized growth rate", side = 2, line = 3.3, cex = 1.3)
mtext("Year", side = 1, line = 2.7, cex = 1.3)
for(g in 1:length(unique(pred$group))){
  g1 <- unique(pred$group)[g]
  p1 <- pred[pred$group == g1,]
  polygon(c(p1$x, rev(p1$x)),
          c(p1$conf.low, rev(p1$conf.high)),
          col = tr(col.gs[g],0.4), border = NA)
}
for(g in 1:length(unique(pred$group))){
  g1 <- unique(pred$group)[g]
  p1 <- pred[pred$group == g1,]
  lines(p1$x, p1$predicted, col = col.gs[g], lwd = 2)
}
legend('topright', c("0.35","0.49", "0.63"), lty = 1, col = col.gs, horiz = F, bty = "n", lwd = 2.5)
text(1905, 0.24, "(a)", font = 2, pos = 4, xpd = T, cex = 1.5)

# Panel b: M6 over time
par(mar = c(4.2,3.1,2.5, 3), xpd = T)
pred <- as.data.frame(ggpredict(mod.gstp, terms = c('year', 'phylum'), type = "fe", fullrange = TRUE))
lines.phy <- c(1, 2)
plot(crecimiento.aic$year, crecimiento.aic$growth.l5, pch = 16, col = tr("grey30"), cex = .5,
     ylim = c(-0.14,0.2), xlim = c(1915,1995),
     type = "n",
     las = 1, ylab = "", xlab = "", xaxt='n', yaxt='n', ann=FALSE)
axis(1, at = seq(1915, 1995, 20),  line = 0, cex = 1.1)
axis(2, at = seq(-0.10, 0.2, 0.05),  line = 0, las  = 2, cex = 1.1)
# mtext("Standardized growth rate", side = 2, line = 3.5, cex = 1.1)
mtext("Year", side = 1, line = 2.7, cex = 1.3)
for(g in 1:length(unique(pred$group))){
  g1 <- unique(pred$group)[g]
  p1 <- pred[pred$group == g1,]
  polygon(c(p1$x, rev(p1$x)),
          c(p1$conf.low, rev(p1$conf.high)),
          col = tr("grey80",0.4), border = NA)
}
for(g in 1:length(unique(pred$group))){
  g1 <- unique(pred$group)[g]
  p1 <- pred[pred$group == g1,]
  lines(p1$x, p1$predicted, lty = lines.phy[g], lwd = 2)
}
legend('topright', c("Angiosperms", "Gymnosperms"), lty = lines.phy, horiz = F, bty = "n", lwd = 2.5)
text(1905, 0.24, "(b)", font = 2, pos = 4, xpd = T, cex = 1.5)
dev.off()


#### Fig3 ###########
plot(ggpredict(mod.gstp, terms = c('suitability.m1915', 'year', 'phylum'),type = "fe", fullrange = TRUE),
     show.legend  = T, show.title = T ) +
  theme_classic() +
  ylab('Basal area growth.l5 rate')  +  labs(title = "temp")

pdf(paste0(path.to.results, "Fig3.pdf"), width = 10, height = 4)
par(mfrow = c(1,2), oma = c(0,0,0,0))
pred <- as.data.frame(ggpredict(mod.gstp, terms = c('suitability.m1915', 'year', 'phylum'), type = "fe", fullrange = TRUE))
pred.angio <- pred[pred$facet == 'Angiosperm',]
pred.gymno <- pred[pred$facet == 'Gymnosperm',]
col.gs <- c('brown4','royalblue3','darkolivegreen4')
# plot a. Angio
par(mar = c(4.2,5.1,2.5, 1), xpd = T)
plot(crecimiento.aic$suitability.m1915, crecimiento.aic$growth.l5, pch = 16, col = tr("grey30"), cex = .5,
     ylim = c(-0.14,0.2), xlim = c(0,1),
     type = "n",
     las = 1, ylab = "", xlab = "", yaxt='n', ann=FALSE)
axis(2, at = seq(-0.10, 0.2, 0.05),  line = 0, las  = 2, cex = 1.1)
mtext("Standardized growth rate", side = 2, line = 3.3, cex = 1.3)
mtext("Climate suitability", side = 1, line = 2.7, cex = 1.3)
for(g in 1:length(unique(pred.angio$group))){
  g1 <- unique(pred.angio$group)[g]
  p1 <- pred.angio[pred.angio$group == g1,]
  polygon(c(p1$x, rev(p1$x)),
          c(p1$conf.low, rev(p1$conf.high)),
          col = tr(col.gs[g],0.4), border = NA)
}
for(g in 1:length(unique(pred.angio$group))){
  g1 <- unique(pred.angio$group)[g]
  p1 <- pred.angio[pred.angio$group == g1,]
  lines(p1$x, p1$predicted, col = col.gs[g], lwd = 2)
}
text(-0.13, 0.24, "(a) Angiosperms", font = 2, pos = 4, xpd = T, cex = 1.5)

# plot b. Gymno
par(mar = c(4.2,3.1,2.5, 3), xpd = T)
plot(crecimiento.aic$suitability.m1915, crecimiento.aic$growth.l5, pch = 16, col = tr("grey30"), cex = .5,
     ylim = c(-0.14,0.2), xlim = c(0,1),
     type = "n",
     las = 1, ylab = "", xlab = "", yaxt='n', ann=FALSE)
axis(2, at = seq(-0.10, 0.2, 0.05),  line = 0, las  = 2, cex = 1.1)
# mtext("Standardized growth rate", side = 2, line = 3.5, cex = 1.1)
mtext("Climate suitability", side = 1, line = 2.7, cex = 1.3)
for(g in 1:length(unique(pred.gymno$group))){
  g1 <- unique(pred.gymno$group)[g]
  p1 <- pred.gymno[pred.gymno$group == g1,]
  polygon(c(p1$x, rev(p1$x)),
          c(p1$conf.low, rev(p1$conf.high)),
          col = tr(col.gs[g],0.4), border = NA)
}
for(g in 1:length(unique(pred.gymno$group))){
  g1 <- unique(pred.gymno$group)[g]
  p1 <- pred.gymno[pred.gymno$group == g1,]
  lines(p1$x, p1$predicted, col = col.gs[g], lwd = 2)
}
text(-0.13, 0.24, "(b) Gymnosperms", font = 2, pos = 4, xpd = T, cex = 1.5)

legend('topright', c("1932","1956", "1980"), lty = 1, col = col.gs, horiz = F, bty = "n", lwd = 2.5)
dev.off()


#### Fig4 ###########

# # we want this fig:
# dotplot(ranef(mod.gstp, condVar = T), strip = T, scales = list(relation='free'))$species

model = mod.gstp 
re <- ranef(mod.gstp, condVar = T)$species
re$phylum <- taxo$phylum[match(rownames(re), taxo$species)]
ranvar <- attr(re, "postVar")
# If you’re willing to assume independence of the conditional variance and the fixed-effect
# sampling variance, then the variance of the intercepts for each group would be the sum
# of the fixed-effect intercept variance and the conditional variance of the intercept for each group:
vcov(model)[1,1]+ranvar[1,1,1]  # for the first group

vcov(model)[2,2]+ranvar[2,2,1]  # this would be for the slopes in the first group
# ranvar is the variance-covariance matrix of the random effects
# vcov(model) is the variance-covariance matrix of the fixed effects
int.var <- c()
slope.var <- c()
for (i in 1:nrow(re)) {
  int.var <- c(int.var, ranvar[1,1,i]) #  [1,1] is the variance of the the intercept
  slope.var <- c(slope.var, ranvar[2,2,i]) # [2,2] is the variance of the the slope
}

re$int.se <- 1.96*sqrt(int.var) 
re$slope.se <- 1.96*sqrt(slope.var) 

### get the CIs
re$int.se.final <- re$int.se
re$slope.se.final <- re$slope.se
names(re)[which(names(re) == "(Intercept)")] <- "intercept"
re$int.lo <- re$intercept - re$int.se.final
re$int.up <- re$intercept + re$int.se.final 
re$slope.lo <- re$suitability.m1915 - re$slope.se.final 
re$slope.up <- re$suitability.m1915 + re$slope.se.final

# intercept
pdf(paste0(path.to.results, "fig4.pdf"), width = 10, height = 8)
par(oma = rep(0,4), mfrow = c(1,2))
re <- re[order(re$intercept),]
par(mar = c(4, 10, 1, 0.1))
plot(0:nrow(re), 0:nrow(re), xlim = c(min(re$int.lo), max(re$int.up)), pch = NA, yaxt = "n", xaxt = "n",
     ylab = "", xlab = "", ylim = c(1, nrow(re)))
abline(h = 1:nrow(re), lty = 1, col = "grey90")
abline(v = 0, lty = 1, col = "grey80")
points(re$intercept, 1:nrow(re), pch = ifelse(re$phylum == "Angiosperm", 15, 17), cex = 1.5)
segments(x0 = re$int.lo, x1 = re$int.up, y0 = 1:nrow(re), lwd = 1.5)
axis(1, at = seq(-0.4, 0.4, 0.2),  line = 0, las  = 1, cex = 1.1)
axis(2, at = 1:nrow(re),  line = 0, las  = 1, cex = 1.1, labels = rownames(re), font = 3)
mtext(side = 1, "Intercept", line = 2.6, cex = 1.3)
# slope
par(mar = c(4, 0.1, 1, 10))
plot(0:nrow(re), 0:nrow(re), xlim = c(min(re$slope.lo), max(re$slope.up)), pch = NA, yaxt = "n", xaxt = "n",
     ylab = "", xlab = "", ylim = c(1, nrow(re)))
abline(h = 1:nrow(re), lty = 1, col = "grey90")
abline(v = 0, lty = 1, col = "grey80")
points(re$suitability.m1915, 1:nrow(re), pch = ifelse(re$phylum == "Angiosperm", 15, 17), cex = 1.5)
segments(x0 = re$slope.lo, x1 = re$slope.up, y0 = 1:nrow(re), lwd = 1.5)
axis(1, at = seq(-0.4, 0.4, 0.1),  line = 0, las  = 1, cex = 1.1)
# axis(2, at = 1:nrow(re),  line = 0, las  = 1, cex = 1.1, labels = rownames(re), font = 3)
mtext(side = 1, "Climate suitability (slope)", line = 2.6, cex = 1.3)

dev.off()



###################################################################
############################ END ##################################
###################################################################
