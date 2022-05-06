#### Growth and suitability in Eastern of North America
# This code cleans, estimates, and exports the growth and suitability dataset
# Authors: Bernal, Zuleta, Feeley

rm(list = ls())

######################################
################ PATHS ###############
######################################

if(grepl("/Users/danielzuleta", getwd())) { 
  path.to.data <- "//Users/danielzuleta/Google Drive/Coauthoring papers/1. ITRDB Project-Metanalisis/growth-suitability-DZ/data/"
  path.to.climate.data <- paste0(path.to.data, "CRUCLIM-30year-average-Bioclimate/CRUCLIM/")
  path.to.results <- "/Users/danielzuleta/Google Drive/Coauthoring papers/1. ITRDB Project-Metanalisis/growth-suitability-DZ/results/"
}


########################################################
############ LOAD PACKAGES AND FUNCTIONALITY ###########
########################################################
library(readxl)
library(dplR) 
library(sp)
library(FedData)
library(rgbif)
library(raster)
library(dismo)
library(rJava)
library(plyr)
options(java.parameters = "-Xmx1g" )


fn <- function(x) as.numeric(as.character(x)) # force numeric

tr <- function(color, alpha=0.5)
{
  # add transparency to a given color
  rgb <- as.vector(col2rgb(color))
  return(rgb(rgb[1], rgb[2], rgb[3], alpha=alpha*255, maxColorValue=255))
}

######################################
############### LOAD DATA ############
######################################

# loads tree ring database
load(paste0(path.to.data, "treeringdatabase2.RData"))

## load names.sp
names.sp <- read.csv(paste0(path.to.data, "names_sp.csv"))

########################################
######## GROWTH DATA CLEANING ##########
########################################

# first thing first - change factor to character variables
str(treeringdatabase)

treeringdatabase$plot <- as.character(treeringdatabase$plot)
treeringdatabase$tree <- as.character(treeringdatabase$tree)
treeringdatabase$species.code <- as.character(treeringdatabase$species.code)

####clean up DBH data
w=which(treeringdatabase$species.code=="PIST")
t=tapply(treeringdatabase$DBH[w], treeringdatabase$tree[w], max, na.rm=T)
n=names(t[t>=900])
w2=which(!is.na(match(treeringdatabase$tree, n)))
treeringdatabase$DBH[w2]=treeringdatabase$DBH[w2]/10

w=which(treeringdatabase$species.code=="QUMA")
t=tapply(treeringdatabase$DBH[w], treeringdatabase$tree[w], max, na.rm=T)
n=names(t[t>=750])
w2=which(!is.na(match(treeringdatabase$tree, n)))
treeringdatabase$DBH[w2]=treeringdatabase$DBH[w2]/10

w=which(treeringdatabase$species.code=="PIPA")
t=tapply(treeringdatabase$DBH[w], treeringdatabase$tree[w], max, na.rm=T)
n=names(t[t>=1000])
w2=which(!is.na(match(treeringdatabase$tree, n)))
treeringdatabase$DBH[w2]=treeringdatabase$DBH[w2]/10

w=which(treeringdatabase$species.code=="PCGL")
t=tapply(treeringdatabase$DBH[w], treeringdatabase$tree[w], max, na.rm=T)
n=names(t[t>=800])
w2=which(!is.na(match(treeringdatabase$tree, n)))
treeringdatabase$DBH[w2]=treeringdatabase$DBH[w2]/10

w=which(treeringdatabase$plot=="cana269")
treeringdatabase=treeringdatabase[-w,]

# to eliminate QUSP from database
w=which(treeringdatabase$species.code=="QUSP")
treeringdatabase=treeringdatabase[-w,]

# to eliminate PISP from database
w=which(treeringdatabase$species.code=="PISP")
treeringdatabase=treeringdatabase[-w,]

# to eliminate PIBA from database
w=which(treeringdatabase$species.code=="PIBA")
treeringdatabase=treeringdatabase[-w,]

# to eliminate CADN from database
w=which(treeringdatabase$species.code=="CADN")
treeringdatabase=treeringdatabase[-w,]

# to eliminate TSCR from database - few data points (only one plot)
w=which(treeringdatabase$species.code=="TSCR")
treeringdatabase=treeringdatabase[-w,]

# to eliminate PIPO from database - It is from the west od North America
w=which(treeringdatabase$species.code=="PIPO")
treeringdatabase=treeringdatabase[-w,]

#### change species.code Q. prinus (QUPR) to Q. michauxii (QUMI)
treeringdatabase[treeringdatabase$species.code=='QUPR',]$species.code <- 'QUMI'

# Merge treeringdatabase and names.sp
treeringdatabase=merge(treeringdatabase,names.sp,by.x = c("species.code"),by.y = c("Code"),all.x = TRUE)

# Order alphabetically by species name
treeringdatabase <- treeringdatabase[order(treeringdatabase$gensp),]

#### calculate Basal Area based on cleaned DBH values
treeringdatabase$BA=pi*(treeringdatabase$DBH/2)^2 

splist=unique(treeringdatabase$species.code)   ##### create species list

#### clean BA values with negative growth and NA's
# removes individuals with negative growth and NA's
BA.rw=treeringdatabase$BA.rw
u=unique(treeringdatabase$tree)
for(i in 1:length(u)){
  w=which(treeringdatabase$tree==u[i])
  o=order(treeringdatabase$year[w])
  BA=treeringdatabase$BA[w][o]
  age=treeringdatabase$age[w][o]
  BA.rw[w][o]=diff(c(0,BA))/diff(c(0,age))
  if(i %% 1000 == 0) print(i)
}
w=which(BA.rw<=0)
BA.rw[w]=NA
treeringdatabase$BA.rw=BA.rw


###############################################
############ GROWTH STANDARDIZATION  ###########
###############################################

####smooth curve plots per species  ##use "supsmu" supersmoother??
splist=unique(treeringdatabase$species.code)   ##### create species list
species.av.l5=list()

spline.parameters.l5 = data.frame()


# removes NA
w=which(is.na(treeringdatabase$BA.rw))
treeringdatabase=treeringdatabase[-w,]

# removes Inf
w=which(!is.finite(treeringdatabase$BA.rw))
treeringdatabase=treeringdatabase[-w,]


library(randomcoloR)
# Based on pre-processing, we excluded additional trees with weird patterns between the log(Ba.rw) ~ BA.
# this are the trees to be excluded
trees.to.exclude <- read.csv(paste0(path.to.results, "Sp_Plot_Ind_to exclude_20210311.csv"))
trees.to.exclude$sp.plot.tree <- paste0(trees.to.exclude$species,"_",trees.to.exclude$plot,"_",trees.to.exclude$ind)

# id
treeringdatabase$sp.plot.tree <- paste0(treeringdatabase$species.code,"_",treeringdatabase$plot,"_",treeringdatabase$tree)
setdiff(trees.to.exclude$sp.plot.tree, treeringdatabase$sp.plot.tree)
intersect(trees.to.exclude$sp.plot.tree, treeringdatabase$sp.plot.tree)
to.exclude <- which(treeringdatabase$sp.plot.tree %in% trees.to.exclude$sp.plot.tree)
# exclude trees with outliers
treeringdatabase <- treeringdatabase[-to.exclude,]
# exclude trees smaller than 10 cm in dbh
(nrow(treeringdatabase[treeringdatabase$DBH>=100,])/nrow(treeringdatabase))*100 # 64% of the entries is for trees with DBH >= 10cm
treeringdatabase <- treeringdatabase[treeringdatabase$DBH>=100,]

#### Data characteristics to report in the methods/result section
# -------
dim(treeringdatabase)
length(unique(treeringdatabase$species.code))
length(unique(treeringdatabase$plot))
length(unique(treeringdatabase$tree))
length(unique(treeringdatabase$sp.plot.tree))

# export clean treeringdatabase to be used in the plots with predictions
saveRDS(treeringdatabase, paste0(path.to.results, "treeringdatabase.rds"))

# this code generates results from the spline ####
for(i in 1:length(splist)){
  
  # i = which(splist == 'BELE')
  w=which(treeringdatabase$species.code==splist[i]) ##### find data for species i
  # plots=unique(treeringdatabase$plot[w]) ##### find sites where species i occurs
  species.av[[i]] <- species.av.l2[[i]] <- species.av.l3[[i]] <- 
    species.av.l4[[i]] <- species.av.l5[[i]] <- matrix(NA, nrow=length(1750:2010), ncol=length(plots)) ##### create empty matrix to hold species average results - timexsite
  
  for(j in 1:length(plots)){
    # j = which(plots == 'ny023')
    w=which(treeringdatabase$plot==plots[j] & treeringdatabase$year>1750) ##### find the data for site j after 1750
    
    if(length(w)>100){ ##### if there are >100 treerings in site j...
      
      ##### fit spline to log basal area ring width vs basal area with different parameters
      l5=smooth.spline(log(treeringdatabase$BA.rw[w])~treeringdatabase$BA[w], spar = 1.2) 
      # lxxx=ffcsaps(y = log(treeringdatabase$BA.rw[w]), f = 0.5)
      
      # save the parameters
      spline.parameters.l5.i.j <- data.frame(species = splist[i], plot = plots[j],
                                         tol = l5$tol, cv.crit = l5$cv.crit, pen.crit = l5$pen.crit,
                                         crit = l5$crit, df = l5$df, spar = l5$spar, ratio = l5$ratio,
                                         lambda = l5$lambda, nk = l5$fit$nk)
     
      ##### calculate residual of log observed vs log predicted basal area ring width
      fitted.ba.l5 = predict(l5, treeringdatabase$BA[w])[[2]]
      resid.ba.l5 = log(treeringdatabase$BA.rw[w]) - fitted.ba.l5 
      
      # # plot spline
      # plot(treeringdatabase$BA[w], log(treeringdatabase$BA.rw[w]), pch = 16, cex = 0.4, col = tr("grey40"),
      #      xlab = "BA", ylab = "log(BA.rw)", main = paste0(splist[i], "; ", plots[j], "; spar = 1.2"), las = 1)
      # lines(l5, col = "darkorchid3")
      
      ##### get years corresponding to ring width residuals
      x=treeringdatabase$year[w] 
      
      resid.ba.year.l5=tapply(resid.ba.l5, x, mean, na.rm=T)
      m.l5 = match(1750:2010, names(resid.ba.year.l5))
      species.av.l5[[i]][,j]=resid.ba.year.l5[m.l5]
      
      # save spline.parameters
      spline.parameters.l5 <- rbind(spline.parameters.l5, spline.parameters.l5.i.j) 
      
    } ##### close j loop
  }
  cat(paste0("i" = i, " species = ", splist[i]))
} ##### close i loop
# save spline.parameters
write.csv(spline.parameters.l5, paste0(path.to.results, "spline.parameters.l5.csv"))

# Set up growth data same format as climatic data
yr=seq(1901, 1981, by=5)
yr2=NA
for(i in 1:length(yr)){
  w=which(seq(1750,2010,by=1)>=yr[i] & seq(1750,2010,by=1)<(yr[i]+30))
  yr2[w]=yr[i]+14
}

gr.year.l5 <- gr.year.l5.MEDIAN <- list()
for(i in 1:length(species.av.l2)){ # the dimensions in the loop are the same for l2, l3, l4, or l5
  if(length(species.av.l2[[i]])>0){
    gr.year.l5[[i]]=matrix(NA, nrow=length(yr), ncol=dim(species.av.l5[[i]])[2])
    gr.year.l5.MEDIAN[[i]]=matrix(NA, nrow=length(yr), ncol=dim(species.av.l5[[i]])[2])
    for(j in 1:dim(species.av.l2[[i]])[2]){
      t.l5 = tapply(species.av.l5[[i]][,j], yr2, mean, na.rm=T)
      t.l5.MEDIAN = tapply(species.av.l5[[i]][,j], yr2, median, na.rm=T)
      m=match(yr+14, names(t.l2)) # is the same for l2, l3, l4, or l5
      
      gr.year.l5[[i]][,j] = t.l5[m]
      gr.year.l5.MEDIAN[[i]][,j] = t.l5.MEDIAN[m]
    }
  }
}

########################################################################
############ OBTAIN SPECIES SUITABILITY AT THE PLOT LOCATION ###########
########################################################################

m=match(splist, names.sp$Code)
splist2=names.sp$gensp[m] # long species name
years=seq(1915,1995, by=5)

# get the suitability: probability of occurrence of a species in a site 
# based on the predictions from the best niche model for each species
# in each time 
suit.m1915 = list()
for(i in 1:length(splist2)){
  cat("species i = ", i, "-->", splist2[i], "\n")
  if(!is.na(splist2[i])){
    w=which(treeringdatabase$gensp==splist2[i])
    p=unique(treeringdatabase$plot[w])
    m=match(p, treeringdatabase$plot)
    x=data.frame(treeringdatabase$Longitude[m], treeringdatabase$Latitude[m])
    suit.m1915[[i]] =  matrix(NA, ncol=length(x[,1]), nrow=17)
    for(j in 1:length(years)){
      cat("   year j = ", j, "-->", years[j], "\n")
      
      # loads best predicted species distribution model for species i in year j
      p.m1915 = raster(paste0(path.to.results, "MIAme.pred.logistic.m1915/", splist2[i],years[j], ".asc"))
      # get suitability for each species year in EACH PLOT c
      suit.m1915[[i]][j,] = extract(p.m1915,x)
    }
  }
}

######################################################
########## MERGE ALL DATA IN A DATAFRAME #############
######################################################

merging_function<- function(localDF, suit){
  
  localSP <- parent.frame()$i[] #Get the position in the list/
  print(paste(i, j))
  print(localSP)
  localSuitability<-suit[[localSP]] # Get the suitability data using the position from the list.
  # print(localSuitability)
  w=which(treeringdatabase$gensp==splist2[localSP]) ##### find data for species i
  plots=unique(treeringdatabase$plot[w]) ##### find sites where species i occurs
  
  years<-seq(1915,1995, by=5)
  
  resultingDF<-cbind(localDF[,1], localSuitability[,1],  as.character(plots[1]), years)  # Get the first column and pair it with the Plot and the Years.
  # If there are more tha one columns, loop through them to add the plot names and the years.
  if(ncol(localDF)>1)
  {
    for(i in 2:ncol(localDF)){
      newPlot<- cbind(localDF[,i], localSuitability[,i], as.character(plots[i]),years)
      resultingDF<-rbind(resultingDF,newPlot)
    }
  }
  
  resultingDF <- cbind(resultingDF, splist2[localSP]) #Add the species.
  resultingDF
}

retorno.l5.m1915 <- lapply(gr.year.l5, merging_function, suit = suit.m1915)
retorno.l5.m1915.MEDIAN <- lapply(gr.year.l5.MEDIAN, merging_function, suit = suit.m1915)

crecimiento.l5.m1915 <-as.data.frame(do.call(rbind, retorno.l5.m1915))
crecimiento.l5.m1915.MEDIAN <-as.data.frame(do.call(rbind, retorno.l5.m1915.MEDIAN))

colnames(crecimiento.l5.m1915)<- c('growth.l5','suitability.m1915','plot','year','species')
colnames(crecimiento.l5.m1915.MEDIAN)<- c('growth.l5.MEDIAN','suitability.m1915','plot','year','species')

crecimiento.l5.m1915$plot.year.species <- paste(crecimiento.l5.m1915$plot, crecimiento.l5.m1915$year, crecimiento.l5.m1915$species, sep = "_")
crecimiento.l5.m1915.MEDIAN$plot.year.species <- paste(crecimiento.l5.m1915.MEDIAN$plot, crecimiento.l5.m1915.MEDIAN$year, crecimiento.l5.m1915.MEDIAN$species, sep = "_")

# creates crecimiento with all the spline types (l2,..,l5) and the two suit predictions
crecimiento <- crecimiento.l5.m1915
crecimiento$growth.l5.MEDIAN <- crecimiento.l5.m1915.MEDIAN$growth.l5.MEDIAN[match(crecimiento$plot.year.species, crecimiento.l5.m1915.MEDIAN$plot.year.species)]

# Exports
write.csv(crecimiento, paste0(path.to.data, "crecimiento-aic.csv"))




##########################################################
###################      END      ########################
##########################################################