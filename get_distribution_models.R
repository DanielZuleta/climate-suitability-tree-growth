# 1-Niche models with 19 climatic variables.
# Bernal, Zuleta, Feeley


rm(list = ls()) 

######################################
################ PATHS ###############
######################################

if(grepl("/My Drive/My Research/", getwd())) { 
  path.to.data <- "G:/My Drive/My Research/1. ITRDB Project-Metanalisis/growth-suitability-DZ/data/"
  path.to.climate.data <- "G:/My Drive/My Research/1. ITRDB Project-Metanalisis/CRUCLIM-30year-average-Bioclimate/CRUCLIM/"
  path.to.results <- "G:/My Drive/My Research/1. ITRDB Project-Metanalisis/"
}

if(grepl("/Users/zuletad", getwd())) { 
  path.to.data <- "//Users/zuletad/Google Drive/Coauthoring papers/1. ITRDB Project-Metanalisis/growth-suitability-DZ/data/"
  path.to.climate.data <- paste0(path.to.data, "CRUCLIM-30year-average-Bioclimate/CRUCLIM/")
  path.to.results <- paste0("//Users/zuletad/Google Drive/Coauthoring papers/1. ITRDB Project-Metanalisis/growth-suitability-DZ/results/")
}


########################################################
############ LOAD PACKAGES AND FUNCTIONALITY ###########
########################################################

# remotes::install_github("johnbaums/rmaxent")
library(remotes)
library(rmaxent)
require(sp)
require(BIEN)
require(raster)
require(rgdal)
require(dismo)
require(ENMeval)
# remotes::install_github("julienvollering/MIAmaxent", force = T)
require(MIAmaxent)
require(rJava)      # Sys.setenv(JAVA_HOME="")
options(java.parameters = "-Xmx1g" )

set.seed(43234248)


get_transformed_variables <- function(x, transformations, model.selected ){
  
  # which variables were selected?
  variables.selected <- names(model.selected$coefficients[-1])
  main.variables <- unique(sub("\\_.*", "", variables.selected))
  
  ### getting transformed vars
  x.var = x[,-1]
  transformed.data <- matrix(NA, nrow = nrow(x.var), ncol = length(transformations) - 1)
  z = 1
  acc.trans.var <- c()
  for(i in 1:ncol(x.var)){
    var <- names(x.var)[i]
    data.var.i <- x.var[,var]
    trans.var <- transformations[grepl(var, names(transformations))]
    
    for(j in 1:length(trans.var)){
      trans.j.var.i <- trans.var[j]
      data.var.i.trans.j <- trans.j.var.i[[1]](data.var.i)
      transformed.data[,z] <- data.var.i.trans.j
      z = z + 1
    } # closes loop for transformations of variable i
    acc.trans.var <- c(acc.trans.var, names(trans.var))
  }  # closees loop for each variable
  
  acc.trans.var <- substr(acc.trans.var, 0, nchar(acc.trans.var)-7)
  colnames(transformed.data) <- acc.trans.var
  all <- cbind(x, transformed.data)
  
  ### convert to lists 
  all.selected.variables <- cbind(RV = all[,"RV"], all[,names(all) %in% variables.selected])
  
  result <- list(RV = all.selected.variables[,"RV"])
  for(mv in 1:length(main.variables)){
    var <- main.variables[mv]
    
    temp <- data.frame(all.selected.variables[,grepl(var, names(all.selected.variables))])
    names(temp) <- names(all.selected.variables)[grepl(var, names(all.selected.variables))]
    result[[mv+1]] <- temp
  }
  names(result) <- c("RV", main.variables)
  
  attr(result, "main.variables") <- main.variables
  return(result)
}


##################################
############ LOAD DATA ###########
##################################

# Loads occurrence data
bien.all.clean <- read.csv(paste0(path.to.data, "occ.bien.clean.csv")) # from BIEN
bien.all.clean <- bien.all.clean[,c('species', 'latitude', 'longitude')] # deleting first two empty columns
gbif.all.clean <- read.csv(paste0(path.to.data, "occ.gbif.clean.csv")) # from GBIF
# merging both data sources gbif and bien
itrdb.all.clean <- rbind(bien.all.clean, gbif.all.clean)

# making it congruent with MIAmaxent format
itrdb.all.clean$Occurrence <- "presence"
# the order
itrdb.all.clean <- itrdb.all.clean[, c(4,3,2,1)]
# Order alphabetically by species name
itrdb.all.clean <- itrdb.all.clean[order(itrdb.all.clean$species),]  

spdata <- itrdb.all.clean
rm(itrdb.all.clean)

# remove species TSCA and PIPO: out of range and low occurrence points
spdata <- spdata[!(spdata$species %in% c('Tsuga caroliniana', 'Pinus ponderosa')),]

years=seq(1915,1995, by=5)  
# loads variables to subsequently thin the data (one for each pixel)
for(w in 1:length(years)){
  MAT=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_01_", years[w], "North.America.asc", sep=""))
  MDR=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_02_", years[w], "North.America.asc", sep=""))
  Isothermality=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_03_", years[w], "North.America.asc", sep=""))
  Tseasonality=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_04_", years[w], "North.America.asc", sep=""))
  Maxtemp=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_05_", years[w], "North.America.asc", sep=""))
  Mintemp=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_06_", years[w], "North.America.asc", sep=""))
  Trange=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_07_", years[w], "North.America.asc", sep=""))
  MeanTempWettestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_08_", years[w], "North.America.asc", sep=""))
  MeanTempDriestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_09_", years[w], "North.America.asc", sep=""))
  MeanTempWarmestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_10_", years[w], "North.America.asc", sep=""))
  MeanTempColdestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_11_", years[w], "North.America.asc", sep=""))
  TAP=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_12_", years[w], "North.America.asc", sep=""))
  PrecWettestMonth=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_13_", years[w], "North.America.asc", sep=""))
  PrecDriestMonth=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_14_", years[w], "North.America.asc", sep=""))
  pseasonality=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_15_", years[w], "North.America.asc", sep=""))
  PrecWettestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_16_", years[w], "North.America.asc", sep=""))
  PrecDriestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_17_", years[w], "North.America.asc", sep=""))
  PrecWarmestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_18_", years[w], "North.America.asc", sep=""))
  PrecColdestQ=raster(paste(path.to.climate.data, years[w], "/CRUCLIM_19_", years[w], "North.America.asc", sep=""))
  predictors=stack(MAT, MDR, Isothermality, Tseasonality, Maxtemp, Mintemp,  Trange, MeanTempWettestQ, MeanTempDriestQ, MeanTempWarmestQ, MeanTempColdestQ, 
                   TAP, PrecWettestMonth, PrecDriestMonth, pseasonality, PrecWettestQ, PrecDriestQ, PrecWarmestQ,  PrecColdestQ)
  names(predictors)=c("MAT", "MDR", "Isothermality", "Tseasonality", "Maxtemp", "Mintemp",  "Trange", "MeanTempWettestQ", "MeanTempDriestQ", 
                      "MeanTempWarmestQ", "MeanTempColdestQ", "TAP", "PrecWettestMonth", "PrecDriestMonth", "pseasonality", "PrecWettestQ", "PrecDriestQ", 
                      "PrecWarmestQ",  "PrecColdestQ")
  projection(predictors)= "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
}

#clean and thin species occurrences data
coordinates(spdata) <- ~longitude+latitude
cell <- cellFromXY(predictors, spdata)
d=data.frame(cell, spdata$species)
dup <- duplicated(d)
spdata2 <- spdata[!dup, ]
predvalues=extract(predictors[[1]], spdata2)
w=which(!is.na(predvalues))
spdata2=spdata2[w,]

# 99.6 % of records were removed after cleaning
round(100-(nrow(spdata2)/nrow(spdata))*100, 1)

# Converting SpatialPoints to csv
spdata2.df <- data.frame(longitude=coordinates(spdata2)[,1], latitude=coordinates(spdata2)[,2], spdata2@data)
splist=unique(spdata2.df$species)
saveRDS(spdata2.df, paste0(path.to.results, "spdata2.df.rds"))

pdf(paste0(path.to.results, "species.occurrence.points.after.filtering.pdf"), width = 7, height = 6)
par(oma = c(0,0,0,0), mar = c(4,8,3,2))
barplot(table(spdata2.df$species)[order(table(spdata2.df$species))], las = 1, cex.names = 0.7, horiz = T, xlab = "Occurrences")
dev.off()

plot(spdata2.df$longitude, spdata2.df$latitude, cex = 0.2, pch  = 19)




##### DISTRIBUTION MODEL ######

# creates one folder to save the occurences for each species - this is the format needed in MIAmaxent
for(i in 1:length(splist)){
  w = which(spdata2.df$species == splist[i])
  # but first eliminates files in case they were previously created
  unlink(paste0(path.to.data, "occurrences/", splist[i],"/"), recursive = T)
  # creates milindroso
  occ.sp = spdata2.df[w,]
  occ.sp.milindroso <- occ.sp[,c(3,1,2)]
  dir.create(paste0(path.to.data, "occurrences/", splist[i],"/"))
  write.csv(occ.sp.milindroso, paste0(path.to.data, "occurrences/", splist[i],"/", "occ.",splist[i],".csv"),row.names = FALSE)
}

# Maxent models and environmental variable selection using cross-validation with k = 5 replicates

# Estimate the best model for the First period (1915 Oldest climate) 

# create empty lists... they are placeholders for the final model output
best.models.1915 = list() # the final list, each element representing one species
best.raw.preds.based.on.1915 = list()## the final list, each element representing one species (each element containing 17 lists)
best.preds.based.on.1915 = list()# the final list, each element representing one species (each element containing 17 lists)
# the obnly two temporal objects
raw.preds = list() #will be a list of 17 for each species
mod.preds = list()#will be a list of 17 for each species


for(i in 1:length(splist)){
  print(paste0(splist[i], " (i=",i,")"))
  w = which(spdata2.df$species == splist[i])
  occ.sp = spdata2.df[w,]  #here paste in your code to produce the occ.sp object (which is the occurrences for species i)
  
  if(dim(occ.sp)[1]>=100){
    
    data.spi.1915.FULL <- readData(
      occurrence = paste0(path.to.data, "occurrences/", splist[i],"/", "occ.",splist[i],".csv"),
      contEV = paste0(path.to.climate.data, '1915'),
      maxbkg=10000, XY = T)
    
    # Since the model selection is only forward we changed a bit the order:
    new.order.EVS.1915 <- paste0(c("CRUCLIM.01", # MAT
                                   "CRUCLIM.12", # TAP
                                   "CRUCLIM.04", # Temp seasonality
                                   "CRUCLIM.15", # Precip seasonality
                                   "CRUCLIM.02", "CRUCLIM.05",
                                   "CRUCLIM.06", "CRUCLIM.13", "CRUCLIM.14",
                                   "CRUCLIM.03", "CRUCLIM.07", "CRUCLIM.08",
                                   "CRUCLIM.09", "CRUCLIM.10", "CRUCLIM.11",
                                   "CRUCLIM.16", "CRUCLIM.17", "CRUCLIM.18",
                                   "CRUCLIM.19"), ".",'1915', "North.America")
    data.spi.1915 <- data.spi.1915.FULL[,c( "RV", new.order.EVS.1915)]
    
    # preparing env variables to predict
    EVfiles.1915 <- c(list.files(paste0(path.to.climate.data, '1915'), 
                                 full.names=TRUE, pattern = "\\.asc$"))
    EVstack.1915 <- raster::stack(EVfiles.1915)
    names(EVstack.1915) <- gsub(".asc", "", basename(EVfiles.1915))
    projection(EVstack.1915) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
    
    # model selection
    w.background <- is.na(data.spi.1915$RV)
    occ.all <- data.spi.1915[!w.background,]
    occ.full.xy <- data.spi.1915.FULL[!w.background,c('x','y')]
    background.all <- data.spi.1915[w.background,]
    aic.val.1915 <- c(rep(NA, 5)) #creates an empty vector for the AICc 
    models.1915 <- transformations <- list()
    # set sample sizes of occ and background
    smp_size.occ.all <- floor(0.75 * nrow(occ.all))
    smp_size.background.all <- floor(0.75 * nrow(background.all))
    for(k in 1:5){
      # set the seed to make the partition reproducible
      set.seed(123*k)
      # train sample
      train.occ.k.ids <- sample(seq_len(nrow(occ.all)), size = smp_size.occ.all)
      train.background.k.ids <- sample(seq_len(nrow(background.all)), size = smp_size.background.all)
      # train dataset
      train.occ.k <- occ.all[train.occ.k.ids, ]
      train.background.k <- background.all[train.background.k.ids, ]
      to.model.k <- data.frame(rbind(train.occ.k,train.background.k))
      
      # Test dataset (observed occrences x and y not used to train the model)
      test.occ.k <- occ.full.xy[-train.occ.k.ids, ]
      
      ### Model selection
      ### Part 1: select the best transformed variables for each explanatory variable
      # Transforming explanatory variables
      spi.1915.DVs <- deriveVars(to.model.k, 
                                 transformtype = c("L","D"))
      # Selecting derived variables
      spi.1915.DVselect <- selectDVforEV(spi.1915.DVs$dvdata, alpha = 0.001, quiet = TRUE)
      # summary(spi.1915.DVselect$dvdata)
      
      ### Part 2: variable selection 
      spi.1915.EVselect <- selectEV(spi.1915.DVselect$dvdata, alpha = 0.001, 
                                    interaction = TRUE, quiet = T)
      # save the model
      models.1915[[k]] <- spi.1915.EVselect$selectedmodel
      
      #### Getting AICc
      spi.1915.Preds <- projectModel(model = spi.1915.EVselect$selectedmodel,
                                     transformations = spi.1915.DVs$transformations,
                                     data = EVstack.1915)
      transformations[[k]] <- spi.1915.DVs$transformations
      
      # GET AICc of the TEST in partition k
      aic.val.1915[k] = calc.aicc(nparam = length(spi.1915.EVselect$selectedmodel$coefficients), 
                                  occ = test.occ.k,
                                  predictive.maps = spi.1915.Preds$output)$AICc
      
      # there are some k for which the model is tooo bad and do not even have variability in the predictions. Then,
      # we are using the following in order to not even consider them (it only happens in 8 cases see growth-suitability-DZ/modelswithoutvariability.xlsx):
      var.all.vals.1915 <- var(spi.1915.Preds$output$PRO@data@values, na.rm = T)
      if(is.na(var.all.vals.1915) | var.all.vals.1915 == 0) aic.val.1915[k] <- NA
    }
    min.AIC.model.1915 <- which.min(aic.val.1915) 
    best.mod.1915 <- models.1915[[min.AIC.model.1915]]
    transformations.best.model.1915 <- transformations[[min.AIC.model.1915]]
    
    # Once the formula of the best model is selected, we use this formula to readjust the best model 
    # using all the data (occurrences and their climatic variables)
    # To do so, we created function get_transformed_variables (see above in Functionality section)
    dv.data.selected.model.1915 <- get_transformed_variables(x = data.spi.1915,
                                                             transformation = transformations.best.model.1915,
                                                             model.selected = best.mod.1915)
    main.variables.1915 <- attr(dv.data.selected.model.1915, "main.variables")
    # 
    formula.temp <- c()
    for (otroi in 1:length(main.variables.1915)) formula.temp <- paste0(formula.temp, main.variables.1915[otroi],ifelse( otroi == length(main.variables.1915),""," + "))
    formula.best.model.1915 <- paste0("~ ", formula.temp)
    readjusted.best.model.1915 <- chooseModel(dvdata = dv.data.selected.model.1915, # 
                                              formula = formula(formula.best.model.1915),  # from the best model
                                              algorithm = "maxent")
    
    # Now, with the best model in 1915, predicts prob of occurrence for the rest of the periods
    for(j in 1:length(years)){
      
      # preparing env variables to predict
      EVfiles.j <- c(list.files(paste0(path.to.climate.data, years[j]), 
                                full.names=TRUE, pattern = "\\.asc$"))
      EVstack.j <- raster::stack(EVfiles.j)
      names(EVstack.j) <- gsub(".asc", "", basename(EVfiles.j))
      projection(EVstack.j) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+towgs84=0,0,0"
      
      # IMPORTANT: we first need to change the names of the variables in the 
      # EV.stack.j because they have the year in the name and those names need to be consistent with the variables
      # in the model objects (readjusted.best.model.1915 and transformations.best.model.1915) for the projectModel function
      names(EVstack.j) <- sub(paste0("_", years[j]), '_1915', names(EVstack.j))
      
      # places the prediction of the best model obtained from 1915 into a list for all the periods
      raw.preds[[j]] = projectModel(model = readjusted.best.model.1915,
                                    transformations = transformations.best.model.1915,
                                    data = EVstack.j,
                                    clamping = TRUE,
                                    raw = TRUE,
                                    rescale = TRUE)$output
      points(occ.sp[,1:2], cex = 0.3, pch = 16)
      mod.preds[[j]] = to_logistic(raw.preds[[j]], "raw", readjusted.best.model.1915[["entropy"]]) # convert raw to logistic prediction
      plot(mod.preds[[j]])
      points(occ.sp[,1:2], cex = 0.3, pch = 16)
    } # closes loop for year j
    
    best.models.1915[[i]] = readjusted.best.model.1915 #places the list of best models per year into another list of each of the species
    best.raw.preds.based.on.1915[[i]] = raw.preds # raw predictions
    best.preds.based.on.1915[[i]] = mod.preds # logistic predictions
  } # closes if n>100
} # closes loop for species i

#save models  and plots of best model predictions
# directories to folders save the files
path.to.best.models.m1915 <- paste0(path.to.results, "MIAmaxent.best.models.m1915/")
path.to.pred.raw.m1915 <- paste0(path.to.results, "MIAme.pred.raw.m1915/")
path.to.pred.logistic.m1915 <- paste0(path.to.results, "MIAme.pred.logistic.m1915/")
path.to.pred.plots.m1915 <- paste0(path.to.results, "species-me-predictions-plots-m1915/")
# first eliminates folders in case they were previously created
unlink(path.to.best.models.m1915, recursive = T)
unlink(path.to.pred.raw.m1915, recursive = T)
unlink(path.to.pred.logistic.m1915, recursive = T)  
unlink(path.to.pred.plots.m1915, recursive = T)
# then, create the empty folders again
dir.create(path.to.best.models.m1915)
dir.create(path.to.pred.raw.m1915)
dir.create(path.to.pred.logistic.m1915)
dir.create(path.to.pred.plots.m1915)
for(i in 1:length(splist)){
  cat("species ",i, "\n")
  # define occurrences of the species to be plotted
  w = which(spdata2.df$species == splist[i])
  occ.sp = spdata2.df[w,]
  
  pdf(paste0(path.to.pred.plots.m1915,  splist[i], ".pdf"), width = 6, height = 5)
  par(oma = c(0,0,0,0), mar = c(4.1,4.1,4.2,3))
  for(j in 1:length(years)){
    cat("year ",j, "\n")
    MIAme.results = best.models.1915[[i]]
    MIAme.pred.best.model.RAW = best.raw.preds.based.on.1915[[i]][[j]]
    MIAme.pred.best.model = best.preds.based.on.1915[[i]][[j]]
    # save
    save(MIAme.results,                 file =    paste0(path.to.best.models.m1915,   splist[i], ".RData"))
    writeRaster(MIAme.pred.best.model.RAW, paste0(path.to.pred.raw.m1915,      splist[i], years[j], ".asc"), format="ascii")
    writeRaster(MIAme.pred.best.model,     paste0(path.to.pred.logistic.m1915, splist[i], years[j], ".asc"), format="ascii")
    # plot best model
    plot(MIAme.pred.best.model, main = paste0(splist[i], "; year = ", years[j]),
         xlab = "Longitude", ylab = "Latitude")
    points(occ.sp[,1:2], cex = 0.3, pch = 16)
  }
  dev.off()
}


############################
############ END ###########
############################


