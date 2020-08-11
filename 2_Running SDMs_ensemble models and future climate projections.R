#Create raster stack of climatic variables (preditors) 
#As mentioned this can be sourced from the package and 
#raster functions can be used to clip them the your area
#of interest

bioclims <- list.files(path="E://GisModelling/ClimateChange_input/Current_1970-2000/tiffs/", pattern = "tif$", full.names = TRUE)
predictors <- raster::stack(bioclims)
#Check the predictors included in the stack:
names(predictors)
#Check presence data file:
head(ibis_fil1)
#Visualise the spatial representation of presence data with any one of the predictors:
plot(predictors, 4)
plot(countries, add=TRUE) + points(ibis_fil1, col='black')

#Assess autocorrelation in response and predictor datasets####
#The blockCV package allows for assessment of the autocorrelation range 
#of predictors using spatial blocks, variograms and random points.
#alternatively pearson_correlation_matrix within the sdmpredictors package can be used.

library(blockCV)
library(sdmpredictors)
library(automap)

auto_predict <- spatialAutoRange(rasterLayer = predictors,sampleNumber = 5000,
                        doParallel = TRUE, showPlots = TRUE)
#assess variogram outputs
summary(auto_predict)
plot(auto_predict$variograms[[1]])

#pearson_correlation_matrix(predictors, cachesize = 20, same_mask = FALSE)

#remove any highly autocorrelated layers (predictors) if needed:
predictors <- dropLayer(predictors, c(1, 3, 4, 6, 9, 11))
#Check predictors remaining in the raster stack:
names(predictors)
plot(predictors)

#Load response variables into workspace
#Merge presence data and backgorund points created earlier into a single file
library(raster)
library(sf)
backpnts <- raster::extract(predictors, background_ibis) 
prespnts <- raster::extract(predictors, presence_ibis) 
pb <- c(rep(1, nrow(prespnts)), rep(0, nrow(backpnts)))
sdmdata <- data.frame(cbind(pb, rbind(prespnts, backpnts)))
str(sdmdata)

#Check number of presence and backgorud points
table(sdmdata$pb)

#visually investigate colinearity in dataset to see if the blockCv autocorrelation 
#excerise was efficient at addressing this
pairs(sdmdata[,2:5], cex=0.1, fig=TRUE)

###Spatial blocking can be used to further control for autocorrelation.
#Spatial blocl grids can be used within both Maxent and Biomod2 modelling frameworks
#Conduct spatial blocking (See blockCv for further details): 
spatialblock1 <- spatialBlock(speciesData = pb_data, # sysematic blocking
                    species = "Species",
                    rasterLayer = predictors,
                    rows = 8,
                    cols = 10,
                    k = 5,
                    selection = "systematic",
                    biomod2Format = TRUE)

spatialblock2 <- spatialBlock(speciesData = pa_data, # checkerboard blocking
                    species = "Species",
                    rasterLayer = predictors,
                    rows = 6,
                    selection = "checkerboard",
                    biomod2Format = TRUE)

enviroblock <- envBlock(rasterLayer = predictors, #Environemental blocking
               speciesData = pb_data,
               species = "Species",
               k = 10,
               standardization = "standard", # rescale variables between 0 and 1
               rasterBlock = FALSE,
               numLimit = 50)

#Alternatively, you can manually scale through block sizes:
rangeExplorer(rasterLayer = predictors,
              speciesData = pb_data, # response data (optional)
              species = "Species", # the responcse column (optional)
              minRange = 30000, # limit the search domain
              maxRange = 100000)

#Building SDM models####
#Using Biomod2 ensemble framework
library(biomod2)
#Merge presence and background DFs created
colnames(presence_ibis) <- c("x", "y")
sbi <- c(rep(1, nrow(presence_ibis)), rep(0, nrow(background_ibis)))
sdmdata <- data.frame(cbind(sbi, rbind(presence_ibis, background_ibis)))
#Prep data format for use in Biomod2
myRespName <- 'sbi'
myResp <- as.numeric(sdmdata[,myRespName])
myRespXY <- sdmdata[,c("x","y")]
myExpl = predictors
BiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)
BiomodData

#View presence-background data
plot(BiomodData)

#Definfing the folds to be used in model runs: 
#use blockCV output here
DataSplitTable <- enviroblock$biomodTable

#running multiple algorithms to build SDMs 
BiomodOption <- BIOMOD_ModelingOptions()
BiomodModelOut <- BIOMOD_Modeling(
  BiomodData,
  models = c('MAXENT.Phillips','GBM','RF'),   #Select models
  models.options = BiomodOption,
  NbRunEval=3,
  DataSplitTable = DataSplitTable, # blocking folds
  DataSplit=75,
  Prevalence=0.5,
  VarImport=3,
  models.eval.meth = c('TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = TRUE,
  do.full.models = FALSE,
  modeling.id = paste(myRespName,"FirstModeling",sep=""))

#View model output summary:
BiomodModelOut

#Model evlauations
BiomodModelEval <- get_evaluations(BiomodModelOut)
BiomodModelEval["ROC","Testing.data",,,]
get_variables_importance(BiomodModelOut)

#Ensemble modelling - creating an ensemble from the indivual model outputs above:
BiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = BiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

#Evlauate emsemble output:
BiomodEM
get_evaluations(BiomodEM)

#Best practice should include ca. 10% data points withheld for 
#model validation. This can be run in the following script.

#Future climate change projections####
#Using the above models built for the species, we can project it from current
#climatic conditions (i.e. "predictors" raster stack) to the same set of climatic 
#varibales simulated under various climate change scenarios (RCPs and GCMs)

#Create raster stack representing furture scenarios
#Although this can also be sourced through the package, newer simulations
#are available from WorldClim directly
bioclims_fut <- list.files(path="E://GisModelling/ClimateChange_input/Future_2020-2040/tiffs/", pattern = "tif$", full.names = TRUE)
predictors_fut <- raster::stack(bioclims_fut)

BiomodProjFuture <- BIOMOD_Projection(
  modeling.output = BiomodModelOut,
  new.env = predictors_fut,
  proj.name = 'GCM_future_2040',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

##The BIOMOD_EnsembleForecasting function can be used to created forecasted plots predicted suitability
#from current to projected future conditions.




