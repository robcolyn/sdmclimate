#Assessing thresholds to convert models from continous to binary ouputs
#Additionally allows for asessment of model fit - AUc, ROC plots, etc. 

library(PresenceAbsence)

#Extract predicted suitability values from SDM for both presence and absence points
#A better approach here is to limit backgorund points to known areas of absence. Using
#the circular window method (see SDm script) is a good method of generating background 
#points in the absence of true absence data.
evaldata <- raster::extract(BiomodEM, sdmdata) 

#Create histogram of predicted suitability to that of presence and absence localities
par(mfrow=c(1,1))
for(i in 1:3)
{presence.absence.hist(DATA = evaldata, na.rm = TRUE, which.model = 1, model.names=c("Southern Bald Ibis"), 
                       N.bars=10, truncate.tallest=TRUE, opt.thresholds=TRUE, xlab = "Predicted probability",
                       ylab = "Number of plots", opt.methods=c("PredPrev=Obs","Sens=Spec","MaxKappa"))}

#Create error plot with various thresholds displayed
error.threshold.plot(DATA = evaldata, na.rm = TRUE, which.model = 1, model.names=c("Southern Bald Ibis"),
                     opt.thresholds=TRUE, opt.methods=c("PredPrev=Obs","Sens=Spec","MaxKappa"), vert.lines=FALSE,
                     color = TRUE)

#Create a AUC ROC plot assessing model fit and threshold selection
auc.roc.plot(DATA = evaldata, color=TRUE, main = "Southern Bald Ibis")

#Replicate the above graphs using the dataset withheld for model
#validation. Subistute presence points to only include validation
#data and assess AUC-ROC plots and histograms of predicted suitability. 

#The threshold selected should ideally limit ommision rates of presence 
#data, produce lower error rates, and higher AUC values. 