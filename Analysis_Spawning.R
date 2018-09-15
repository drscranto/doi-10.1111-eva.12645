################
################
## Created: 11/17/16
## Author: Katherine Scranton scranton.kt@gmail.com
## Updated: 04/26/18
## 
## Code for Spawning analysis in:
## Evaluating the potential for pre-zygotic isolation and hybridization between landlocked and anadromous alewife (Alosa pseudoharengus) following secondary contact. 2018. Katherine A. Litrell, David Ellis, Stephen R. Gephard, Andrew D. MacDonald, Eric P. Palkovacs, Katherine Scranton, David M. Post. Evolutionary Applications.
##
################

# load the local source file that calculates model likelihoods
source("TimeToEvent_Lik.R")

# load the local spawning data
sp.data <- read.csv("Alewife_Spawning_Data.csv",stringsAsFactors=FALSE)

# ensure factors are treated correctly
sp.data$Day <- as.numeric(sp.data$Day)
sp.data$Count <- as.numeric(sp.data$Count)
sp.data$Year <- as.numeric(sp.data$Year)
sp.data$Temp <- as.numeric(sp.data$Temp)
sp.data$I.Form <- as.numeric(sp.data$Form=="L")

# calculate likelihoods for all subsets of the model 
# rate parameter ~  exp( ln(rate0) + form + year + lake)
# year was treated as a factor because of the limited sample (3 years for 1 lake, 2 years for other 4 lakes)
# lake was also treated as a fixed effect because of the limited sampling within body form (2 lakes for anadromous fish, 3 for landlocked fish)

fixed.models <- list()
par.inits <- c(1/160,10,0,0,0,0,0,0,0,0,0)
fixed.models[[1]] <- max.lik.fixed.null(sp.data,par.inits)
par.inits <- c(1/140,10,-0.3,0,0,0,0,0,0,0,0)
fixed.models[[2]] <- max.lik.fixed.f(sp.data,par.inits)
par.inits <- c(1/160,10,0,0,0,0,0,0,0,0,0)
fixed.models[[3]] <- max.lik.fixed.y(sp.data,par.inits)
par.inits <- c(1/160,10,0,0,0,0,0,0,0,0,0)
fixed.models[[4]] <- max.lik.fixed.y2013(sp.data,par.inits) 
par.inits <- c(1/160,10,0,0,0,0,0,0,0,0,0)
fixed.models[[5]] <- max.lik.fixed.y2014(sp.data,par.inits)
par.inits <- c(1/160,10,0,0,0,0,0,0,0,0,0)
fixed.models[[6]] <- max.lik.fixed.l(sp.data,par.inits)
par.inits <- c(1/160,10,0,0,0,0,0,0,0,0,0)
fixed.models[[7]] <- max.lik.fixed.lP(sp.data,par.inits) 

par.inits <- c(1/140,10,-0.3,0.06,0.02,0,0,0,0,0,0)
fixed.models[[8]] <- max.lik.fixed.fy(sp.data,par.inits)
par.inits <- c(1/140,10,-0.3,0.06,0,0,0,0,0,0,0)
fixed.models[[9]] <- max.lik.fixed.fy2013(sp.data,par.inits)
par.inits <- c(1/140,10,-0.3,0,0,0,0.23,0.25,0.07,-0.02,0)
fixed.models[[10]] <- max.lik.fixed.fl(sp.data,par.inits)
par.inits <- c(1/140,10,-0.3,0,0,0,0,0,0.07,0,0)
fixed.models[[11]] <- max.lik.fixed.flP(sp.data,par.inits)
par.inits <- c(1/160,10,0,0.06,0.02,0,0.23,0.25,0.07,-0.02,0)
fixed.models[[12]] <- max.lik.fixed.yl(sp.data,par.inits)
par.inits <- c(1/160,10,0,0.06,0,0,0.23,0.25,0.07,-0.02,0)
fixed.models[[13]] <- max.lik.fixed.y2013l(sp.data,par.inits) 
par.inits <- c(1/160,10,0,0.06,0.02,0,0,0,0.07,0,0)
fixed.models[[14]] <- max.lik.fixed.ylP(sp.data,par.inits)
par.inits <- c(1/160,10,0,0.06,0,0,0,0,0.07,0,0)
fixed.models[[15]] <- max.lik.fixed.y2013lP(sp.data,par.inits) 

par.inits <- c(0.007,17,-0.25,0.08,0.002,0,0.02,0.04,0.06,-0.02,0)
fixed.models[[16]] <- max.lik.fixed.fyl(sp.data,par.inits)
par.inits <- c(0.007,17,-0.25,0.08,0,0,0.02,0.04,0.06,-0.02,0)
fixed.models[[17]] <- max.lik.fixed.fy2013l(sp.data,par.inits)
par.inits <- c(0.007,17,-0.25,0.08,0.002,0,0,0,0.06,0,0)
fixed.models[[18]] <- max.lik.fixed.fylP(sp.data,par.inits)
par.inits <- c(0.007,17,-0.25,0.08,0,0,0,0,0.06,0,0)
fixed.models[[19]] <- max.lik.fixed.fy2013lP(sp.data,par.inits)

# make vector of the number of parameters estimated in each model
ks <- 2+c(0,1,2,1,1,4,1,3,2,5,2,6,5,3,2,7,6,4,3)

# using the likelihoods and parameters, calculate aic values for each model
all.aics <-c()
for (i in 1:19){
	all.aics <- c(all.aics,2*ks[i]+2*fixed.models[[i]]$value)
}


