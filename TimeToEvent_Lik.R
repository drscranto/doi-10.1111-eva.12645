################
################
## Created: 11/17/16
## Author: Katherine Scranton scranton.kt@gmail.com
## Updated: 04/26/18
## 
## Code for Spawning analysis in:
## Evaluating the potential for pre-zygotic isolation and hybridization between landlocked and anadromous alewife (Alosa pseudoharengus) following secondary contact. 2018. Katherine A. Litrell, David Ellis, Stephen R. Gephard, Andrew D. MacDonald, Eric P. Palkovacs, Katherine Scranton, David M. Post. Evolutionary Applications.
##
## Likelihood calculations for a fixed effects time to event model of spawning
## 3 fixed effects: form, lake, year
## We can resolve spawning time to the day: interval censored on days
## We use a Weibull distribution because of the flexibility of the hazard
## We assume Normal distributions for the errors
##
################

# main likelihood function
#####################

# calculate likelihood of the data set for a given set of parameter values

log.lik.fixed <- function(sp.data,w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q,beta.2014.P){
	
	all.lakes <- c("Bride","Dodge","Pattagansett","Quonnipaug","Rogers")
	log.lik <- 0
	for (i.lake in 1:length(all.lakes)){
		lake.effect <- c(beta.B,beta.D,beta.P,beta.Q,0)[i.lake]
		lake.year.int <- as.numeric(i.lake==3)
		lake.data <- sp.data[sp.data$Lake == all.lakes[i.lake],]
		lake.years <- unique(lake.data$Year)
		
		lake.lik <- 0
		for (i.year in 1:length(lake.years)){
			
			year.effect <- beta.2013*as.numeric(lake.years[i.year]==2013) + beta.2014*as.numeric(lake.years[i.year]==2014)
			lake.year.int <- beta.2014.P*lake.year.int*as.numeric(lake.years[i.year] == 2014)
			
			ly.data <- lake.data[lake.data$Year == lake.years[i.year],]
			all.effects <- beta.form*ly.data$I.Form + year.effect + beta.temp*ly.data$Temp + lake.effect + lake.year.int
			
			# exact prob
			#(w.shape*exp(log(w.rate)+all.effects))*((ly.data$Day*exp(log(w.rate)+all.effects))^(w.shape-1))*exp(-1*((exp(log(w.rate)+all.effects)*(ly.data$Day))^w.shape))
			## S = exp(-(lambdas*t)^shape)
			## lambdas = exp(lograte + fixed effects + random effects)
			S.L <- exp(-1*((exp(log(w.rate)+all.effects)*(ly.data$Day))^w.shape))
			S.H <- exp(-1*((exp(log(w.rate)+all.effects)*(ly.data$Day+1))^w.shape))
			
			lake.lik <- lake.lik + sum(log((S.L - S.H)^ly.data$Count))

		}
		log.lik <- log.lik + lake.lik
	}
	log.lik
}

# optimize the likelihood function for certain subsets of the model
#####################

# rate parameter ~  exp( ln(rate0) )
max.lik.fixed.null <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,0,0,0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[1:2],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form )
max.lik.fixed.f <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],0,0,0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[1:3],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year )
max.lik.fixed.y <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,x[3],x[4],0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,4,5)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year13 )
max.lik.fixed.y2013 <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,x[3],0,0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,4)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year14 )
max.lik.fixed.y2014 <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,0,x[3],0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,5)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + lake)
max.lik.fixed.l <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,0,0,0,x[3],x[4],x[5],x[6],0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,7,8,9,10)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + lakePattagenset)
max.lik.fixed.lP <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,0,0,0,0,0,x[3],0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,9)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + year )
max.lik.fixed.fy <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],x[4],x[5],0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,4,5)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + year13 )
max.lik.fixed.fy2013 <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],x[4],0,0,0,0,0,0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,4)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + lake)
max.lik.fixed.fl <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],0,0,0,x[4],x[5],x[6],x[7],0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,7,8,9,10)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + lakePattagansett)
max.lik.fixed.flP <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],0,0,0,0,0,x[4],0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,9)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year + lake)
max.lik.fixed.yl <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,x[3],x[4],0,x[5],x[6],x[7],x[8],0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,4,5,7,8,9,10)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year13 + lake)
max.lik.fixed.y2013l <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,x[3],0,0,x[4],x[5],x[6],x[7],0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,4,7,8,9,10)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year + lakePattagansett)
max.lik.fixed.ylP <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,x[3],x[4],0,0,0,x[5],0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,4,5,9)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + year13 + lakePattagansett)
max.lik.fixed.y2013lP <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],0,x[3],0,0,0,0,x[4],0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,4,9)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + year + lake)
max.lik.fixed.fyl <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],x[4],x[5],0,x[6],x[7],x[8],x[9],0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,4,5,7:10)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + year13 + lake)
max.lik.fixed.fy2013l <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],x[4],0,0,x[5],x[6],x[7],x[8],0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,4,7:10)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + year + lakePattangansett)
max.lik.fixed.fylP <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],x[4],x[5],0,0,0,x[6],0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,4,5,9)],pen.lik,control=list(maxit=2000))
}

# rate parameter ~  exp( ln(rate0) + form + year13 + lakePattagansett)
max.lik.fixed.fy2013lP <- function(sp.data,w.pars){
	#w.pars c(w.rate,w.shape,beta.form,beta.2013,beta.2014,beta.temp, beta.B,beta.D,beta.P,beta.Q)
	pen.lik <- function(x){
		all.pens <- rep(0,length(x))
		if (x[1] <= 0){
			all.pens[1] <- (x[1])^2
			x[1] <- 0.00001
		}		
		## Think about how to worry if rate gets too big?
		if (x[2] <= 0){
			all.pens[2] <- (x[2])^2
			x[2] <- 0.00001
		}
		unpen.lik <- log.lik.fixed(sp.data,x[1],x[2],x[3],x[4],0,0,0,0,x[5],0,0)
		-1*unpen.lik + sum(all.pens)
	}
	optim(w.pars[c(1,2,3,4,9)],pen.lik,control=list(maxit=2000))
}

