### Gut Model ###

library(deSolve) #loads DE solver package

## mainFunc contains the differential equations for use by deSolve.
mainFunc <- function(Time, State, Pars) {
	with(as.list(c(State)), { #converts the State variable vector to a list
	
	#Splits state vector in to bacteria and substrate components for easier reading of dB and dS
	B <- State[1:Pars$Bnum]
	S <- State[(Pars$Bnum+1):Pars$Total]
	
	##Adjusts Infantis Monomer production by leakage factor
	leakers <- c("Lac","Fuc","Gal","Glc")
	Pars$rho[leakers, ,1] = Pars$rho[leakers, ,1]*Pars$L
	##Calculates vector P of total substrate produced for each S
	P <- rowSums(apply(aperm(as.vector(S*t(Pars$alpha*B))*aperm(Pars$rho,c(2,3,1)),c(3,1,2)),c(1,2),sum))
	
	#Calculates the change in B and S at a given time
	#Vector notation is more efficient, but also harder to read
	dB <- B_in(Time) + colSums(t(Pars$theta*Pars$alpha*B)*S) - Pars$fOut*B #+ rowSums(t(Pars$mu*B)*B)
	dS <- S_in(Time,Pars$Order) + P - rowSums(S*t(Pars$alpha*B)) -Pars$sigma*S - Pars$fOut*S 

	return(list(c(dB, dS)))
	})
}

## Function for converting each matrix in a XxYxZ array to lower tri
low.tri <- function(A){
	M <- length(A[1,1, ])
	for (i in 1:M){
		a <- A[ , ,i]
		a[upper.tri(a,diag=TRUE)]=0
		A[ , ,i] <- a
	}
	return(A)
}

## Runs the simulation, displays numerical and graphical results
sim <- function(Yini = yini, Parms = pars, Times = times){
	if (length(which(Yini[1:Parms$Bnum] == 0)) > 0){ 
		out<-ode(func = mainFunc, y = Yini, parms = Parms, times = Times, hmax=0.02)
	}else{
		out<-ode(func = mainFunc, y = Yini, parms = Parms, times = Times)
	}
	col.vect <- c("black",colors()[c(254,26,34,51,62,67,500,93,116,84,508,655,656)])
	windows()
	if(logscale == 1){
		par(ps=16)
		matplot(out[,"time"], out[,2:(NB+1)], log ="y", type = "l", xlab = "Time (Days)", ylab = "Bacterial Concentration [Mega CFU/ml]",
			main = "Bacterial Concentration", lwd = 3, col = col.vect[1:NB], lty=1)
		legend("bottomleft", c(paste("Infantis",Inf.index), "Bifidum", "Longum", "Breve", "Bacteroides", "Clostridium",
		   "Veillonella","Streptococcus","Enterococcus","Escherichia"), col = col.vect[1:NB], cex = 0.6, lwd=3, lty=1)
	}else{
		par(ps=16)
		matplot(out[,"time"], out[,2:(NB+1)], type = "l", xlab = "Time (Days)", ylab = "Bacterial Concentration [Mega CFU/ml]",
			main = "Bacterial Concentration", lwd = 3, col = col.vect[1:NB], lty = 1 )
		legend("topright", c(paste("Infantis",Inf.index), "Bifidum", "Longum", "Breve", "Bacteroides", "Clostridium",
		   "Veillonella","Streptococcus","Enterococcus","Escherichia"), col = col.vect[1:NB], lty = 1, cex = 0.6, lwd=3)
	}

	windows()
	par(ps=16)
	matplot(out[10:length(Times),"time"], out[10:length(Times),(NB+2):(NB+NS+1)], type = "l", xlab = "Time (Days)", ylab = "Substrate Concentration [g/l]",
		main = "Substrate Depletion", lwd = 3, col = col.vect[1:NS], lty=1)
	legend("topright", names(mol.weight), col = col.vect[1:NS], cex = 0.6, lwd=3, lty=1)

	windows() #Creates a new window to plot results of this optim run against the data
	par(mfrow=c(2,2), oma = c(0,0,2,0)) #Divides window for bacteria and substrate plots
	matplot(out[,"time"],out[,2:5], type = "l", xlab = "Time (Days)", ylab = "Bacterial Conc. [Mega CFU/ml]",
			main = "Bifidobacterial Concentration", lwd = 3, col = col.vect[1:4], lty = 1 )
	legend("topright", c(paste("Infantis",Inf.index), "Bifidum", "Longum", "Breve"), col = col.vect[1:4], lty = 1, cex = 0.6, lwd=3)
	matplot(out[,"time"],out[,6], type = "l", xlab = "Time (Days)", ylab = "Bacterial Conc. [Mega CFU/ml]",
			main = "Bacterodiales Concentration", lwd = 3, col = col.vect[5:5], lty = 1 )
	legend("topright", c("Bacteroides"), col = col.vect[5:5], lty = 1, cex = 0.6, lwd=3)
	matplot(out[,"time"],out[,7:(NB+1)], type = "l", xlab = "Time (Days)", ylab = "Bacterial Conc. [Mega CFU/ml]",
			main = "Other Concentration", lwd = 3, col = col.vect[6:NB], lty = 1 )
	legend("topright", c("Clostridium","Veillonella","Streptococcus","Enterococcus","Escherichia"), 
			col = col.vect[6:NB], lty = 1, cex = 0.6, lwd=3)
	matplot(out[10:length(Times),"time"],out[10:length(Times),(NB+2):(NB+NS+1)], type = "l", xlab = "Time (Days)", ylab = "Substrate Conc. [g/l]",
			main = "Substrate Depletion", lwd = 3, col = col.vect[1:NS], lty=1)
	legend("topright", c("LNDFH I", "LNDFH II", "LNFP I", "LNFP II/III", "LNnT", "LNT", "LDFT", "2' FL", "3' FL", 
			"LNB", "Lac", "Fuc","Gal","Glu"), col = col.vect[1:NS], cex = 0.6, lwd=3, lty=1)
	return(out)
}

##Asakuma plotting function.  Completely general and for use with Olists that contain their own Akasuma data
Aplot <- function(OL,I){
	Gdata <- OL$Gdata
	Gtimes <- Gdata[ ,1]
	Hdata <- OL$Hdata
	Htimes <- Hdata[ ,1]
	both.times <- union(Gtimes,Htimes)
	S.num <- length(OL$Hdata[1 ,]) - 1
	
	windows() #Creates a new window to plot results of this optim run against the data
	par(mfrow=c(2,2), oma = c(0,0,2,0)) #Divides window for bacteria and substrate plots
	
	matplot(Gtimes,Gdata[ ,2],type="o",col=1, pch=1,lty=1)
	points(Gtimes,OL$oderesults[both.times %in% Gtimes,2],type="o",col=1, pch=2,lty=2)

	matplot(Htimes,Hdata[ ,2:5],type="o",col=1:4, pch=1,lty=1)
	matpoints(Htimes,OL$oderesults[both.times %in% Htimes,3:6],type="o",col=1:4, pch=2,lty=2)
	legend("topright",colnames(Hdata)[2:5],col=1:4,lty=rep(1,4),cex=0.6)
	
	matplot(Htimes,Hdata[ ,6:10],type="o",col=5:9, pch=1,lty=1)
	matpoints(Htimes,OL$oderesults[both.times %in% Htimes,7:11],type="o",col=5:9, pch=2,lty=2)
	legend("topright",colnames(Hdata)[6:10],col=5:9,lty=rep(1,4),cex=0.6)
	
	matplot(Htimes,Hdata[ ,11:(S.num+1)],type="o",col=10:S.num, pch=1,lty=1)
	matpoints(Htimes,OL$oderesults[both.times %in% Htimes,12:(S.num+2)],type="o",col=10:S.num, pch=2,lty=2)
	legend("topright",colnames(Hdata)[11:(S.num+1)],col=10:S.num,lty=rep(1,(S.num-10)),cex=0.6)
	
	title(paste("Run:",I,"Score:",OL$score,"Converg:",OL$convergence,"Wt:",OL$weight), outer=TRUE)
}

##Old Asakuma plotting function for use with Olists that don't contain the growth and substrate data
##You must manually run the code to generate the appropriate AGtimes,AHtimes,AGdata2, and AHdata2 first
Aplot.old <- function(OL,I){
	windows() #Creates a new window to plot results of this optim run against the data
	par(mfrow=c(2,2), oma = c(0,0,2,0)) #Divides window for bacteria and substrate plots
	
	matplot(AGtimes,AGdata2[ ,2],type="o",col=1, pch=1,lty=1)
	points(AGtimes,OL$oderesults[times %in% AGtimes,2],type="o",col=1, pch=2,lty=2)

	matplot(AHtimes,AHdata2[ ,2:5],type="o",col=1:4, pch=1,lty=1)
	matpoints(AHtimes,OL$oderesults[times %in% AHtimes,3:6],type="o",col=1:4, pch=2,lty=2)
	legend("topright",colnames(AHdata2)[2:5],col=1:4,lty=rep(1,4),cex=0.6)
	
	matplot(AHtimes,AHdata2[ ,6:10],type="o",col=5:9, pch=1,lty=1)
	matpoints(AHtimes,OL$oderesults[times %in% AHtimes,7:11],type="o",col=5:9, pch=2,lty=2)
	legend("topright",colnames(AHdata2)[6:10],col=5:9,lty=rep(1,4),cex=0.6)
	
	matplot(AHtimes,AHdata2[ ,11:(NS+1)],type="o",col=10:NS, pch=1,lty=1)
	matpoints(AHtimes,OL$oderesults[times %in% AHtimes,12:(NS+2)],type="o",col=10:NS, pch=2,lty=2)
	legend("topright",colnames(AHdata2)[11:(NS+1)],col=10:NS,lty=rep(1,(NS-10)),cex=0.6)
	
	title(paste("Run:",I,"Score:",OL$score,"Converg:",OL$convergence,"Wt:",OL$weight), outer=TRUE)
}

## Function plots the best trial for each weight in an Olist
plot.best <- function(OL){
	scores <- getScores(OL)
	for (i in 1:length(scores[ ,1])){
		best <- which(scores[i, ] == min(scores[i,which(scores[i, ] != 0)]))[1]
		Aplot(OL[[i,best]],best)
	}
}

## Fitness Functions  
fit.withdeath <- function(Est, Fixed){
	Fixed$alpha[alpha.index] = t(exp(Est[1:num.A]))
	Fixed$theta[theta.index] = t(exp(Est[(num.A+1):num.I]))
	Fixed$delta = exp(Est[(num.I+1)])
	colnames(Fixed$alpha)=NULL
	colnames(Fixed$theta)=NULL
	names(Fixed$delta)=NULL
	out <- ode(func = mainFunc, y = yini, parms = Fixed, times = times)
	value <- sum((out[times %in% AGtimes,2:(NB+1)]-AGdata2[,2:(NB+1)])^2)
	value <- value + sum((fit.weight*t(weight.dist*t(out[times %in% AHtimes,(NB+2):(pars$Total+1)]-AHdata2[,2:(NS+1)])))^2)
	return(value)
} 
fit.nodeath <- function(Est, Fixed){
	thetaTemp <- rep(0,Fixed$Snum)
	thetaTemp[theta.index] <- t(exp(Est[(num.A+1):num.I]))
	thetaTemp[c(2,4,9,14)] <- thetaTemp[c(1,3,8,13)]
	if (B.index == 1){
		Fixed$L = exp(Est[(num.I+1)])
		thetaTemp[12] <- thetaTemp[13]
	}
	Fixed$alpha[alpha.index] = t(exp(Est[1:num.A]))
	Fixed$theta = thetaTemp
	colnames(Fixed$alpha)=NULL
	colnames(Fixed$theta)=NULL
	out <- ode(func = mainFunc, y = yini, parms = Fixed, times = times)
	value <- sum((out[times %in% AGtimes,2:(NB+1)]-AGdata2[,2:(NB+1)])^2)
	value <- value + sum((fit.weight*t(weight.dist*t(out[times %in% AHtimes,(NB+2):(pars$Total+1)]-AHdata2[,2:(NS+1)])))^2)
	return(value)
}
fit.singlesub <- function(Est, Fixed){
	Fixed$alpha[1,alpha.index] = exp(Est[1])
	Fixed$theta[1,theta.index] = exp(Est[2])
	colnames(Fixed$alpha)=NULL
	colnames(Fixed$theta)=NULL
	out <- ode(func = mainFunc, y = yini, parms = Fixed, times = times)
	value <- sum((out[times %in% AGtimes,2:(NB+1)]-AGdata2[,2:(NB+1)])^2)
	value <- value + sum((fit.weight*(out[times %in% AHtimes,(theta.index+2)]-AHdata2[,(theta.index+1)]))^2)
	return(value)
}

## Function for creating matrix of fitness scores
getScores <- function(OL){
	row.num <- length(OL[ ,1])
	col.num <- length(OL[1, ])
	M <- matrix(0,row.num,col.num)
	for(i in 1:row.num){
		for(j in 1:col.num){
			M[i,j] <- OL[[i,j]]$score
		}
	}
	return(M)
}

import.data <- function(i){
	if(i == 1){
		rho.read <- read.csv("InfantisTestA_rho.csv",header=TRUE)
		AG.read <- read.csv("Asakuma_Infantis_growth.csv",header=TRUE)
		AH.read <- read.csv("Infantis_HMO.csv")
		AG.read <- AG.read[-1:-2, ]
		AH.read <- AH.read[-1:-2, ]
		AG.read[ ,1] <- AG.read[ ,1] - AG.read[1,1]
		AH.read[ ,1] <- AH.read[ ,1] - AH.read[1,1]	
	}else if(i == 2){
		rho.read <- read.csv("Bifidum_rho.csv",header=TRUE)
		AG.read <- read.csv("Asakuma_Bifidum_growth.csv",header=TRUE)
		AH.read <- read.csv("Bifidum_HMO.csv")
	}else if(i == 3){
		rho.read <- read.csv("Longum_rho.csv",header=TRUE)
		AG.read <- read.csv("Asakuma_Longum_growth.csv",header=TRUE)
		AH.read <- read.csv("Longum_HMO.csv")		
	}else if(i == 4){
		rho.read <- read.csv("Breve_rho.csv",header=TRUE)
		AG.read <- read.csv("Asakuma_Breve_growth.csv",header=TRUE)
		AH.read <- read.csv("Breve_HMO.csv")		
	}else{
		print("Invalid Bacterial Index.  Select 1,2,3, or 4 for Infantis,Bifidum,Longum, or Breve")	
	}
	return(list(rho.read = rho.read, AG.read = AG.read, AH.read = AH.read))
}

getEst <- function(P){
	## Create vector of names for Est
	if (B.index == 1){
		est.names <- rep("0",(num.I+1))	
		est.names[(num.I+1)] <- "L"			
	}else{
		est.names <- rep("0",num.I)
	}
	for (j in 1:num.A){
		est.names[j] <- paste("alpha",alpha.index[j],sep="")
	}
	for (j in (num.A+1):num.I){
		est.names[j] <- paste("theta",theta.index[(j-num.A)],sep="")
	}
	
	odeCheck <- -1
	while(odeCheck != 2){
		# Fills Est vector
		if (B.index == 1){
			Est <- rep(0,(num.I+1))
			Est[(num.I+1)] <- runif(1,log(L.est - L.est*range[(num.I+1)]),log(L.est + L.est*range[(num.I+1)]))			
		}else{
			Est <- rep(0,num.I)
		}
		for(j in 1:num.A){
			Est[j] <- runif(1,log(alpha.est[j] - alpha.est[j]*range[j]),log(alpha.est[j] + alpha.est[j]*range[j]))		
		}
		for (j in (num.A+1):num.I){
			Est[j] <- runif(1,log(theta.est[j-num.A] - theta.est[j-num.A]*range[j]),log(theta.est[j-num.A] + theta.est[j-num.A]*range[j]))
		}
		
		# Fills alpha and theta vectors for use in ode solver
		alphaTemp <- rep(0,pars$Snum)		
		thetaTemp <- rep(0,pars$Snum)
		alphaTemp[alpha.index] <- exp(Est[1:num.A])
		thetaTemp[theta.index] <- exp(Est[(num.A+1):num.I])
		thetaTemp[c(2,4,9,14)] <- thetaTemp[c(1,3,8,13)]
		if (B.index == 1){
			thetaTemp[12] <- thetaTemp[13]
		}
		
		# Calls ode solver using Est values
		P$alpha = t(alphaTemp) 
		P$theta = t(thetaTemp) 
		if (B.index == 1){
			P$L <- exp(Est[(num.I+1)])
		}
		# Resets odeCheck.  A 2 value indicates successfull integration by deSolve
		odeCheck <- diagnostics.deSolve(sim(Parms = P))$istate[1]	
	}	
	names(Est)=est.names
	return(Est)
}

#Calculates coefficients based on values from Sin matrix
Pfit <- function(x,y,Psize){
	model <- lm(y ~ poly(x, Psize, raw=TRUE))
	return(model$coefficients)
}

#Call Sin to calculate substrates at time evalTime
S_in <- function(evalTime,Order){
	SubPerTime <- matrix(0,NS,length(evalTime))
	for (i in 1:length(evalTime)){
		for (j in 0:polyDegreeS){
			SubPerTime[,i] <- SubPerTime[,i]+coeffS[,(j+1)]*evalTime[i]^j
		}
		SubPerTime[which(SubPerTime[,i] < 0),i] = 0
	}
	
	SubPerTime <- SubPerTime[Order,]
	if(PlotSin){
		col.vect <- c("black",colors()[c(254,26,34,51,62,67,500,93,116,84,508,655,656)])
		matplot(evalTime,t(SubPerTime),type="l",pch=1,col = col.vect[1:NS],lty=1,lwd=3,
			xlab="SDays",ylab="HMO Concentration [g/l]",main="Changes in HMO in time")
		legend("topright",c("LNDFH I","LNDFH II","LNFP","LNFP II/III","LNnT","LNT",
			"LDFT","2' FL","3' FL","LNB","Lac","Fuc","Gal","Glc"),
			col = col.vect[1:NS],lty=1,cex=0.6)
		for(i in 1:NS){
			points(SDays,pars$Sin[i,],col=col.vect[i])
		}
	}
	return(SubPerTime)
}

#Call Sin to calculate substrates at time evalTime
B_in <- function(evalTime){
	BacPerTime <- matrix(0,NB,length(evalTime))
	for (i in 1:length(evalTime)){
		for (j in 0:polyDegreeB){
			BacPerTime[,i] <- BacPerTime[,i]+coeffB[,(j+1)]*evalTime[i]^j
		}
		BacPerTime[which(BacPerTime[,i] < 0),i] = 0
	}	
	if(PlotBin){
		col.vect <- c("black",colors()[c(254,26,34,51,62,67,500,93,116,84,508,655,656)])
		matplot(evalTime,t(BacPerTime),type="l",pch=1,col = col.vect[1:NB],lty=1,lwd=3,
			xlab="BDays",ylab="Bacterial Concentration [M CFU/ml]",main="Changes in Bin through Time")
		legend("topright",c("Infantis", "Bifidum", "Longum", "Breve", "Bacteroides", "Clostridium",
				"Veillonella","Streptococcus","Enterococcus","Enterobacteria"),
			col = col.vect[1:NB],lty=1,cex=0.6)
		for(i in 1:NB){
			points(BDays,pars$Bin[i,],col=col.vect[i])
		}
	}
	return(BacPerTime)
}