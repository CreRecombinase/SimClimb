### Code to run multispecies multisubstrate trial ###
#Make sure to have run SimpleMainCode
#Also need to have downloaded the CSVs from googledocs and "Infantis_Sep4.2" from Jason's email
#This code will allow you to update the alpha and theta CSVs with the alphas and thetas from the
#best runs for the bifidos. You can choose from three types of infantis, each with a unique leakage 
#scenario. In the modifications section you can create different cases of Bin and Sin changing 
#through time as well as change the initial conditions and the other parameters. 
#Three graphs will appear. (1) A four pane window that shows bifido, bacteroides, others and substrates  
#on separate graphs all of these are on regular scale. (2) substrates on regular scale. (3) all 
#bacteria on regular or log-scale depending if logscale=TRUE/FALSE

#################################### Please Specify The Following #################################
NS = 14
NB = 10

##Give working directory for folder holding the CSVs and Olists
CSV.wd = "C:/Users/nknoblau/Google Drive/Current Model Files/"
##Index to decide which Infantis to use. Half Lactose Leaked=1, Lactose Leaked=2, No Lactose Leaked=3
Inf.index = 1
##Set TRUE if want the most recent parameters stored in Olists, FALSE if want parameters from CSVs
importOlist = TRUE
#Plot bacteria on logscale? TRUE/FALSE
logscale = TRUE

#################################################################################################
##List the file names that hold the Olists for the bifs indicated above. Pls put in same order
##Infantis has 3 Olists for each test type
Inf.Olists = c("Infantis_Sep4.3","Infantis2Olist","Infantis3Olist")
Olist.FileName = c(Inf.Olists[Inf.index],"BifidumOlist","LongumOlist","BreveOlist")	

## Volume of the Infant Colon in ml and milk volume consumed in ml
colonVol <- 100
milkConsumed <- 761.5184 

## Substrate weights in g/mol
mol.weight <- c(LNDFH1=1001.378,LNDFH2=1001.381,LNFP1=855.323,
				LNFP23=855.320,LNnT=709.265,LNT=709.266,
				LDFT=636.248,FL2=490.191,FL3=490.189,LNB=383.35,
				Lac=342.3,Fuc=164.16,Gal=180.16,Glc=180.16)

## Declares parameter list for use in the simulation
pars <- list(alpha = matrix(0,NB,NS),
		fOut = 0.1,
		theta = matrix(0,NB,NS),
		sigma = rep(0,NS),
		rho = low.tri(array(0,dim=c(NS,NS,NB))),
		mu = matrix(0,NB,NB),
		delta = rep(0,NB),
		L = 1,
		Sin = matrix(0,NS,7),  #if use a different breastmilk concentration CSV, change this.
		Bin = matrix(0,NB,10),
		Bnum = NB,
		Snum = NS,
		Total = (NB+NS),
		Order = 1:14
		)
		

### Imports Parameters ###
setwd(CSV.wd)
alpha.read <- read.csv("Alpha.csv",header=TRUE)[,-1]
theta.read <- read.csv("Theta.csv",header=TRUE)[,-1]

#Trim alpha and theta CSVs to exclude the infanti we're not interested in
#Also read in all of the rho matrices to fill the rho array
rho.read  <-  low.tri(array(0,dim=c(NS,NS,NB)))
if (Inf.index == 1){
	alpha.read <- alpha.read[-c(2,3,4), ]
	theta.read <- theta.read[-c(2,3,4), ]
	rho.read[,,1] <- as.matrix(read.csv("InfantisTest1_rho.csv",header=TRUE)[1:NS,2:(NS+1)])	
}else if(Inf.index == 2){
	alpha.read <- alpha.read[-c(1,3,4), ]
	theta.read <- theta.read[-c(1,3,4), ]
	rho.read[,,1] <- as.matrix(read.csv("InfantisTest2_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
}else if(Inf.index ==3){
	alpha.read <- alpha.read[-c(1,2,4), ]
	theta.read <- theta.read[-c(1,2,4), ]
	rho.read[,,1] <- as.matrix(read.csv("InfantisTest3_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
}
rho.read[,,2] <- as.matrix(read.csv("Bifidum_rho.csv",header=TRUE)[1:NS,2:(NS+1)])			
rho.read[,,3] <- as.matrix(read.csv("Longum_rho.csv",header=TRUE)[1:NS,2:(NS+1)])			  
rho.read[,,4] <- as.matrix(read.csv("Breve_rho.csv",header=TRUE)[1:NS,2:(NS+1)])			
rho.read[,,5] <- as.matrix(read.csv("Bacteroidales_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
rho.read[,,6] <- as.matrix(read.csv("Clostridiales_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
rho.read[,,7] <- as.matrix(read.csv("Clostridiales_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
rho.read[,,8] <- as.matrix(read.csv("Lactobacilales_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
rho.read[,,9] <- as.matrix(read.csv("Lactobacilales_rho.csv",header=TRUE)[1:NS,2:(NS+1)])
rho.read[,,10] <- as.matrix(read.csv("Enterobacteriales_rho.csv",header=TRUE)[1:NS,2:(NS+1)])

### Converting Rho from mol/1mol to g/1g and normalizing
for (k in 1:NB){
	for (j in 1:NS){
		#Converts from mol/1mol to g/1g
		pars$rho[ ,j,k] <- pars$rho[ ,j,k] * mol.weight/mol.weight[j]
	}
}

#if want Olist parameters, read in the alphas and thetas from the best saved run
if(importOlist == TRUE){
	for (i in 1:3){
		load(Olist.FileName[i])
		
		num.runs = dim(Olist)[2] 
		score = rep(0,num.runs)
		for (j in 1:num.runs){score[j]=Olist[[1,j]]$score}
		broken = max(which(score != 0))
		best = which(score==min(score[1:broken]))[1]
		
		alpha.read[i,] <- Olist[[1,best]]$alpha
		theta.read[i,] <- Olist[[1,best]]$theta
		
		#if infantis then update the leakage constant
		if (i==1){pars$L <- Olist[[1,best]]$L}
	}
#if don't want Olist parameters, set the leakage constant to these predetermined values
}else{
	if (Inf.index == 1){
		pars$L <- 0.1716188
	}else if(Inf.index == 2){
		pars$L <- 0.08617654
	}else if(Inf.index == 3){
		stop("Need to set importOlist must be TRUE!!! <3 Nick and Kim")
	}
}

rownames(pars$Bin) <- c("Infantis", "Bifidum", "Longum", "Breve", "Bacteroides", "Clostridium",
			"Veillonella","Streptococcus","Enterococcus","Enterobacteria")
pars$alpha <- alpha.read*24
pars$theta <- theta.read
pars$rho <- rho.read
rownames(pars$rho) <- names(mol.weight)

########################################## MODIFICATIONS #######################################################
#Substrate Order: (1)LNDFH I, (2)LNDFH II, (3)LNFP I, (4)LNFP II/III, (5)LNnT, (6)LNT, (7)LDFT, (8)2'FL, (9)3'FL
#					(10)LNB, (11)Lac, (12)Fuc, (13)Gal, (14)Glu
#Bacteria Order: (1)Infantis, (2)Bifidum, (3)Longum, (4)Breve, (5)Bacteroides, (6)Clostridium, (7)Veillonella
#					(8)Streptococcus, (9)Enterococcus, (10)Escherichia


#####CHANGE Sin:
	#This section allows you to change the way substrate is entering the system. 
	#Sin is a matrix whose rows are substrates and whose columns are 7 time points (days 3,8,15,22,30,60,90)
	#Each row shows how the concentration of that substrate varies through time. 
	#The matrix is used by the function S_in( ) in the SimpleMainCode which fits a polynomial to the 7 time points
	#for each of the substrates so that substrates can change gradually from one concentration to the next. 
	#In the case that the substrates stay constant, the polynomials are also constant through time.
	
	PlotSin = FALSE			#If TRUE then the polynomials that fit to the Sin's per time will be displayed when S_in( ) is called.
							#To see the graph of the polynomials run this section of code then type in the console: "S_in(1:90)"
							#If running entire simulation set to FALSE to speed up simulation
	
	VarySin = TRUE			#If TRUE and no other changes are made then only HMOs will enter the system each day
							#and the concentrations of the HMOs will vary according to the observations of a recent study. 
							#If FALSE then all substrates will enter the system at a constant rate specified by ConstSinValue.
	ConstSinValue = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1) #Corresponds to the constant Sin for each of the 14 substrates
	
	lacPercent <- 0 ## Percentage of Lac in milk that enters colon
	if(VarySin){
		SinRead <- as.matrix(read.csv(file="Breastmilk_composition_over_time.csv",header=FALSE))
		SinRead[2:11, ] <- SinRead[2:11, ] * milkConsumed / colonVol
		pars$Sin[1:9, ] <- SinRead[2:10, ]
		pars$Sin[11, ] <- SinRead[12, ]*lacPercent ##Adds the percentage of Lac that makes it to the colon
		SDays <- SinRead[1, ]
	}else{
		pars$Sin <- matrix(ConstSinValue, NS, 7)
		SDays <- as.matrix(read.csv(file="Breastmilk_composition_over_time.csv",header=FALSE))[1, ]
	}
	
	#pars$Sin[11,] <- c(0,0,0,4,0,0,0)		#This is just an example of how substrate can be further modified
											#You can change any substrate at any of the 7 days by assigning the rows of Sin this way. 
	#pars$Sin[6,7] <- rep(3,7)
	
	#Create coefficients for Sin polynomials
	polyDegreeS <- 3
	coeffS <- matrix(0,NS,(polyDegreeS+1))
	for (i in 1:NS){
		coeffS[i,] <- Pfit(SDays,pars$Sin[i,],polyDegreeS)
	}
	
###SCRAMBLE Sin 
	#This section scrambles the abundances of the substrates which was set in the CHANGE Sin section
	#For example pars$Order <- c(1,2,3,4,6,5,7,8,9,10,11,12,13,14) will give LnNT the concentrations of LNT
	#and LNT will have the concentrations of LnNT because the 5,6 are swapped. This will test if the relative 
	#abundances of HMO in breastmilk alone can impact the bacterial community since the total amount of HMO doesn't change.
	
	#pars$Order <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)	#This will rearrange the abundances of the substrates, 1:14 is default
	#pars$Sin <- pars$Sin[pars$Order, ]
	
#####CHANGE Bin:
	#This section allows you to vary the bacteria coming in at 10 possible time points.
	#pars$Bin is a matrix whose rows are the bacteria and whose columns are time points.
	#This matrix is used by the function B_in( ) which is very similar to S_in( ) (see above). 
	#The times at which the bacterial concentration coming in can vary is determined by the
	#Bdays vector. By changing the days in this vector, you can have most of the Bin variation
	#occur at the beginning of the simulation, or the end, or you can have the Bin variation 
	#span the whole simulation. NOTE: The simulation runs to day 90. 

	PlotBin = FALSE		#If TRUE then the polynomials that fit to the Bin's per time will be displayed when B_in( ) is called.
						#To see the graph of the polynomials run this section of code then type in the console: "B_in(1:90)"
						#If running entire simulation set to FALSE to speed up simulation

	
	VaryBin = FALSE		#If want constant Bin then set FALSE and enter the constant value for Bin specified by ConstBinValue
	ConstBinValue = c(0,0,0,0,0,0,0,0,0,0) #Corresponds to the constant Bin for each of the 10 bacteria
	if(VaryBin){
		BDays <- 						c(  0, 10, 20, 30, 40, 50, 60, 70, 80, 90)  #Days the Bin can vary
		pars$Bin["Infantis",] <- 		c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)	#Enter the way you'd like to vary Bin for each 
		pars$Bin["Bifidum",] <- 		c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)	#bacteria over the specified 10 time points. 
		pars$Bin["Longum",] <- 			c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)	#units: [M CFU/ ml / day]
		pars$Bin["Breve",] <- 			c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)
		pars$Bin["Bacteroides",] <- 	c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)	#If want any row to be const then may use
		pars$Bin["Clostridium",] <- 	c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)	#rep(1,NB) command to fill row with 1's, or
		pars$Bin["Veillonella",] <- 	c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)	#rep(2,NB) to fill row with 2's etc. 
		pars$Bin["Streptococcus",] <- 	c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)
		pars$Bin["Enterococcus",] <- 	c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)
		pars$Bin["Enterobacteria",] <-	c(.01,.01,.01,.01,.01,.01,.01,.01,.01,.01)
	}else{
		pars$Bin <- matrix(ConstBinValue,NB,10)
		BDays <- seq(0,90,length=10)
	}

#####CHANGE who's there:
	#This section allows you to specify what bacteria are present in the simulation.
	#The "present" vector has the numbers corresponding to the bacteria you want to include.
	#(This change has to occur before the coeffs are calculated for pars$Bin)
	present <- c(1,2,3,4,5,6,7,8,9,10)
	pars$Bin[setdiff(1:NB,present),] <- rep(0,10)
	
	#create coefficients for Bin polynomials 
	polyDegreeB <- 3	#This must be less than the number of unique Bin points for each of the bacteria
	coeffB <- matrix(0,NB,(polyDegreeB+1))
	for (i in 1:NB){
		coeffB[i,] <- Pfit(BDays,pars$Bin[i,],polyDegreeB)
	}

#####TO MODIFY LNB alphas AND thetas FOR INFANTIS AND BREVE (b/c not confident in values)
	#pars$alpha[1,10] <- pars$alpha[1,11]*2				#Give Infantis LNB alpha 2x the lactose alpha
	#pars$alpha[4,10] <- pars$alpha[4,11]*2				#Give Breve LNB alpha 2x the lactose alpha
	pars$theta[1,10] <- (pars$theta[2,10]+pars$theta[3,10])/2		#Give Infantis LNB theta the average of longum and bifidum theta
	pars$theta[4,10] <- (pars$theta[2,10]+pars$theta[3,10])/2		#Give Breve LNB theta the average of longum and bifidum theta
	pars$alpha[1,10] <- (pars$alpha[2,10]+pars$alpha[3,10])/2		#Give Infantis LNB alpha the average of longum and bifidum alpha
	pars$alpha[4,10] <- (pars$alpha[2,10]+pars$alpha[3,10])/2		#Give Breve LNB alpha the average of longum and bifidum alpha

#####INITIAL CONDITIONS & OTHER PARAMETERS
pars$fOut <- 0.1
pars$sigma <- 0
#pars$mu[3,1] <- -.0001	
TotalBactStart <- 1 #Total MCFU/ml in initial inoculation
RelativeBactStart <- c(1,1,1,1,1,1,1,1,1,1) #Relative abundances of Bacteria in initial inoculation
BactStart <- TotalBactStart*RelativeBactStart
BactStart[setdiff(1:NB,present)] <- 0
SubStart <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
yini <- c(BactStart,SubStart)
times <- seq(1,2000,by=1)

# pars$L <- 0
### Runs Simulation ###

graphics.off()
Rprof(filename="Profiling.txt")
out <- sim()
Rprof(filename=NULL)
summaryRprof("Profiling.txt")
######### Saving your run ###########
### To save, uncomment this code, enter the file name in quotes for saveName
### and run the code.  Remember to also create a txt file of the same name that briefly
### describes the run and contains a copy of the RunTrial code that you used.
### Naming convention:  Begin each file name with the first 3 letters of your first name
### For example:  jasinfantisrun.txt

# saveName <- "savetest"
# write.csv(out,file=saveName)


