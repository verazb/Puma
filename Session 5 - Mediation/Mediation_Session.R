### AEA2020 ###
##  Mediation ##

# remove old objects from R working space
rm(list=ls())

####################################################################
#   LOAD PACKAGES  
####################################################################

# Load Packages
# Define packages that you need
packages_vector <- c("tidyverse", "haven", "dplyr", "tidyr",  "sandwich", " arsenal",
     				  "xtable", "data.table", "stargazer",  "causalweight")
#install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) 

# Set working directory
work_dir <- "Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/Mediation/Job Corps/"
setwd(work_dir)

# load Job Corps data
data = read_dta("causalmech_clean.dta")

# Store each variable in own R object
attach(data)

# treat                 1=in program group; 0=in control group
# health30              1=exc health 2=good 3=fair 4=poor, 30 months after assignment
# work2year2q           employed in 1st half of second year after assignment
# female                1 if female; 0 if male
# age_cat               age at application in years 16-24
# hhsize                household size at assignment
# hhsizemis             household size at assignment missing
# eduhigh               higher education
# edumis                education missing
# welf1                 once on welfare while growing up
# welf2                 twice on welfare while growing u
# publicass             public assistance in yr before assignment
# fdstamp               received foodstamps in yr before assignment
# rwhite                white
# health0mis            general health at assignment missing
# health012             good or very good health at assignment
# pe_prb0               physical/emotional problems at assignment
# pe_prb0mis            missing - physical/emotional problems at assignment
# everalc               ever abused alcohol before assignment
# everilldrugs          ever took illegal drugs before assignment
# hhinc12               low household income at assignment
# hhinc8                high household income at assignment
# schobef               in school 1yr before eligibility
# jobyrbef              job in year before job corps
# multjobyrbef          number of jobs before Job Corps since reference date
# jobeverbef            ever had a job before Job Corps
# trainyrbef            training in year before Job Corps
# everarr               ever arrested before Job Corps
# pe_prb12              1=phys/emot probs at 12 mths 0=no prob
# vocq4                 in vocational training 9-12 months after assignment



####################################################################
#  DEFINE VARIABLES  
####################################################################

# We focus on the subsample of females

y = health30[female==1]   	# outcome: 1=excellent health after 30 months, 4=poor health
d = treat[female==1]  		# treatment: randomization into Job Corps
m = work2year2q[female==1] 	# meadiator:  employment in the first half of the second year after randomization 

####################################################################
#  MODEL SPECIFICATION  
####################################################################

	# Around this point, might want to... 
	# - Drop observations with no common support, and check how this changes the sample by running descriptive statistics once more
	# - Generate any additional covariates 
	#   (e.g. categories or polynomials for more flexible functional form)

	# Define pre-treatment confounders
	x = cbind( schobef, trainyrbef, jobeverbef, jobyrbef, health012, health0mis,pe_prb0, pe_prb0mis, everalc,
	           everilldrugs, age_cat, edumis, eduhigh, rwhite, everarr, chibef, hhsize, hhsizemis, hhinc12, hhinc8, fdstamp, 
	           welf1, welf2, publicass)[female==1,] 
	x_names = colnames(x)

	# discuss what might be missing (children?)

####################################################################
#  DESCRIPTIVE STATISTICS   
####################################################################
	mydata <-data.frame(y,d,m,x)

	# Summarize outcome and mediator by treatment - suggestive of total effects
	summary(tableby(d~y+m, data = mydata, control=tableby.control(numeric.stats=c("meansd")), total = FALSE))

	# Summarize X by D

	# Define a function estimating the differences in variables across D
	balance_check.model <- function(x){
	    
	    # Conditional means
	    mean_d0 <- mean(x[d==0])
	    mean_d1 <- mean(x[d==1])
	    
	    # Difference in means
	    diff_d <- lm(x ~ d)
	    cov <- vcovHC(diff_d, type = "HC")
	    robust.se <- sqrt(diag(cov))
	    
	    list(mean_d0 = mean_d0, mean_d1 = mean_d1,
	        coeff = diff_d$coefficients[2], 
	        robust.se = robust.se[2], 
	        pval = 2*pnorm(-abs(diff_d$coefficients[2]/robust.se[2])) )             
	}

	diff_output <- apply(x, 2, balance_check.model)

	# convert list to table
	diff_output<-rbindlist(diff_output)
	rownames(diff_output)<- x_names
	colnames(diff_output)<- c("E(X|D=0)", "E(X|D=1)", 
	                          "Difference", "s.e.", 
	                          "p-value")

	# plot table
	print("Difference in means by treatment status")
	xtable(diff_output, digits=3)

	# Summarize X by M 
	# Define a function estimating the differences in variables across M
	balance_check.model.m <- function(x){
	    
	    # Conditional means
	    mean_m0 <- mean(x[m==0])
	    mean_m1 <- mean(x[m==1])
	    
	    # Difference in means
	    diff_m <- lm(x ~ m)
	    cov <- vcovHC(diff_m, type = "HC")
	    robust.se <- sqrt(diag(cov))
	    
	    list(mean_m0 = mean_m0, mean_m1 = mean_m1,
	        coeff = diff_m$coefficients[2], 
	        robust.se = robust.se[2], 
	        pval = 2*pnorm(-abs(diff_m$coefficients[2]/robust.se[2])) )             
	}

	diff_output_m <- apply(x, 2, balance_check.model.m)

	# convert list to table
	diff_output_m<-rbindlist(diff_output_m)
	rownames(diff_output_m)<- x_names
	colnames(diff_output_m)<- c("E(X|M=0)", "E(X|M=1)", 
	                          "Difference", "s.e.", 
	                          "p-value")

	# plot table
	print("Difference in means by mediator status")
	xtable(diff_output_m, digits=3)

####################################################################
#  ESTIMATION OF THE DIRECT AND INDIRECT EFFECTS WITH DIFFERENT
#  ESTIMATORS
####################################################################


### 1. Linear Structural Equation Model (LSE) ###

	###############
	#  Estimation #
	###############

	effects.LSE<-function(y,d,m,x){

	  eq1 =lm(y~d)
	  eq2 =lm(m~cbind(d,x))
	  eq3 =lm(y~cbind(d,m,x))

	  # Total Effect (ATE)
	  te = eq1$coefficients[2]
	  # Direct effects
	  de = eq3$coefficients[2]
	  # Indirect effect(s)
	  ie.1 = eq1$coefficients[2]-eq3$coefficients[2]
	  ie.2 = eq2$coefficients[2]*eq3$coefficients[3]
	  # choose one version

	 list(te=te, de.treat=de, de.control=de, 
	      ie.treat=ie.1, ie.control=ie.1)
	}

	est.LSE<-effects.LSE(y=y,d=d,m=m,x=x)

	print("Output of Function est.LSE")
	est.LSE


	##############
	#  Inference #
	##############

	###  BOOTSTRAP ###
	# Define a function for the bootstrap
	bootstrap.mediation<-function(y,d,m,x,estimator,boot=99){

	obs<-length(y) # vector storing the number of observations
	mat=c() # empty matrix for storing estimators
	temp=c() # empty vector - length counting bootstrap replications
	    
	# The bootstrap loop starts here:
	# as long as temp < number of bootstrap replications that are planned
	while(length(temp)<boot){ 
	    
	# draw a sample    
	sboot<-sample(x=1:obs, # observations that are drawn from y (id)
	              size=obs, # number of obs as in original data
	              replace=TRUE) # with replacement
	    
	# define y, d, m and x in this sample    
	yb=y[sboot] 
	db<-d[sboot]
	mb=m[sboot]
	if (length(x)==length(y)) xb<-x[sboot] 
	if (length(x)!=length(y)) xb<-x[sboot,] # in case of x being a matrix (length = row x columns)

	# estimate the effects using the bootstrap sample 
	# pre-defined estimator!
	est<-c(estimator(yb,db,mb,xb)) 
	    
	# collect the estimated (non-zero) coefficients as additional rows in the matrix
	# one column per effect
	if (sum(is.na(est))==0) mat<-rbind(mat, est) 

	# increase length by 1
	temp<-c(temp,1)
	}

	# store standard deviationa of the estimated effects
	list(sd.te=sd(as.numeric(mat[,1])), 
	     sd.de.treat=sd(as.numeric(mat[,2])), 
	     sd.de.control=sd(as.numeric(mat[,3])),  
	     sd.ie.treat=sd(as.numeric(mat[,4])), 
	     sd.ie.control=sd(as.numeric(mat[,5])) )
	}

	inf.LSE<-bootstrap.mediation(y=y,d=d,m=m,x=x,estimator=effects.LSE,boot=99)
	inf.LSE

	############
	#  Results #
	############

	results.LSE<-rbind(cbind(est.LSE$te, est.LSE$de.treat, est.LSE$de.control, est.LSE$ie.treat, est.LSE$ie.control), 
					cbind(inf.LSE$sd.te, inf.LSE$sd.de.treat, inf.LSE$sd.de.control, inf.LSE$sd.ie.treat, inf.LSE$sd.ie.control),
					cbind(2*pnorm(-abs(est.LSE$te/inf.LSE$sd.te)), 2*pnorm(-abs(est.LSE$de.treat/inf.LSE$sd.de.treat)), 2*pnorm(-abs(est.LSE$de.control/inf.LSE$sd.de.control)),  2*pnorm(-abs(est.LSE$ie.treat/inf.LSE$sd.ie.treat)), 2*pnorm(-abs(est.LSE$ie.control/inf.LSE$sd.ie.control)) )) 
	colnames(results.LSE) <- c("ATE", "de.treat", "de.control", "ie.treat", "ie.control")
	rownames(results.LSE) <- c("effect", "se", "p-val")
	xtable(results.LSE, digits=3)


### 2. Flexible parametric estimator based on the controlled direct effect ###

	###############
	#  Estimation #
	###############

effects.CDE<-function(y,d,m,x){

	mydata <-data.frame(y,d,m,x)

	# Estimate the outcomes model in subsamples o f m and d
	m.d1.m1<-lm(y~.-m-d, data=subset(mydata, d==1 & m==1))
	m.d1.m0<-lm(y~.-m-d, data=subset(mydata, d==1 & m==0))
	m.d0.m1<-lm(y~.-m-d, data=subset(mydata, d==0 & m==1))
	m.d0.m0<-lm(y~.-m-d, data=subset(mydata, d==0 & m==0))

	# Predict yhat out of sample
	yhatd1m1=predict(m.d1.m1, newdata = mydata)
	yhatd1m0=predict(m.d1.m0, newdata = mydata)
	yhatd0m1=predict(m.d0.m1, newdata = mydata)
	yhatd0m0=predict(m.d0.m0, newdata = mydata)

	# Construct CDEs
	de.d1m1 = mean((yhatd1m1-yhatd0m1)[d==1 & m==1])
	de.d1m0 = mean((yhatd1m0-yhatd0m0)[d==1 & m==0])
	de.d0m1 = mean((yhatd1m1-yhatd0m1)[d==0 & m==1])
	de.d0m0 = mean((yhatd1m0-yhatd0m0)[d==0 & m==0])

	# Shares of mediated in d=1 and d=0
	prm1d1=mean(m[d==1])
	prm1d0=mean(m[d==0])

	# Total effect
	te=mean(y[d==1])-mean(y[d==0])

	# Direct effect under treatment
	de.treat= de.d1m0*(1-prm1d1) + de.d1m1*prm1d1
	# Direct effect under control
	de.control= de.d0m0*(1-prm1d0) + de.d0m1*prm1d0
	# Indirect effect under treatemnt
	ie.treat= te - de.control
	# Indirect effect under control
	ie.control= te - de.treat

list(te=te, de.treat=de.treat, de.control=de.control, ie.treat=ie.treat, ie.control=ie.control)
}

	est.CDE <- effects.CDE(y=y,d=d,m=m,x=x)
	est.CDE

	##############
	#  Inference #
	##############

	inf.CDE<-bootstrap.mediation(y=y,d=d,m=m,x=x,estimator=effects.CDE,boot=99)
	inf.CDE

	#############
	#  Results  #
	#############

	results.CDE<-rbind(cbind(est.CDE$te, est.CDE$de.treat, est.CDE$de.control, est.CDE$ie.treat, est.CDE$ie.control), 
					cbind(inf.CDE$sd.te, inf.CDE$sd.de.treat, inf.CDE$sd.de.control, inf.CDE$sd.ie.treat, inf.CDE$sd.ie.control),
					cbind(2*pnorm(-abs(est.CDE$te/inf.CDE$sd.te)), 2*pnorm(-abs(est.CDE$de.treat/inf.CDE$sd.de.treat)), 2*pnorm(-abs(est.CDE$de.control/inf.CDE$sd.de.control)),  2*pnorm(-abs(est.CDE$ie.treat/inf.CDE$sd.ie.treat)), 2*pnorm(-abs(est.CDE$ie.control/inf.CDE$sd.ie.control)) )) 
	colnames(results.CDE) <- c("ATE", "de.treat", "de.control", "ie.treat", "ie.control")
	rownames(results.CDE) <- c("effect", "se", "p-val")
	xtable(results.CDE, digits=3)

	#######################################
	#  Combine parametric Results 		 #
	#######################################

	results.par<-rbind(results.LSE, results.CDE)
	colnames(results.par) <- c("ATE", "de.treat", "de.control", "ie.treat", "ie.control")
	rownames(results.par) <- c("effect LSE", "se.LSE", "p-val.LSE", "effect CDE", "se.CDE", "p-val.CDE")

	print("Results based on the parametric estimators")
	xtable(results.par, digits=3)


### 3. Inverse probability weighting ###

	# Estimate pscores and compare CS for D=1 versus D= 0

	# Estimate the  propensity score p(M,X)
	    mydata <-data.frame(y,d,m,x)
	    mydata$pscore1 = predict(glm(d~m+x, family=binomial(probit)), type="response" )

	# Check for common support in propensity score 
	    # add factor variable for displaying
	    mydata$treat_f <- factor(d, levels = c(0,1), label = c("D=0", "D=1")) 
	    # plot separately by D
	    ggplot(mydata, aes(x = pscore1, fill = treat_f)) + 
	           geom_density(alpha=0.4) + scale_fill_grey() + 
	          theme_bw(base_size = 20) +
	          xlim(0, 1)

	# Estimate the  propensity score p(X)
	    mydata$pscore2 = predict(glm(d~x, family=binomial(probit)), type="response" )

	#  Check for common support in propensity score
	    ggplot(mydata, aes(x = pscore2, fill = treat_f)) + 
	           geom_density(alpha=0.4) + scale_fill_grey() + 
	          theme_bw(base_size = 20) +
	          xlim(0, 1)

# Estimate indirect effect (causalweight package)  - cprrects for small sample imbalances by estimating P(D|X)

		IPW<-medweight(y=y,d=d,m=m,x=x, boot=19, trim = 0.05)
		IPW

		# manual IPW for those who are interested

		  # pscore1 = predict(glm(d~m+x, family=binomial(probit)), type="response" )
		  # pscore2 = mean(d) #w/o small sample correction
		  # # pscore2=glm(d~x,family=binomial(probit),data=mydata)$fitted #with small sample correction
		  
		  # y1m1<-sum(y*d/pscore2)/sum(d/pscore2)
		  # y0m0<-sum(y*(1-d)/(1-pscore2))/sum((1-d)/(1-pscore2))
		  # y1m0<-(sum(y*d*(1-pscore1)/((1-pscore2)*pscore1))/sum(d*(1-pscore1)/((1-pscore2)*pscore1)))
		  # y0m1<-(sum(y*(1-d)* pscore1/(pscore2*(1-pscore1)))/sum((1-d)* pscore1/(pscore2*(1-pscore1))))
		  
		  # de.treat=y1m1-y0m1
		  # ie.control=y0m1-y0m0
		  # de.control=y1m0-y0m0
		  # ie.treat=y1m1-y1m0

	##############################
	#  Combine all Results 		 #
	##############################

		results.all<-rbind( results.LSE, results.CDE, 
				rbind(cbind(IPW$results[1,1], IPW$results[1,2], IPW$results[1,3], IPW$results[1,4], IPW$results[1,5]),
		         cbind(IPW$results[2,1], IPW$results[2,2], IPW$results[2,3], IPW$results[2,4], IPW$results[2,5]),
		         cbind(IPW$results[3,1], IPW$results[3,2], IPW$results[3,3], IPW$results[3,4], IPW$results[3,5])) )
		colnames(results.all) <- c("ATE", "de.treat", "de.control", "ie.treat", "ie.control")
		rownames(results.all) <- c("effect LSE", "se.LSE", "p-val.LSE", 
		                   "effect CDE", "se.CDE", "p-val.CDE",
		                   "effect IPW", "se.IPW", "p-val.IPW")

		print("Results based on all estimators")
		xtable(results.all, digits=3)


####################################################################
 # Small Assignment
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

### Effects for males ###

y = health30[female==0]   	# outcome: 1=excellent health after 30 months, 4=poor health
d = treat[female==0]  		# treatment: randomization into Job Corps
m = work2year2q[female==0] 	# meadiator:  employment in the first half of the second year after randomization 
x = cbind( schobef, trainyrbef, jobeverbef, jobyrbef, health012, health0mis,pe_prb0, pe_prb0mis, everalc,
	       everilldrugs, age_cat, edumis, eduhigh, rwhite, everarr, hhsize, hhsizemis, hhinc12, hhinc8, fdstamp, 
	       welf1, welf2, publicass)[female==0,] 
x_names = colnames(x)

# - CDE - #
est.CDE <- effects.CDE(y=y,d=d,m=m,x=x)
inf.CDE <- bootstrap.mediation(y=y,d=d,m=m,x=x,estimator=effects.CDE,boot=99)

results.CDE.m<-rbind(cbind(est.CDE$te, est.CDE$de.treat, est.CDE$de.control, est.CDE$ie.treat, est.CDE$ie.control), 
				cbind(inf.CDE$sd.te, inf.CDE$sd.de.treat, inf.CDE$sd.de.control, inf.CDE$sd.ie.treat, inf.CDE$sd.ie.control),
				cbind(2*pnorm(-abs(est.CDE$te/inf.CDE$sd.te)), 2*pnorm(-abs(est.CDE$de.treat/inf.CDE$sd.de.treat)), 2*pnorm(-abs(est.CDE$de.control/inf.CDE$sd.de.control)),  2*pnorm(-abs(est.CDE$ie.treat/inf.CDE$sd.ie.treat)), 2*pnorm(-abs(est.CDE$ie.control/inf.CDE$sd.ie.control)) )) 
colnames(results.CDE.m) <- c("ATE", "de.treat", "de.control", "ie.treat", "ie.control")
rownames(results.CDE.m) <- c("effect", "se", "p-val")
xtable(results.CDE.m, digits=3)

# - IPW - #
medweight(y=y,d=d,m=m,x=x, boot=99, trim = 0.05)
