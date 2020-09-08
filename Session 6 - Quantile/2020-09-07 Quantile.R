### AEA2020 ###

# Regression discontinuity design

# CIT refers to Cattaneo Idrobo Titunik, 2019 


####################################################################
#   PHASE 0: LOAD PACKAGES AND READ IN DATA  
####################################################################

rm(list=ls())


packages_vector <- c("haven", "dplyr", "tidyr", "sandwich", "expss",
    				"fBasics", "xtable", "data.table", "stargazer", "mfx", 
    				"jtools", "ggplot2", "counterfactual")
#install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) 


# Quantile-specific packages 
packaged_vector_rdd <- c("quantreg")
# install.packages(packages_vector)
lapply(packaged_vector_rdd, require, character.only = TRUE) 


# List loaded packages 
(.packages())


# Working Directory
work_dir <- "Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/Quantile/Melly 2005/"
setwd(work_dir)


# Read in the data
data <- as.data.frame(read_dta("data.dta")) 



####################################################################
#   PHASE 1: INSPECT YOUR DATA AND DEFINE KEY VARIABLES
####################################################################

# Inspect data
head(data)


# Vector with all variable names
varnames <- colnames(data)


# Check lnwage positive
sum(data$lnwage>0)


# store each variable in own R object
attach(data)


####################################################################
#   PHASE 2: DESCRIPTIVE STATISTICS AND DISTRIBUTION
####################################################################


# Number of treated and control
cro(public)


# Unconditional wage distribution



# Summary statistics of wages by sector
# ----------------------------------------------------

mean_d0 <- mean(lnwage[public==1])
mean_d1 <- mean(lnwage[public==0])

# Difference in means
diff_d <- lm(lnwage ~ public)
cov <- vcovHC(diff_d, type = "HC")
robust.se <- sqrt(diag(cov))




# Compute quantiles
# ----------------------------------------------------

quantiles_all <- quantile(lnwage, probs = seq(0, 1, 0.2), na.rm = FALSE,
         names = TRUE, type = 7)

quantiles_public <- quantile(lnwage[public==1], probs = seq(0, 1, 0.2), na.rm = FALSE,
         names = TRUE, type = 7)

quantiles_private <- quantile(lnwage[public==0], probs = seq(0, 1, 0.2), na.rm = FALSE,
         names = TRUE, type = 7)


quantiles_all
quantiles_public
quantiles_private


# Summary statistics of covariates by sector
# ----------------------------------------------------

# Selected covariates
covariates <- dplyr::select(data, alter, frau) 
covariates_names  <- colnames(covariates)


# Define a function estimating the differences in variables across sectors 
balance_check.model <- function(x){
    
    # Conditional means
    mean_d0 <- mean(x[data$public==1])
    mean_d1 <- mean(x[data$public==0])
    
    # Difference in means
    diff_d <- lm(x ~ data$public)
    cov <- vcovHC(diff_d, type = "HC")
    robust.se <- sqrt(diag(cov))
    
    list(mean_d0 = mean_d0, mean_d1 = mean_d1,
        coeff = diff_d$coefficients[2], 
        robust.se = robust.se[2], 
        pval = 2*pnorm(-abs(diff_d$coefficients[2]/robust.se[2])) )             
}

diff_output <- apply(covariates, 2, balance_check.model)

# convert list to table
diff_output<-rbindlist(diff_output)
rownames(diff_output)<- covariates_names
colnames(diff_output)<- c("E(X|Public)", "E(X|Private)", "Difference", "s.e.", "p-value")

# plot table
print("Difference in means by treatment status")
xtable(diff_output, digits=3)







# Basic quantile plots
# ----------------------------------------------------

# https://mgimond.github.io/ES218/Week05a.html


# # All - probably don't need that 

# 	lnwage_sorted <- sort(lnwage) 
# 	length <- length(lnwage_sorted)

# 	i <- 1:length
# 	perc <- (i - 1)/(length - 1) # compute percentiles
# 	data_plot <- data.frame(lnwage_sorted, perc)

# 	qplot <- ggplot(data_plot,                
# 	       aes(x = perc, y = lnwage_sorted)) +
# 	 	   geom_line() +
# 	 	   ylab("Sample quantile") +
# 	 	   xlab("Sample fraction") +
# 	 	   theme_bw(base_size = 20) 
# 	qplot


# Density plots
# ----------------------------------------------------

# Distribution

	# Factor public variable for plots
	data$public_factor <- factor(data$public, 
	                        levels = c(0,1), label = c("Private", "Public")) 


	density_plot <- ggplot(data, aes(x=lnwage, group=public_factor, color=public_factor, fill=public_factor)) +
						geom_density(alpha=0.4) +
						theme_bw(base_size = 20)

	density_plot



# Cumulative Density Functions
	cdfx_pub <- c(1:length(lnwage[public==1]))/length(lnwage[public==1])
	cdfx_priv <- c(1:length(lnwage[public==0]))/length(lnwage[public==0])
	plot(c(6,14),range(c(0,1)), xlim =c(6, 14), type="n", xlab = "log wage", ylab="Probability")
	lines(sort(lnwage[public==1]), cdfx_pub)
	lines(sort(lnwage[public==0]), cdfx_priv, lwd = 2, col = "grey70")
	legend(6, .2, c("Public", "Private"), lwd=c(1,2),bty="n", col=c(1,"grey70"))




# Quantile Functions

	# Public
	lnwage_sorted <- sort(lnwage[public==1]) 
	length <- length(lnwage_sorted)
	i <- 1:length
	perc <- (i - 1)/(length - 1) # Explain where this comes from 
	data_plot_public <- data.frame(lnwage_sorted, perc)
	data_plot_public$public <- 1

	# Private
	lnwage_sorted <- sort(lnwage[public==0]) 
	length <- length(lnwage_sorted)
	i <- 1:length
	perc <- (i - 1)/(length - 1) 
	data_plot_private <- data.frame(lnwage_sorted, perc)
	data_plot_private$public <- 0


	data_plot <- rbind(data_plot_public, data_plot_private)
	data_plot$public_factor <- factor(data_plot$public, 
	                        levels = c(0,1), label = c("Private", "Public")) 

	# Plot
	qplot_public_private <- ggplot(data_plot, aes(x = perc, y = lnwage_sorted, group=public_factor, color=public_factor)) +
			geom_line() + 
	 	   	ylab("Sample quantile of ln(wage)") +
	 	   	xlab("Sample fraction") +
	 	   	theme_bw(base_size = 20) 
	qplot_public_private




####################################################################
#   PHASE 3: Quantile effect estimation and robustness checks
####################################################################


# Unconditional OLS
# ----------------------------------------------------

ols1 <- lm(lnwage ~ public)
out.ols1 <-summ(ols1, , robust = "HC1")
out.ols1


# Conditional OLS 
# ----------------------------------------------------

ols2 <- lm(lnwage ~ public + covariates)
out.ols2 <-summ(ols2, , robust = "HC1")
out.ols2



# Unconditional median regression
# ----------------------------------------------------

rq1 <- rq(lnwage ~ public, tau=.5, data=data, method="br", model = TRUE) 
out.rq1 <- summary(rq1)
out.rq1

median_public <- quantile(lnwage[public==1], .5)
median_public
median_private <- quantile(lnwage[public==0], .5)
median_private


# Conditional quantile 
# ----------------------------------------------------

covariates <- as.matrix(covariates) # make sure covariates are a matrix 
rq2 <- rq(lnwage ~ public + covariates, tau=.5, data=data, method="br", model = TRUE) 
out.rq2 <- summary(rq2)
out.rq2



# Conditional quantile regression, across different quantiles
# ----------------------------------------------------

# Vector of quantiles to be estimated has to be strictly within (0,1) interval 
quants <- as.matrix(seq(0.1, 0.9, by=0.2)) 

rq3 <- rq(lnwage ~ public + covariates, tau=quants, data=data, method="br", model = TRUE) 
out.rq3 <- summary(rq3)
out.rq3

# Plot quantile effects
plot(out.rq3)

# see Koenker (2001, JEP) for interpretation


####################################################################
#   PHASE 4: Counterfactuals
####################################################################

# either: counterfactual population is an artificial transformation of a reference population
# or: reference and counterfactual populations correspond to two different groups (e.g. treatment and control group)

#cf_function <- counterfactual(formula, data, weights, na.action = na.exclude,
							 # group, treatment = FALSE, decomposition = FALSE,
						## for transformation of reference population
							 # transformation = FALSE, counterfactual_var,
							 # quantiles, 

						## method (quantreg):
						    # method = "qr", 
							 # trimming = 0.005, nreg = 100, 

						## for other methods (loc, locsca, cox, lgit, probit): 
							# scale_variable,
							 # counterfactual_scale_variable,
							 # censoring = 0, right = FALSE, nsteps = 3,
							 # firstc = 0.1, secondc = 0.05, 
						## for inference: 
							 # noboot = FALSE,
							 # weightedboot = FALSE, seed = 8, robust = FALSE,
							 # reps = 100, alpha = 0.05, 
							 # first = 0.1, last = 0.9,	# subset of quantile indexes of interest, tails should not be used
							 # cons_test = 0, 
						# options for parallel computing
							 # printdeco = TRUE,
							 # sepcore = FALSE, ncore=1)

taus <- c(1:99)/100
first <- sum(as.double(taus <= .10))
last <- sum(as.double(taus <= .90))
rang <- c(first:last)

# tails of the distribution should not be used because standard asymptotic does not apply to these parts.

rq.4 <-  counterfactual(lnwage ~ covariates, data=data, group=public, quantiles=taus,
						treatment = TRUE, decomposition = TRUE,  reps = 3,  sepcore = TRUE, ncore=2) 



duqf_SE <- (rq.4$resSE)[,1]
l.duqf_SE <- (rq.4$resSE)[,3]
u.duqf_SE <- (rq.4$resSE)[,4]

duqf_CE <- (rq.4$resCE)[,1]
l.duqf_CE <- (rq.4$resCE)[,3]
u.duqf_CE <- (rq.4$resCE)[,4]

duqf_TE <- (rq.4$resTE)[,1]
l.duqf_TE <- (rq.4$resTE)[,3]
u.duqf_TE <- (rq.4$resTE)[,4]

range_x <- min(c(min(l.duqf_SE[rang]), min(l.duqf_CE[rang]),min(l.duqf_TE[rang])))
range_y <- max(c(max(u.duqf_SE[rang]), max(u.duqf_CE[rang]),max(u.duqf_TE[rang])))

par(mfrow=c(1,3))

plot(c(0,1), range(c(range_x, range_y)), xlim = c(0,1), type = "n",
xlab = expression(tau), ylab = "Difference in Wages", cex.lab=0.75, 
main = "Total")
polygon(c(taus[rang],rev(taus[rang])),
c(u.duqf_TE[rang], rev(l.duqf_TE[rang])), density = -100, border = F, col = "grey70", lty = 1, lwd = 1)
lines(taus[rang], duqf_TE[rang])
abline(h = 0, lty = 2)
plot(c(0,1), range(c(range_x, range_y)), xlim = c(0,1), type = "n",
xlab = expression(tau), ylab = "", cex.lab=0.75, main = "Structure")

polygon(c(taus[rang],rev(taus[rang])),
c(u.duqf_SE[rang], rev(l.duqf_SE[rang])), density = -100, border = F,
col = "grey70", lty = 1, lwd = 1)
lines(taus[rang], duqf_SE[rang])
abline(h = 0, lty = 2)
plot(c(0,1), range(c(range_x, range_y)), xlim = c(0,1), type = "n",
xlab = expression(tau), ylab = "", cex.lab=0.75, main = "Composition")

polygon(c(taus[rang],rev(taus[rang])),
c(u.duqf_CE[rang], rev(l.duqf_CE[rang])), density = -100, border = F,
col = "grey70", lty = 1, lwd = 1)
lines(taus[rang], duqf_CE[rang])
abline(h = 0, lty = 2)


p1 <- ggplot(mp, aes(year, wow))+
    geom_point()+
    geom_line(data=predframe)+
    geom_ribbon(data=predframe,aes(ymin=lwr,ymax=upr),alpha=0.3)

####################################################################
#  Solution to small assignment
####################################################################

# Men and women ??? as in Melly 2005 



