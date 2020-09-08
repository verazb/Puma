### AEA2020 ###

# Regression discontinuity design

# CIT refers to Cattaneo Idrobo Titunik, 2019 


####################################################################
#   PHASE 0: LOAD PACKAGES AND READ IN DATA  
####################################################################

rm(list=ls())


packages_vector <- c("haven", "dplyr", "tidyr", "sandwich", "expss",
    "fBasics", "xtable", "data.table", "stargazer", "mfx", "jtools", "ggplot2")
# install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) 


# RDD-specific packages 
packaged_vector_rdd <- c("grid", "lpdensity", "rddensity", "rdlocrand", "rdrobust", "TeachingDemos")
# install.packages(packages_vector)
lapply(packaged_vector_rdd, require, character.only = TRUE) 


# List loaded packages 
(.packages())


# Working Directory
work_dir <- "Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/RDD/Meyersson 2014/Prep/"
setwd(work_dir)


# Read in the data
data <- as.data.frame(read_dta("data_meyersson.dta")) 



####################################################################
#   PHASE 1: INSPECT YOUR DATA AND DEFINE KEY VARIABLES
####################################################################

# Inspect data
head(data)


# Vector with all variable names
varnames <- colnames(data)


# store each variable in own R object
attach(data)


# X               Islamic Margin of Victory
# Y               Female High School percentage 
# D               Islamic rule 
# ageshr19        Percentage of population below 19 in 2000
# ageshr60        Percentage of population above 60 in 2000
# buyuk           Metro center
# hischshr1520m   Percentage of men aged 15-20 with high school
#                                                 education
# i89              Islamic Mayor in 1989
# lpop1994         Log population in 1994
# merkezi          District center
# merkezp          Province center
# partycount       Number of parties receiving votes 1994
# prov_num         Province number 
# sexr             Gender ratio in 2000
# shhs             Household size in 2000
# subbuyuk         Sub-metro center
# vshr_islam1994   Islamic percentage of votes in 1994


# Options for RD plots
options(width=280)
par(mar = rep(2, 4))
xlabel <- "Islamic Margin of Victory"
ylabel <- "Female High School percentage"
dlabel <- "Islamic rule"



####################################################################
#   PHASE 2: DESCRIPTIVE STATISTICS AND VALIDITY CHECKS
####################################################################


# Number of treated and control
cro(D)

# Table 1 in Meyersson
# Summary statistics - Pooled 
desc <- fBasics::basicStats(data) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))

# Summary statistics - Treated 
desc_treat <- dplyr::filter(data, D==1) %>% fBasics::basicStats() %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc_treat, digits=2))

# Summary statistics - Control 
desc_control <- dplyr::filter(data, D==0) %>% fBasics::basicStats() %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc_control, digits=2))

# Replace this with balancedness table !!!


ttest <- t.test(Y ~ D, data = data)
ttest$estimate[[1]]
ttest$estimate[[2]]
ttest$p.value



# Existence of discontinuity: Plot D against X
# ----------------------------------------------------

ggplot(data,                
       aes(x = X, y = D)) +
 	   geom_line() +
 	   ylab("Effect on employment rate") +
 	   theme_bw(base_size = 20) 



# Descriptive evidence of effect: plot Y against X
# ----------------------------------------------------

# Figure 3a in CIT
# Raw comparison of means (no polynomial)
rdplot(Y, X, nbins = c(2500, 500), p = 0, title = "", x.label = xlabel, y.label = ylabel)


# Figure 3b in CIT
# Local comparison of means within 50 percentage points bandwidth, but same number of bins
rdplot(Y[abs(X) <= 50], X[abs(X) <= 50], nbins = c(2500, 500), p = 4, title = "", x.label = xlabel, y.label = ylabel)




# No sorting in X: plot histogram and density of X 
# ----------------------------------------------------

# Histogram 
# Figure 19a in CIT 2019, Figure 2(a) in Meyersson 2014

# Specify bandwidth manually 
# Useful to plot across the full range of possible values of the running variable 
bw <- as.numeric(100)

plot1 = ggplot(data=data, aes(X)) + 
  geom_histogram(data = data, aes(x = X, y= ..count..), breaks = seq(-bw, 0, 2), fill = "blue", col = "black", alpha = 1) +
  geom_histogram(data = data, aes(x = X, y= ..count..), breaks = seq(0, bw, 2), fill = "red", col = "black", alpha = 1) +
  labs(x = xlabel, y = "Number of Observations") + 
  geom_vline(xintercept = 0, color = "black") +
  theme_bw(base_size = 20) 
plot1


# Estimated Density
# Figure 19b in CIT, Figure 2(b) in Meyersson 

# Specify bandwidth manually, within a bandwidth that makes sense depending on histogram (i.e. where is mass)
bw = as.numeric(30)


est1 = lpdensity(data = X[X < 0 & X >= -bw], grid = seq(-bw, 0, 0.1), bwselect = "IMSE",
                 scale = sum(X < 0 & X >= -bw) / length(X))
est2 = lpdensity(data = X[X >= 0 & X <= bw], grid = seq(0, bw, 0.1), bwselect = "IMSE",
                 scale = sum(X >= 0 & X <= bw) / length(X))
plot2 = lpdensity.plot(est1, est2, CIshade = 0.2, lcol = c(4, 2), CIcol = c(4, 2), legendGroups = c("Control", "Treatment"))+
  labs(x = xlabel, y = "Density") + 
  geom_vline(xintercept = 0, color = "black") +
  theme_bw(base_size = 20) + 
  theme(legend.position = c(0.8, 0.85))
plot2



# No sorting in X: test for discontinuities in density 
# ----------------------------------------------------

# RD Manipulation Test using local polynomial density estimation.
out = rddensity(X)
summary(out)



# No sorting in covariates: plot covariates against X
# ----------------------------------------------------

# Figure 16 in CIT, Figure 3 in Meyersson
# Just plots, do not serve as formal inference

bw <- as.numeric(100)

plot_covariate <-rdplot(vshr_islam1994, X, h = bw, x.label = xlabel, y.label = "vshr_islam1994", title = "")
plot_ageshr60 <-rdplot(ageshr60, X, h = bw, x.label = xlabel, y.label = "ageshr60", title = "")
plot_ageshr19 <-rdplot(ageshr19, X, h = bw, x.label = xlabel, y.label = "ageshr19", title = "")
plot_lpop1994 <-rdplot(lpop1994, X, h = bw, x.label = xlabel, y.label = "lpop1994", title = "")
plot_sexr <-rdplot(sexr, X, h = bw, x.label = xlabel, y.label = "sexr", title = "")
plot_partycount <-rdplot(partycount, X, h = bw, x.label = xlabel, y.label = "partycount", title = "")



# No sorting in covariates: covariates as outcome
# ----------------------------------------------------

# This part serves for formal inference -- see RDD estimation part below 


# Select covariates (see Meyersson p. 246)

# Create dummies for province fixed effects -- drop these??
# prov_fe <- as.matrix(dummy_cols(data["prov_num"], select_columns = c("prov_num"), remove_selected_columns = TRUE, remove_first_dummy = FALSE, remove_most_frequent_dummy = TRUE))

covariates <- as.matrix(dplyr::select(data, lpop1994, vshr_islam1994, partycount, ageshr60, ageshr19, sexr, shhs, merkezi, merkezp, subbuyuk, buyuk))
# covariates <- cbind(covariates, prov_fe)
covariates <- cbind(covariates)
covariates_names <- colnames(covariates)

covariate_rd <- function(covariate) {

	out <- rdrobust(covariate, X)
	summary(out)
	list(effect=out$Estimate[[1]], se=out$Estimate[[3]])  

}

covariate_out <- apply(covariates, 2, covariate_rd)


# convert list to table
covariate_table <- rbindlist(covariate_out)
covariate_table


# Export table with all estimates 

rownames(covariate_table) <- covariates_names
xtable(covariate_table, digits =3)



####################################################################
#   PHASE 3: RDD effect estimation and robustness checks
####################################################################

# Create matrix to store effects
effects <- matrix(, nrow = 5, ncol = 8)


# Can run different combinations of all of these
# Optimal bandwidth choice algorithm differs from Meyersson


# Table II Panel A (Women), Meyersson

# Global OLS estimation without covariates (raw difference in means using all observations) 
# Col (1)
# ----------------------------------------------------

out <- lm(Y ~ D)
out <- summ(out)
out
effects[1,1] <- out$coeftable[2,1] 		# effect
effects[2,1] <- out$coeftable[2,2] 		# se
effects[3,1] <- 100 					# bandwidth left 
effects[4,1] <- 100 					# bandwidth right 
effects[5,1] <- nrow(out$model$model) 	# total number of observations



# Parametric estimation with covariates (using all observations, i.e. global estimation)
# ----------------------------------------------------


# Col (2)
out <- lm(Y ~ D + covariates)
out <- summ(out)
out
effects[1,2] <- out$coeftable[2,1]
effects[2,2] <- out$coeftable[2,2] 
effects[3,2] <- 100
effects[4,2] <- 100
effects[5,2] <- nrow(out$model$model)



# Local linear regression (nonparametric as in CIT, vs. Meyersson, parametric estimation)
# ----------------------------------------------------

# using all observations
out <- rdrobust(Y, X, kernel = "triangular", scaleregul = 1, p = 1, h = 100)
summary(out)


# Col (3) with optimal bandwidth selection 
out <- rdrobust(Y, X, kernel = "triangular", scaleregul = 1, p = 1, bwselect = "mserd")

summary(out)
effects[1,3] <- out$Estimate[[1]]
effects[2,3] <- out$Estimate[[3]]
effects[3,3] <- out$bws[[1]]
effects[4,3] <- out$bws[[2]]
effects[5,3] <- out$N[[1]] + out$N[[2]]


# Col (4) with covariates
out <- rdrobust(Y, X, covs = covariates, kernel = "triangular", scaleregul = 1, p = 1, bwselect = "mserd")

summary(out)
effects[1,4] <- out$Estimate[[1]]
effects[2,4] <- out$Estimate[[3]]
effects[3,4] <- out$bws[[1]]
effects[4,4] <- out$bws[[2]]
effects[5,4] <- out$N[[1]] + out$N[[2]]


# Vary sample around cutoff 
# ----------------------------------------------------

# Store h from previous estimation
bw_left <- as.numeric(out$bws[[1]])
bw_right <- as.numeric(out$bws[[2]])


# Col (5) with h/2
out = rdrobust(Y, X, h = c(bw_left/2, bw_right/2), kernel = "triangular", p = 1)

summary(out)
effects[1,5] <- out$Estimate[[1]]
effects[2,5] <- out$Estimate[[3]]
effects[3,5] <- out$bws[[1]]
effects[4,5] <- out$bws[[2]]
effects[5,5] <- out$N[[1]] + out$N[[2]]


# Col (6) with 2h
out = rdrobust(Y, X, h = c(bw_left*2, bw_right*2), kernel = "triangular", p = 1)

summary(out)
effects[1,6] <- out$Estimate[[1]] 
effects[2,6] <- out$Estimate[[3]] 
effects[3,6] <- out$bws[[1]] 
effects[4,6] <- out$bws[[2]] 
effects[5,6] <- out$N[[1]] + out$N[[2]] 



# Higher order polynomials of (Xi - x) 
# ----------------------------------------------------

# Col (7), quadratic
out = rdrobust(Y, X, covs = covariates, kernel = "triangular", p = 2, bwselect = "mserd")

summary(out)
effects[1,7] <- out$Estimate[[1]]
effects[2,7] <- out$Estimate[[3]]
effects[3,7] <- out$bws[[1]]
effects[4,7] <- out$bws[[2]]
effects[5,7] <- out$N[[1]] + out$N[[2]]



# Col (8), cubic
out = rdrobust(Y, X, covs = covariates, kernel = "triangular", p = 3, bwselect = "mserd")

summary(out)
effects[1,8] <- out$Estimate[[1]]
effects[2,8] <- out$Estimate[[3]]
effects[3,8] <- out$bws[[1]]
effects[4,8] <- out$bws[[2]]
effects[5,8] <- out$N[[1]] + out$N[[2]]



# Export table with all estimates 

rownames(effects) <- c("Effect", "SE", "Bw left", "Bw right", "Total obs")
xtable(effects, digits =3)





####################################################################
#   EXTENSION 1: EFFECT HETEROGENEITY
####################################################################

# Table V Panel C, Meyersson 

# Restrict sample to bandwidth of 25
# Only after, split at median of Islam vote share
median <- as.numeric(median(vshr_islam1994[abs(X)<=25]))
median

# Col (1) above median 
out = rdrobust(Y[vshr_islam1994 >= median], X[vshr_islam1994 >= median], kernel = "triangular", p = 1, h = 25)
summary(out)


# Col (2) below median 
out = rdrobust(Y[vshr_islam1994 < median], X[vshr_islam1994 < median], kernel = "triangular", p = 1, h = 25)
summary(out)




####################################################################
#  EXTENSION 2: DONUT HOLE APPROACH 
####################################################################


# Test for excluding different ranges of observations right next to the cutoff 
# ----------------------------------------------------

donut_range <- as.matrix(c(0, 0.1, 0.2, 0.3, 0.4, 0.5))

donut_tests <- function(donut) {          
	
	if (is.numeric(donut)) {
		out = rdrobust(Y[abs(X) >= donut], X[abs(X) >= donut], kernel = "triangular", p = 1)
		list(effect=out$Estimate[[1]], se=out$Estimate[[3]])  
	} else {
		print("ERROR - Donut range must be numeric")
	}

}

donut_out <- apply(donut_range, 1, donut_tests)


# convert list to table
donut_table <- rbindlist(donut_out)
donut_table




####################################################################
#  Small assignment solution Question 3 
####################################################################


# Test for jumps in Y at non-discontinuity points (placebo cutoffs)
# ----------------------------------------------------

placebo_cutoffs <- as.matrix(c(-3, -2, -1, 0, 1, 2, 3))

placebo_tests <- function(cutoff) {          
	
	if (is.numeric(cutoff)) {
		if (cutoff < 0) {
			out = rdrobust(Y[X < 0], X[X < 0], c = cutoff)
		} else if (cutoff == 0) {
			out = rdrobust(Y, X, c = cutoff)
		} else if (cutoff > 0) {
			out = rdrobust(Y[X >= 0], X[X >= 0], c = cutoff)
		}
			list(effect=out$Estimate[[1]], se=out$Estimate[[3]]) 
	} else {
		print("ERROR - Cutoff must be numeric")
	}

}

placebo_out <- apply(placebo_cutoffs, 1, placebo_tests)


# convert list to table
placebo_table <- rbindlist(placebo_out)
placebo_table





