### AEA2020 ###


# Also provide Stata code 


####################################################################
#   PHASE 0: LOAD PACKAGED AND READ IN DATA  
####################################################################


# remove old objects from R working space
rm(list=ls())

# Define packages that you need
packages_vector <- c("tidyverse", "foreign", "haven", "Hmisc", "dplyr", "tidyr", "stringr", "dummies", "sandwich", "lmtest",
     "expss", "fastDummies", "psych", "fBasics", "knitr", "xtable", "sjlabelled", "data.table", "lubridate", "arsenal", 
     "stargazer", "mfx", "jtools")
# install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) 
# List loaded packages 
(.packages())


# # List packages separately
# library(foreign) # only up until stata version 12
# library(haven)	# older than 13 
# library(Hmisc)
# library (dplyr)
# library(tidyr)
# library(stringr)


# require(knitr)
# require(survival)



# Working Directory
getwd()
work_dir <- "Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/First Session"
setwd(work_dir)


# Read in the data
data <-read.dta("sobs_All_v13.dta") # foreign, for Stata versions 5-12
data <-read_dta("sobs_All_v13.dta") # haven
data <- as.data.frame(data)

# Attach variables
# attach(data)

####################################################################
#   PHASE 1: INSPECT YOUR DATA   
####################################################################

#1. Understand the Structure of your Data

# Check the class of the dataset. 
class(data)

# View the column names to get a first sense of what you have. 
names(data)

# See a compact summary of the data
glimpse(data) 
# str(data) # less compact

# Note R encodes missing values as NAs. Look for other common 
# missing values encodings in your data (-1, 9999, N/A).

# Have a look at the first  or last 6 rows of a dataframe. - alternative to viewing the data frame 
head(data) 
# tail() 

# do you have labels?  - only for single variables
var_lab(data$id)
var_lab(data$sex)

# Browse data by condition
head(data[data$age>59,])

# Make list of categorical and continuous variables  
xcont_names <- c("age", "insured_earn", "contr_5y", "unempl_r") 
xcat_names <- c("sex", "marits", "region", "swiss", "educ", "full_time", "lastj_occpt", "lastj_fct")

# Summary of the distribution of each variable
# continuous variables
desc <- fBasics::basicStats(data[xcont_names]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))
# t() is to transpose the table
# does not work for dates
# shows how many missings you have
# for certain functions you might have to specify the package name before the command name package::command()
# here we use piping to apply a number of functions

# categorical variables
desc <- fBasics::basicStats(data[xcat_names]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))
# here we overwrite the stored descriptives

####################################################################
#   PHASE 2: CLEANING  
####################################################################

#########################################
#   Duplicates, missings and identifiers
#########################################

# Check the number of observations
nobs <- nrow(data)
nobs

# Sort rows by values of a column or columns 
data <- arrange(data, data$id, data$date_start) 

# Remove rows with duplicate values in all variables
data <- distinct(data) 
nrow(data)

# Check missings (redundant w.r.t. summary stats)
sum(is.na(data$treat))

# Drop rows containing NA’s in important variables (tidyr)
data <- drop_na(data, id, treat, date_start)
nrow(data)

# Replace missing end date with a maximum duration of 2 years after program start
maxdur <- 24
data$date_end[is.na(data$date_end)] <- data$date_start[is.na(data$date_end)] + maxdur*30

# Drop end dates that are before start dates 
sum(data$date_end < data$date_start)
data <- dplyr::filter(data, date_end> date_start)

# Create observation identifier (id is)
data$idobs <- rep(1:nrow(data))

# Create spell identifier within person identifier 
data <- data %>% group_by(id) %>% mutate(spell = sequence(n())) %>% ungroup()
data <- data %>% group_by(id) %>% mutate(num_spells = max(spell)) %>% ungroup()

# Order columns
data_order1 <- dplyr::select(data, id, idobs, spell, num_spells) 
data_order2 <- dplyr::select(data, -id, -idobs, -spell, -num_spells)
data <- cbind(data_order1, data_order2)
rm(list="data_order1", "data_order2") # Drop these frames


##################################
#  Clean continuous covariates 
##################################

# Summary of the distribution of each variable
desc <- fBasics::basicStats(data[xcont_names]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))

# Create dummy flagging observations with missing 
data$missing_unempl_r <- ifelse(is.na(data$unempl_r), 1, 0)

# Code missings as mean of variable by group - continuous, remember to specify that missings should be ignored
data <- data %>% group_by(region) %>% mutate(mean_unempl_r = mean(unempl_r, na.rm=TRUE)) %>% ungroup()
is.grouped_df(data)
data$unempl_r[is.na(data$unempl_r)] <- data$mean_unempl_r[is.na(data$unempl_r)]

# Histogram - see if outliers, if values are capped, min max etc. 
ggplot(data, aes(x = age)) +
        geom_histogram(aes(y = ..count..), binwidth = 5, col="red", fill="blue", alpha = .2)
ggplot(data, aes(x = insured_earn)) +
        geom_histogram(aes(y = ..count..), col="red", fill="blue", alpha = .2)
ggplot(data, aes(x = contr_5y)) +
        geom_histogram(aes(y = ..count..), col="red", fill="blue", alpha = .2)
ggplot(data, aes(x = unempl_r)) +
        geom_histogram(aes(y = ..count..), col="red", fill="blue", alpha = .2)

#Find outliers and check 
outlier <- which.max(data$age) #Which row contains the max value for age
data[outlier,]

# Drop those above 60 - Extract rows that meet logical criteria.
sum(data$age<18 & data$age >= 60)
data <- dplyr::filter(data, age>=18 & age < 60)	

# Age groups
data$agegr<-NA
data$agegr[data$age>=18 & data$age < 30] <- 1
data$agegr[data$age>=30 & data$age < 40] <- 2
data$agegr[data$age>=40 & data$age < 50] <- 3
data$agegr[data$age>=50 & data$age < 60] <- 4
data$agegr[data$age>=60] <- 5

##################################
#  Clean categorical covariates
##################################

# Update categorical variables names, including categorical created from continuous
xcat_names <- c(xcat_names, "agegr")

# Summary statistics - repeat 
desc <- fBasics::basicStats(data[xcat_names]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))
# Emphasize that mean does not make sense, except for dummies, mostly look at min max 

# Tabulate with labels (if exist)
cro(data$sex)
cro(data$agegr)
cro(data$marits)
cro(data$region)
cro(data$swiss)
cro(data$educ) 
cro(data$full_time)
cro(data$lastj_occpt)
cro(data$lastj_fct)

# Cross-tabulate
cro(data$sex, data$agegr)

# Check missings by treatment group, example with education
sum(is.na(data$educ[data$treat==1]))
sum(is.na(data$educ[data$treat==0]))

# Code missings as separate category with a specified value, if many observations missing - categorical
data<-replace_na(data, replace=list(educ=99))
data<-replace_na(data, replace=list(lastj_fct=99))

# Create dummies for all categorical variables
xcat_desc <-dummy_cols(data[xcat_names], select_columns = c(xcat_names), remove_selected_columns = TRUE, remove_first_dummy = FALSE, remove_most_frequent_dummy = FALSE)
# Leave one dummy out (most frequent category) for regressions
xcat_reg  <-dummy_cols(data[xcat_names], select_columns = c(xcat_names), remove_selected_columns = TRUE, remove_first_dummy = FALSE, remove_most_frequent_dummy = TRUE)

# Labelling example -- Age group
var_lab(data$agegr) = "Age group"
val_lab(data$agegr) = num_lab("
            1 Age 18-29
            2 Age 30-39
            3 Age 40-49
            4 Age 50-59
            5 Age 60+ ")

var_lab(data$agegr)

# Bar plot example -- Gender by age group 
barplot(
  table(
        sjlabelled::as_label(data$sex),
        sjlabelled::as_label(data$agegr)),
  legend.text = TRUE, beside=TRUE
)

# ???
# drop_val_labs(data$agegr)
# check missing labels in xaxis

# Same thing
# data$educ<-as.factor(data$educ)
# plot(data$educ) # nice histogram for factor vars


##################################
#  Generate outcomes
##################################

# Total unemployment spell duration
data$duration <- as.numeric(data$date_end - data$date_start + 1)

# Set maximum horizon to look at after start of unemployment, = 2 years and censor duration for those with long spells and missing end dates (never found a job)
maxdur <- 24
data$duration[data$duration>maxdur*30] <- maxdur*30
data$duration[is.na(data$duration)] <- maxdur*30

# Employed within 1 year after program start/unemployment
data$employed1y <- ifelse(data$duration<366,1,0)

####################################################################
#   PHASE 3: CHECK CLEANING  
####################################################################

# Regroup variables
dataid <- dplyr::select(data, id, idobs, date_end, date_start, spell, num_spells) 
treat <- dplyr::select(data, treat)
outcomes <- dplyr::select(data, duration, employed1y)
xcont_clean <- dplyr::select(data, insured_earn, contr_5y, unempl_r)

# Create vector of categorical variable names for descriptives and regression
xcat_names_desc <- colnames(xcat_desc)
xcat_names_reg <- colnames(xcat_reg) 
xcat_names_desc
xcat_names_reg

# Update other vectors of names
xcont_names <- colnames(xcont_clean)

# Create new clean data frames - 1 for descriptives with all dummies, 1 for regression without most frequent one
data_desc <- cbind(treat, outcomes, xcont_clean, xcat_desc) 
data_reg <- cbind(dataid, treat, outcomes, xcont_clean, xcat_reg)

# Table with Summary Statistics  
desc <- fBasics::basicStats(data_desc) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2)) 

# Export as Latex-file
# xtable, kable, stargazer

# Histogram for duration 
ggplot(data_desc, aes(x = duration)) +
        geom_histogram(aes(y = ..density..), binwidth = 30, col="red", fill="blue", alpha = .2)

# Check data is ungrouped - remove later
# is.grouped_df(data_reg)
# is.grouped_df(data_desc)

save(data_desc, xcat_names_desc, xcont_names, file="data_desc.RData")
save(data_reg, xcat_names_reg, xcont_names, file="data_reg.RData")



####################################################################
 # DATA CLEANING PHASE ENDS HERE, APPLICATION STARTS HERE
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################



####################################################################
#   PHASE 4: DESCRIPTIVE STATISTICS   
####################################################################


load("data_desc.RData")

data_desc$treat <- as.logical(data_desc$treat)
# Specify formulas for summary statistics by treatment group
f_sum_cov= as.formula(paste("treat", paste(xcont_names, paste(xcat_names_desc, collapse = " + "), sep = " + "), sep=" ~ " )) 
f_sum_cov

# Covariates
sumstats_byd_cov<-summary(tableby(f_sum_cov, data = data_desc, control=tableby.control(numeric.stats=c("meansd")), total = FALSE))
sumstats_byd_cov<-summary(tableby(f_sum_cov, data = data_desc,  control=tableby.control(numeric.stats=c("meansd")), total = FALSE))
# Outcomes
sumstats_byd_out<-summary(tableby(treat ~ duration, data = data_desc, control=tableby.control(numeric.stats=c("meansd")), test=FALSE))

sumstats_byd_cov
sumstats_byd_out

desc <- fBasics::basicStats(data$duration[data$treat==1]) %>% t() %>% as.data.frame() %>% 
  dplyr::select(Mean, Stdev, Minimum, Maximum, nobs)
print(round(desc, digits=1))

desc <- fBasics::basicStats(data$duration[data$treat==0]) %>% t() %>% as.data.frame() %>% 
  dplyr::select(Mean, Stdev, Minimum, Maximum, nobs)
print(round(desc, digits=1))



####################################################################
#   PHASE 5: MODEL SPECIFICATION WITH CROSS-SECTIONAL TREATMENT
####################################################################


load("data_reg.RData")

# Create application-specific variables that might be needed

# Quarter of entry into unemployment 
data_reg$quarter <- as.numeric(quarter(data_reg$date_start))
data_reg<-dummy_cols(data_reg, select_columns = c("quarter"), remove_most_frequent_dummy = TRUE)

# save final data set
save(data_reg, xcat_names_reg, xcont_names, file="data_reg_final.RData")


# Group variables, adding newly created ones - need matrix!
y1 <- as.matrix(data_reg$duration)
y2 <- as.matrix(data_reg$employed1y)
treat <- as.matrix(data_reg$treat)
x1 <- as.matrix(dplyr::select(data_reg, xcat_names_reg, xcont_names, starts_with("quarter_")))

dim(y1)
dim(y2)
dim(x1)

# Alternatively, specify formula for regression
# x1_reg <- paste(colnames(x1), collapse=" + ")
# x1_reg


##################################
# Single outcome
##################################

# Linear regression

# Linear model - effect on duration
lm.model.y1 <- lm(y1 ~ treat + x1)
cov <- vcovHC(lm.model.y1, type = "HC")
robust.se.y1 <- sqrt(diag(cov))

# Linear model - effect on employment in year 1
lm.model.y2 <- lm(y2 ~ treat + x1)
cov <- vcovHC(lm.model.y2, type = "HC")
robust.se.y2 <- sqrt(diag(cov))

# Probit model - effect on employment in year 1
probit.model <- glm(y2 ~ treat + x1, family = binomial(link = "probit") )
cov <- vcovHC(probit.model, type = "HC")
robust.se.probit <- sqrt(diag(cov))
# Marginal effect
probit.model.2 <- probitmfx(y2 ~ treat + x1, data=data_reg, atmean = TRUE, robust = TRUE)

# Output Coefficients
stargazer(lm.model.y1,lm.model.y2,probit.model, se=list(robust.se.y1, robust.se.y2, robust.se.probit), 
            keep=c("treat"), keep.stat = c("n", "rsq"),
              column.labels=c("Duration","Emply1", "Emply1"), align=TRUE, dep.var.labels.include = FALSE)

# Output ME
stargazer(lm.model.y1, lm.model.y2, probit.model.2$fit, coef = list(NULL, NULL, probit.model.2$mfxest[,1]), se=list(robust.se.y1, robust.se.y2 , probit.model.2$mfxest[,2]), 
            keep=c("treat"), keep.stat = c("n", "rsq"),
              column.labels=c("Duration","Emply1", "Emply1"), align=TRUE, dep.var.labels.include = FALSE)




############################################
# Small assignment (first part)
############################################

# Table with Summary Statistics  by gender
desc_m <- fBasics::basicStats(data_desc[data_desc$sex_0==1, -c(7:8)]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc_m, digits=2)) 
desc_f <- fBasics::basicStats(data_desc[data_desc$sex_0==0, -c(7:8)]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc_f, digits=2)) 

# or
xcat_bysex <- xcat_desc[c(-1,-2)]
xcat_names_bysex <- colnames(xcat_bysex)
f_sum_cov= as.formula(paste("sex_1", paste(xcont_names, paste(xcat_names_bysex, collapse = " + "), sep = " + "), sep=" ~ " )) 
f_sum_cov

sumstats_bysex_cov<-summary(tableby(f_sum_cov, data = data_desc, control=tableby.control(numeric.stats=c("meansd")), test=FALSE))
sumstats_bysex_cov

# Regression for females
# Linear model - effect on duration
x1 <- as.matrix(dplyr::select(data_reg, xcat_names_reg, xcont_names, -sex_1, starts_with("quarter_")))
lm.model.y1.f <- lm(duration[females==1] ~ treat[females==1] + x1[females==1,] )

data_test <-cbind(data_reg$duration, data_reg$duration)
lm.model.y1.f <- lm(duration ~ treat + x1, data=data_reg[females==1] )

cov <- vcovHC(lm.model.y1.f, type = "HC")
robust.se.y1.f <- sqrt(diag(cov))

# if variables are put in separately and not as matrices you can also specify a subset of the data 

# Output Coefficients
stargazer(lm.model.y1.f, se=list(robust.se.y1.f), 
            keep=c("treat"), keep.stat = c("n", "rsq"),
            column.labels=c("Duration (Females)"), align=TRUE, dep.var.labels.include = FALSE)
