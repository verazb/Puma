### AEA2020 ###

# Difference-in-Differences

# HP refers to Hjort and Poulsen, 2019 


####################################################################
#   PHASE 0: LOAD PACKAGES AND READ IN DATA  
####################################################################

rm(list=ls())


packages_vector <- c("haven", "dplyr", "tidyr", "sandwich", "expss",
    "fBasics", "xtable", "data.table", "stargazer", "mfx", "jtools", "ggplot2")
# install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) 


# DiD-specific packages 
packaged_vector_did <- c("fixest")
# install.packages(packages_vector)
lapply(packaged_vector_did, require, character.only = TRUE) 


# List loaded packages 
(.packages())


# Working Directory
work_dir <- "Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/DiD/Prep"
setwd(work_dir)


# Read in the data
# In this replication, use the QLFS only!
data <- as.data.frame(read_dta("data.dta")) 



####################################################################
#   PHASE 1: INSPECT YOUR DATA AND DEFINE KEY VARIABLES
####################################################################

# Inspect data
head(data)

# Treated group variable
radius <- 556.6 # preferred connection radius of 556.6m, i.e. 0.005 decimal degrees
data$connected <- ifelse(data$distance <= radius, 1, 0)

# Treatment period variable (after 2009Q3)
data$submarines <- ifelse(data$quarter>=20093, 1, 0)

# Treatment variable (treated x post)
data$treatment <- data$connected * data$submarines

# Attach data 
attach(data)


####################################################################
#   PHASE 2: DESCRIPTIVE STATISTICS AND VALIDITY CHECKS
####################################################################

did_vars <- c("connected", "submarines", "treatment")

# Number of treated and control, before and after treatment
desc <- fBasics::basicStats(data[did_vars]) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))

# Overall
cro(connected, treatment)

# By quarter
cro(quarter, treatment)


# Table 1 - Summary statistics of outcomes by treatment group in 2008Q1
# ----------------------------------------------------

outcomes <- c("employed", "skilled", "hoursworked", "workmore", "formal", "informal")

# Summary statistics - All 
desc <- fBasics::basicStats(data[outcomes], quarter==20081) %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc, digits=2))

# Summary statistics - Treated 
desc_treat <- dplyr::filter(data[outcomes], connected==1 & quarter==20081) %>% fBasics::basicStats() %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc_treat, digits=2))

# Summary statistics - Control 
desc_control <- dplyr::filter(data[outcomes], connected==0 & quarter==20081) %>% fBasics::basicStats() %>% t() %>% as.data.frame() %>% dplyr::select(Mean, Stdev, Minimum, Maximum, nobs, NAs)
print(round(desc_control, digits=2))

# Replace this with balancedness table !!!


# Common trends graph (unadjusted means of outcome by quarter)
# Figure 6 of HP
# ----------------------------------------------------

# Create a group-means data set
common_trends <- data %>% group_by(time, connected) %>% summarise(empl_mean = mean(employed))
is.grouped_df(data)

ggplot(data = common_trends, aes(x = time, y = empl_mean, group = connected, color = connected)) + 
  geom_line() +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_bw(base_size = 20)

  # color legend
  # has 2 additional periods before compared to paper


# Alternative outcome - hours worked
common_trends <- data %>% group_by(time, connected) %>% summarise(loghours_mean = mean(loghours))
is.grouped_df(data)

ggplot(data = common_trends, aes(x = time, y = loghours_mean, group = connected, color = connected)) + 
  geom_line() +
  geom_vline(xintercept = 0, linetype="dashed") +
  theme_bw(base_size = 20)



####################################################################
#   PHASE 3: DiD effect estimation 
####################################################################

# Using the fixest package
# For more details, see https://cran.r-project.org/web/packages/fixest/vignettes/fixest_walkthrough.html


# Standard DiD estimation (raw difference without fixed effects)
# ----------------------------------------------------

did1 <- feols(employed ~ connected + submarines + treatment, data)
did1 <- summary(did1, cluster = "location")
did1


# Standard DiD estimation (with quarter fixed effects only)

# Quarter fixed effects absorb the submarine variable
# ----------------------------------------------------

did2 <- feols(employed ~ connected + treatment | time, data)
did2 <- summary(did2, cluster = "location")
did2


# Fixed effects DiD estimation (with location fixed effects)
# Location fixed effects absorb the connected variable
# ----------------------------------------------------

# Table 3 HP
# Panel A, column (3)
did3 <- feols(employed ~ treatment | time + location, data)
did3 <- summary(did3, cluster = "location")
did3


# Panel B, column (1)
did4 <- feols(loghours ~ treatment | time + location, data)
did4 <- summary(did4, cluster = "location")
did4


# Event study with interaction between time and treatment
# Quarter x treatment interaction terms, with the pre-treatment period (time = 0) the reference 'category'
# ----------------------------------------------------

did5 <- feols(employed ~ connected::time(0) | time + location, data)
did5 <- summary(did5, cluster = "location")
did5

coefplot(did5)


# View all estimations together -- check that this works in Jupyter
etable(did1, did2, did3, did4, did5, 
    cluster="location", subtitles = c("Standard DiD", "Quarter FE", "Location + quarter FE", "Log hours", "Event study"))



####################################################################
#   PHASE 4: Extensions
####################################################################

# Varying the assumed connection radius 
# Fig 4 in HP
# ----------------------------------------------------

# Define vector of radii
radius_range <- as.matrix(seq(400, 2000, by=100))

# Define function to estimate effect for each radius
radius_vary <- function(r) {          
  
  if (is.numeric(r)) {

    data$newconnected <- ifelse(distance <= r, 1, 0)
    data$newtreatment <- data$newconnected * submarines

    out <- feols(employed ~ newtreatment | time + location, data)
    list(radius=r, effect=out$coeftable[[1,1]], se=out$coeftable[[1,2]])

  } else {
    print("ERROR - Radius must be numeric")
  }

}

# Run function over rows (1) of possible radii vector
radius_out <- apply(radius_range, 1, radius_vary)


# Plot estimates

# convert list to table and create confidence intervals
radius_table <- rbindlist(radius_out)
radius_table$cil <- radius_table$effect - 1.96*radius_table$se
radius_table$cih <- radius_table$effect + 1.96*radius_table$se
radius_table$sig <- ifelse(abs(radius_table$effect)/radius_table$se>1.64, radius_table$effect, NA)


  ggplot(radius_table, aes(x = radius, y = effect)) +
    geom_line() +
    geom_point(aes(x = radius, y = sig), shape = 18, size  = 3) +
    geom_errorbar(aes(ymin  = cil, ymax  = cih), width = 0.2, size  = 0.7) +
    geom_hline(yintercept = 0, linetype="dashed") + 
    theme_bw(base_size = 20) 



# Placebo effect of road and electricity infrastructure (no data on 3G)
# Table 4, column (5) in HP
# ----------------------------------------------------

# Create road and electricity connected variables
data$connected_elec <- ifelse(distance_elec <= radius, 1, 0)
data$connected_road <- ifelse(distance_road <= radius, 1, 0)

# Create placebo treatment variable (treated x post)
data$treatment_elec <- data$connected_elec * submarines
data$treatment_road <- data$connected_road * submarines

# Fixed effects DiD estimation (with location fixed effects)
out <- feols(employed ~ treatment + treatment_elec + treatment_road | time + location, data)
out <- summary(out, cluster = "location")
out



####################################################################
#  Small assignment solution  
####################################################################


# Question 2
# Common trends using an alternative definition of the treatment
# ----------------------------------------------------

# Create treated group variable
data$connected_q2 <- ifelse(data$distance <= 1500, 1, 0)

# Create treatment variable (treated x post)
data$treatment_q2 <- data$connected_q2 * data$submarines


# Create a group-means data set 
common_trends_q2 <- data %>% group_by(time, connected_q2) %>% summarise(empl_mean_q2 = mean(employed))

common_trends_q2$connected_q2 <- factor(common_trends_q2$connected_q2, 
                        levels = c(0,1), label = c("Connected = 0", "Connected = 1")) 

ggplot(data = common_trends_q2, aes(x = time, y = empl_mean_q2, group = connected_q2, color = connected_q2)) + 
  geom_line() +
  geom_vline(xintercept = 0, linetype="dashed") + 
    theme_bw(base_size = 20)

    
# Question 3
# Effect of fast Internet on employment across space 
# Fig 7 in HP
# ----------------------------------------------------

# Define vector of alternative treatments defined by distance to backbone network
data$connected2 <- ifelse(556.6 < distance & distance <= 1500, 1, 0)
data$connected3 <- ifelse(1500 < distance & distance <= 2500, 1, 0)
data$connected4 <- ifelse(2500 < distance & distance <= 3500, 1, 0)

data$treatment2 <- data$connected2 * submarines
data$treatment3 <- data$connected3 * submarines
data$treatment4 <- data$connected4 * submarines


# Run regression with multiple treatments 
out <- feols(employed ~ treatment + treatment2 + treatment3 + treatment4 | time + location, data)
out <- summary(out, cluster = "location")
out


# Plot estimates

space_effects <- as.data.frame(matrix(, nrow=4, ncol=1))
space_effects$radius <- c(1:4)
space_effects$effect <- out$coeftable[,1]
space_effects$se <- out$coeftable[,2]
space_effects$cil <- space_effects$effect - 1.96*space_effects$se
space_effects$cih <- space_effects$effect + 1.96*space_effects$se
space_effects$sig <- ifelse(abs(space_effects$effect)/space_effects$se>1.64, space_effects$effect, NA)
space_effects


  ggplot(space_effects, aes(x = radius, y = effect)) +
    geom_line() +
    scale_x_continuous(name = "Radius", 
                       breaks = c(1,2,3,4), 
                       labels=c("0-500m", "500-1500m", "1500-2500m", "2500-3500m")
                     ) + 
    geom_point(aes(x = radius, y = sig), shape = 18, size  = 3) +
    geom_errorbar(aes(ymin  = cil, ymax  = cih), width = 0.2, size  = 0.7) +
    geom_hline(yintercept = 0, linetype="dashed") + 
    theme_bw(base_size = 20)


