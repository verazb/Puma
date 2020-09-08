### AEA2020 ###
# - Session 2 -#


####################################################################
#   LOAD PACKAGES AND READ IN DATA  
####################################################################

# remove old objects from R working space
rm(list=ls())

# Define packages that you need
packages_vector <- c("tidyverse", "haven", "Hmisc", "dplyr", 
                     "tidyr", "stringr", "sandwich", "lmtest", 
                     "jtools", "fBasics", "knitr", "xtable", 
                     "data.table", "stargazer", "AER", 
                     "causalweight", "np")

# install.packages(packages_vector)
lapply(packages_vector, require, character.only = TRUE) 
# List loaded packages 
(.packages())

# LOAD FUNCTIONS
  source("Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/First Session/radiusmatch.R")
  source("Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/First Session/radiusatet.R")
  source("Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/First Session/inferenceweights.R")
  source("Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/First Session/inferenceweights2.R")


# Working Directory
work_dir <- "Q:/Arbeitsmarktökonomie/Data-Lehre/Methods/2020 Applied Empirical Analysis/Applications/First Session"
setwd(work_dir)

# Read in Data
load("data_reg_final.RData")
attach(data_reg)

##################################
# Define variables
##################################

# Group variables, adding newly created ones - need matrix!
y1 <- as.matrix(duration)
y2 <- as.matrix(employed1y)
treat <- as.matrix(treat)
x1 <- as.matrix(dplyr::select(data_reg, xcat_names_reg, xcont_names, starts_with("quarter_")))
x_names<- colnames(x1)

####################################################################
#   Semiparametric Estimation
####################################################################


##################################
# P-score  
##################################


# 1) Estimate the treatment probability
pscore.model <- glm(treat ~ x1, family = binomial(link = "probit"))
summ(pscore.model, , robust = "HC1")

# 2) Estimate the propensity score
data_reg$pscore <- pscore.model$fitted.values 
data_reg$treat_f <- factor(data_reg$treat, levels = c(0,1), label = c("D=0", "D=1")) 


# 3) Check for common support in propensity score and trim
ggplot(data_reg, aes(x = pscore, fill = treat_f)) + geom_density(alpha=0.4) + scale_fill_grey() + theme_classic() + xlim(0, 1)


# Trim: delete treated observations with a p-score larger as the maximum p-score in the control population
nrow(data_reg)
pscore_max0 <-max(data_reg$pscore[treat==0])
insample<- ifelse(data_reg$pscore <= pscore_max0, 1, 0)
data_reg <- dplyr::filter(data_reg, insample ==1)
nrow(data_reg)

# refine variables
y1 <- y1[insample==1]
y2 <- y2[insample==1]
treat <- treat[insample==1]
x1 <- x1[insample==1,]


####################################################################
#  IPW
####################################################################

  # using the causalweight package by Huber et al.
  ipw_atet <- treatweight(y=y1, d=treat, x=x1, ATET =TRUE, trim=0.05, boot = 199)
  ipw_atet

  # Manually
  # generate weights based on estimated pscores
  data_reg$weight[treat==1] <- 1
  data_reg$weight[treat==0] <- data_reg$pscore[treat==0]/(1-data_reg$pscore)[treat==0]
  weight <- data_reg$weight

  Y1<-sum(y1*treat)/sum(treat)
  Y0<-sum(y1*(1-treat)*weight)/sum((1-treat)*weight)
  ipw_atet_manual <- Y1-Y0
  ipw_atet_manual

# You could impose some additional trimming large values of the pscore by conditioning on some indicator:
# ind=as.numeric((data_reg$pscore<0.05) | (data_reg$pscore>(0.95))) 
# Y1<-sum(y1[ind==0]*treat[ind==0])/sum(treat[ind==0])
# Y0<-sum(y1[ind==0]*(1-treat[ind==0])*weight[ind==0])/sum((1-treat[ind==0])*weight[ind==0])

####################################################################
# Check balancedness in covariates
####################################################################


balanced.model <- function(x){
                        weighted_diff <- lm(x[insample==1] ~ treat,weight=data_reg$weight )
                        cov <- vcovHC(weighted_diff, type = "HC")
                        robust.se <- sqrt(diag(cov))
                        list( coeff = weighted_diff$coefficients[2], robust.se = robust.se[2], pval = 2*pnorm(-abs(weighted_diff$coefficients[2]/robust.se[2])) )
                        }
diff_output <- apply(x1,2,balanced.model)

# convert list to table
diff_output<-rbindlist(diff_output)
rownames(diff_output)<- x_names

# plot table
print("Weighted differences in X between control and treatment group")
xtable(diff_output, digits=3)


####################################################################
#  Radiusmatching with bias correction
####################################################################

  # using the radiusmatching package by Huber et al.
  rmatch_atet <- radiusmatch(y=data_reg$duration, d=data_reg$treat, x=x1, biascorr=1, commonsup = 1, radius = 3, estimand ="ATET", boot=2 )
  rmatch_atet$effect
  rmatch_atet$se.boot

  # radius is defined by the multiplier of the maximum distance in pair matching  E.g., setting radius=3 implies a radius of 3 times the maximum 
  # distance (or a particular quantile).

###################################################
# Dynamic outcomes
###################################################


# # Exit rates into employment

# exit <- matrix(0, nrow=nrow(data_reg), ncol=(maxdur))
# exit_list <- paste("exit", 1:maxdur, sep="_")
# colnames(exit) <- exit_list

# for (i in 1:maxdur) {
# 	exit[,i] = ifelse(data_reg$date_end < data_reg$date_start + 30*i & data_reg$date_end >= data_reg$date_start + 30*(i-1), 1, 0)
# }

# exit <- as.data.frame(exit)
# data_reg <- cbind(data_reg, exit)


# # Graph

# exit_monthly_treat <- apply(X=t(exit[data_reg$treat==1,]), MARGIN=1, FUN=mean)
# exit_monthly_contr <- apply(t(exit[data_reg$treat==0,]),1,mean)

# month <- rep(1:maxdur)
# exit_monthly <- cbind(exit_monthly_treat, exit_monthly_contr, month)
# colnames(exit_monthly) <- c("rate_treat", "rate_contr", "month")
# exit_monthly <- as.data.frame(exit_monthly)


# ggplot() +
#   geom_line(data = exit_monthly, aes(x = month, y = rate_treat), color = "blue") +
#   geom_line(data = exit_monthly, aes(x = month, y = rate_contr), color = "red") +
#   theme_classic() +
#   theme(axis.title = element_text(face = "bold"),
#         axis.title.x = element_blank()) +
#   ylab("Exit rate") +
#   xlab("Months") +
#   # ylim(0, 1) + 
#   geom_hline(yintercept = 0, linetype="dotted") 


    # Set maximum horizon to look at after start of unemployment, = 2 years
    maxdur <- 24

    # Monthly employment status
    empl <- matrix(0, nrow=nrow(data_reg), ncol=(maxdur))
    empl_list <- paste("empl", 1:maxdur, sep="_")
    colnames(empl) <- empl_list

    for (i in 1:maxdur) {
      empl[,i] = ifelse(data_reg$date_end < data_reg$date_start + 30*i, 1, 0)
    }

    empl <- as.data.frame(empl, data_reg$treat)
    data_reg <- cbind(data_reg, empl)



    # Graph

    empl_monthly_treat <- apply(t(empl[treat==1,]),1,mean)
    empl_monthly_contr <- apply(t(empl[treat==0,]),1,mean)

    month <- rep(1:maxdur)
    empl_monthly <- cbind(empl_monthly_treat, empl_monthly_contr, month)
    colnames(empl_monthly) <- c("rate_treat", "rate_contr", "month")
    empl_monthly <- as.data.frame(empl_monthly)


    ggplot() +
      geom_line(data = empl_monthly, aes(x = month, y = rate_treat), color = "blue") +
      geom_line(data = empl_monthly, aes(x = month, y = rate_contr), color = "red") +
      theme_classic() +
      theme(axis.title = element_text(face = "bold"),
            axis.title.x = element_blank()) +
      ylab("Employment rate") +
      xlab("Months") +
      # ylim(0, 1) + 
      geom_hline(yintercept = 0, linetype="dotted") 



      ############################################
      # Monthly outcomes regressions (parametric)
      ############################################


      reg_monthly <- function(y) {
      	out <- lm(y ~ treat + x1)
      	# out2 <- summary(out)
      	out2 <- coeftest(out, vcov = vcovHC(out))
      	list(effect=out2[2,1], se=out2[2,2])
      }

      # Graphs with confidence intervals and dots for significance

      reg_monthly_out <- apply(empl, 2, reg_monthly)
      reg_monthly_out2 <- rbindlist(reg_monthly_out)
      reg_monthly_out2$month <- rep(1:maxdur)
      reg_monthly_out2$cil <- reg_monthly_out2$effect - 1.96*reg_monthly_out2$se
      reg_monthly_out2$cih <- reg_monthly_out2$effect + 1.96*reg_monthly_out2$se
      reg_monthly_out2$sig <- ifelse(abs(reg_monthly_out2$effect)/reg_monthly_out2$se>1.64, reg_monthly_out2$effect, NA)


      ggplot(reg_monthly_out2,                
             aes(x = month, y = effect)) +
        geom_line() +
        geom_point(aes(x = month, y = sig), shape = 18,
                   size  = 3) +
        geom_errorbar(aes(ymin  = cil,
                          ymax  = cih),
                      width = 0.2,
                      size  = 0.7) +
        theme_classic() +
        theme(axis.title = element_text(face = "bold"),
              axis.title.x = element_blank()) +
        ylab("Effect on employment rate (OLS)") +
        # ylim(0, 1) + 
        geom_hline(yintercept = 0, linetype="dotted") 



####################################################################
 # Small Assignment (first part)
####################################################################
####################################################################
####################################################################
####################################################################
####################################################################

  # Danymic effects with IPW instead of OLS

  reg_monthly_a <- function(y) {
  ipw_atet <- treatweight(y=y, d=treat, x=x1, ATET =TRUE, trim=0.05, boot = 19)
  list(effect=ipw_atet$effect, se=ipw_atet$se)
}

# Graphs with confidence intervals and dots for significance

reg_monthly_out <- apply(empl, 2, reg_monthly_a)
reg_monthly_out3 <- rbindlist(reg_monthly_out)
reg_monthly_out3$month <- rep(1:maxdur)
reg_monthly_out3$cil <- reg_monthly_out3$effect - 1.96*reg_monthly_out3$se
reg_monthly_out3$cih <- reg_monthly_out3$effect + 1.96*reg_monthly_out3$se
reg_monthly_out3$sig <- ifelse(abs(reg_monthly_out3$effect)/reg_monthly_out3$se>1.64, reg_monthly_out3$effect, NA)


ggplot(reg_monthly_out3,                
       aes(x = month, y = effect)) +
  geom_line() +
  geom_point(aes(x = month, y = sig), shape = 18,
             size  = 3) +
  geom_errorbar(aes(ymin  = cil,
                    ymax  = cih),
                width = 0.2,
                size  = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"),
        axis.title.x = element_blank()) +
  ylab("Effect on employment rate (IPW)") +
  # ylim(0, 1) + 
  geom_hline(yintercept = 0, linetype="dotted") 


