

## Launch a R session
## Choose the appropriate working directory avec setwd("patway/to/directory")

###########################
# Packages and workspace  #
###########################

# import the survival package
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)


############################
# create accdata dataframe #
############################

accdata <- read.table("6-GSE10927-norm.txt", header=TRUE, sep="\t", dec=".", row.names=1)

# The following commands select lines with known survival data and indicate that these data are numbers.
survidata <- subset(accdata, accdata$time != "unknown")
survidata$time <- as.numeric(as.character(survidata$time))
survidata$status <- as.numeric(as.character(survidata$status))
survidata$clusterforACCs <- factor(survidata$clusterforACCs)
survidata$MitoticRate <- factor(survidata$MitoticRate)
survidata$MitoticRate <- as.numeric(as.character(survidata$MitoticRate))


#########################
# Kaplan Meier analysis #
#########################

### Build the standard "survival" object (which links time and events)

km <- with(survidata, Surv(time, status))
# km contains periods (time to last follow-up) and the status (+) indicates censored data.


### Draw the Kaplan Meier plot

jpeg("Kaplan-Meier_curves.jpg", width = 700, height = 700)
km_fit <- survfit(Surv(time, status) ~ survidata$clusterforACCs, data=survidata)
autoplot(km_fit)
dev.off()


### Log-rank test for the difference between the two groups

survdiff(Surv(time, status) ~ survidata$clusterforACCs, data=survidata)

