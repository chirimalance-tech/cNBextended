# SHINY APP
Along with the development of a R package A R shiny app was developed to explore the behaviours of the cNB-RM.

The app is avaliable here   [![Launch Shiny App](https://img.shields.io/badge/Launch%20App-Click%20Here-brightgreen)](https://chirimalance-tech.shinyapps.io/cNBextendedApp/)


# cNBextended
Useful resources developed in R that can be used to understand and implement an extesion of the Contaminated Negative Binomial Regression Model (cNB-RM). 


# Example
Using the package.
Code to reproduce Table 4(b) "Number of visits to the physicians office"

```
To install the package, use the following code in R
install.packages("devtools")
library(devtools)
install_github("chirimalance-tech/cNBextended")

# Alternatively
install.packages("remotes")
library(remotes)
install_github("chirimalance-tech/cNBextended")

library(AER)
library(cNB)
library(cNBextended)

# Some information about the mle function and how it works.
?ml.ex.cnb

# Load NMES1988 dataset
data("NMES1988")
head(NMES1988)

# Construct design matrix
x <- model.matrix(~ visits + nvisits + ovisits + novisits + emergency + hospital +
                    health + chronic + adl + region + afam + gender + married +
                    school + income + employed + insurance + medicaid,
                  data = NMES1988)
x <- as.data.frame(x)

options(contrasts = c("contr.treatment", "contr.poly"))


# Regression on mean, dispersion, delta, and eta
cnb_ex_trips_mu_alpha_delta_eta <- ml.ex.cnb(formula = visits ~ hospital + insuranceyes +
                                                          healthpoor + chronic + adllimited +
                                                          gendermale + school + afamyes,
                                            alpha.formula = ~ adllimited,
                                            delta.formula = ~ nvisits,
                                            eta.formula = ~ chronic,
                                            data = x,
                                            method = "BFGS")

# Parameter estimates
round(cnb_ex_trips_mu_alpha_delta_eta$results, 4)

# Information criteria and log-likelihood
round(cnb_ex_trips_mu_alpha_delta_eta$AIC, 4)
round(cnb_ex_trips_mu_alpha_delta_eta$BIC, 4)
round(cnb_ex_trips_mu_alpha_delta_eta$HQIC, 4)
round(cnb_ex_trips_mu_alpha_delta_eta$loglike, 4)
```
