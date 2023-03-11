library(tidyverse)
library(sjlabelled)
library(huxtable)
library(pubh)
library(jtools)
library(sandwich)
library(msm)

### load data
data(diet, package = "Epi")

### explore
diet <- diet %>%
   var_labels(
    chd = "Coronary Heart Disease",
    job = "Type of job"
  )

diet %>% group_by(job, chd) %>%
  count() %>%
  as_hux() %>% theme_pubh(1)

### fit model 
model_rr <- glm(chd ~ job, data = diet, family = poisson(link = "log"))

### Use the summ function from the jtools package to get Robust SE. 
### You can specify the different types of Robust SE: 
### HC0 - Original formula. No small sample bias adjustment.
### HC1 - Used by Stata, applies a degrees of freedom-based correction, (n − 1)/(n − k)
###       where n is the number of observations and k is the number of explanatory or predictor 
###       variables in the model. Most used, less effective than HC2 & HC3 when the number of 
###       clusters is small.
### HC2 - Designed to be used with GLMs.
### HC3 - Also for GLMS, better for small sample sizes and outperforms all other options in the 
###       presence of heteroskedacity.
### HC4 - Better than HC3 if there are outliers
### HC4m - A different HC4 specification to correct for outliers.
### HC5 - Also, corrected for presence of influential outliers.
### You can also specify clustered-robust SE by providing the clustering variable in the 
### cluster option (default set to NULL)

model_rr %>%
  summ(confint = T, model.info = F, exp = T, robust = "HC3", cluster = NULL)

### Use function glm_coef from package pubh to obtain Robust SE. 
model_rr %>%
  glm_coef(labels = c("Constant", "Conductor", "Bank worker"), se_rob = T) %>%
  as_hux() %>% 
  set_align(everywhere, 2:3, "right") %>%
  theme_pubh(1) %>%
  add_footnote(get_r2(model_rr), font_size = 9)


### Obtain Robust SE using the sandwich package
cov.m1 <- vcovHC(model_rr, type="HC0")

### Another way of doing it is by hand code it. It gives you slighlty different results, 
### because the type of Robust SE the following approach gives you is HC0 robust SE.

### Use the covariance matrix from before to calculate SE's, p-values, and CI's on the log scale
std.err <- sqrt(diag(cov.m1))
r.est <- cbind(Estimate= coef(model_rr), "Robust SE" = std.err,
               "Pr(>|z|)" = 2 * pnorm(abs(coef(model_rr)/std.err), lower.tail=FALSE),
               LL = coef(model_rr) - 1.96 * std.err,
               UL = coef(model_rr) + 1.96 * std.err)
r.est

### Using delta method to compute the SE's of the exponentiated RR
s <- deltamethod(list(~ exp(x1), ~ exp(x2), ~ exp(x3)), 
                 coef(model_rr), cov.m1)

## exponentiate old estimates dropping the p values
rexp.est <- exp(r.est[, -3])
## replace SEs with estimates for exponentiated coefficients
rexp.est[, "Robust SE"] <- s

rexp.est
