###BRMS mvbind() practice
#Outside of brms environment, mvbind() works just like cbind()
#Basically mvbind is for multivariate preds

#mvbind((y1, y2) ~ x)
#> y1 ~ x
#> y2 ~ x

library(brms)
data("BTdata", package = "MCMCglmm")
head(BTdata)

#Simple multivariate normal model
#note brm() and brms() are both within the brms library
#brm is just for NON-linear models (I think)
fit1 <- brm(
  mvbind(tarsus, back) ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam),
  data = BTdata, chains = 2
)
#the extra |p| (or q) in the grouping part of the formula is apparently to indicate varying effects of 
#fosternest should be modeled as correlated. I.e. the random intercept of fosternest for y1 should be correlated to
#the random intercept of fosternest on y2 (i think)
#Though if the case, I would expect you would have to do that for every usage of mvbind?
#key point is that |<ID>|, <ID> can be anything p, q whatever. The important bit is having the double bar ||
# from the vignette CRAN on advanced BMRS usage (page 4), this || for correlated grouping seems to be an important
#feature of non-linear models or distributional models (idk what those latter are) since they have multiple parameters
#each with their own population(fixed) and group(random) level effects, therefore you need multiple formulas, therefore you need to 
#run into the problem of specifying those group-level effects of the same grouping var are correlated ACROSS those formulas is tough (the || is BRMS' workaround of that issue)

#Group level terms with the same <ID> in the || are modeled as correlated if they share the same grouping factor
#with a multi-level model of mvbind, that's a given since its basically running the model (same grouping vars) twice with different DVs
#TLDR: (1| group) is standard as you'd expect, one formula
#(1|ID|group) means in two formulas, models the varying intercepts of group across those two model formulations as correlated
#so yeah, I think we'd use that for pretty much every mvbind situation, like in the tutorial there

fit1 <- add_criterion(fit1, "loo")
#loo, means Leave-one-out cross validation. It's a way to evaluate the performance of a ml algorithm
#it's functionally the logical extreme of K-fold cross validation
#cross validation as a term refers to a model evaluation method that is better than residual evaluations
#residual based evaluation is entirely within data (about the errors/residuals from your model, in your data)
#crossvalidation is about the test/train splits to evaluate the model based on performance to "new" data
#the adding criterion here is so that we can compare it to another model later if we want
#e.g. loo(fit1, fit2)
summary(fit1)

#I would use (1 + time |p| ID), which would explicitly model that the intercept of sat and frus
#within an individual are correlated, which is what I want.
#it would also give me the time effect which i don't think really matters that much, unless I care about the 
#plasticity of sat or frus over time generally (which ain't the point of this study)

#you also get a nice residual correlation at the bottom
#this tells you how much unmodeled dependency still exists between your y1 and y2
#if magnitude is big, that means residuals of y1 and residuals of y2 are correlated
#i.e. there is outstanding dependency (the residuals correlation indicating some relationship) between y1 and y2 that we aren't accounting for

#now we can start looking at model fit
#posterior-predictive checks
library(bayesplot)
pp_check(fit1, resp = "tarsus") #compares observed y, to simulated datasets yrep from posterior predictive distribution
pp_check(fit1, resp = "back")

library(rstantools) #or in brms too
bayes_R2(fit1)

#will want to check if mvbind is literally pulling to lms together, or whether it does make some use of correlating the DVs together
#if just two lms at once, will defos need the |<ID>| thing to specify IV correlations
#mvbind((y1, y2) ~ x)
#> y1 ~ x
#> y2 ~ x

fit2 <- brm(
  tarsus ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam),
  data = BTdata, chains = 2, cores = 2)

fit3 <- brm(
  back ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam),
  data = BTdata, chains = 2, cores = 2)

fit2 <- add_criterion(fit2, "loo")
fit3 <- add_criterion(fit3, "loo")
loo(fit1, fit2, fit3)
summary(fit1)
summary(fit2)
summary(fit3)
#just comparing fixed population effects estimates AND CIs mvbind(y1,y2) _ x literally does just mean running two lms, y1 ~ x, y2 ~x
#also in terms of "family specific" parameters, group-level effects
bayes_R2(fit1)
bayes_R2(fit2)
bayes_R2(fit3)
#R2 also all the same too

#loo can't really help me here to tell which is better since not all the models have the same y variable
#So yes, as far as I can see, mvbind is literally JUST running two lms at once. Does give extra correlation stats baked in, though they don't seem to affect parameter estimation at all

###What about comparing the use of |<ID>| or not?
fit1 <- brm(
  mvbind(tarsus, back) ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam),
  data = BTdata, chains = 2, cores = 2,
  save_pars = save_pars(all = TRUE))
fit4 <- brm(
  mvbind(tarsus, back) ~ sex + hatchdate + (1|fosternest) + (1|dam),
  data = BTdata, chains = 2, cores = 2,
  save_pars = save_pars(all = TRUE)) #same as fit1 model, but without the random group intercepts marked as correlated
#fit1 <- add_criterion(fit1, "loo")
#fit4 <- add_criterion(fit4, "loo")
#may not actually need add loo criterion here, perhaps b/c of saving all parameters in the model itself?
summary(fit1)
summary(fit4) 
#random effect parameter estimates, est errors, and CIs all look the same between fit1(correlations) and fit4
#population and family specific parameters are all basically identical. MAYBE correlated version is tiny bit more precise, but by incredibly small margins and likely not nearly as large of precision that would come from simply more data/priors

library(bayestestR)
loo(fit1, fit4)#according to this, I THINK, fit1 is marginally better than fit4 but not by much at all
loo(fit4, fit1)#confirming, it's not just the ordering of models--fit4 does keep showing as -7.4
bayesfactor_models(fit1, fit4)
bayesfactor_models(fit4, fit1)

####
####
####
#So we learned:
#mvbind just combines two lms into one. Does nothing else fancy or change the estimates of the models themselves in any meaningful way
#|<ID>| seems to be marginally better than not specifying the correlation via tiiiiiiny bit more precision and better loo, but not by a lot (at least in this data)
#bayesfactor is messy, probs cuz these are really small model samplings



###Trying the prios in mvbind
brm(
  mvbind(tarsus, back) ~ sex + hatchdate + (1|p|fosternest) + (1|q|dam),
  data = BTdata, chains = 2,
  prior = c(
    set_prior("normal(0, 5)", class = "b", resp = c("tarsus", "back"),
              coef = "sexMale"),
    set_prior("normal(0, 5)", class = "b", resp = c("tarsus", "back"),
              coef = "sexUNK"),
    set_prior("normal(0, 5)", class = "b", resp = c("tarsus", "back"),
              coef = "hatchdate"),
    set_prior("normal(0, 5)", class = "Intercept", resp = c("tarsus", "back"),
              coef = "")),
  sample_prior = "only")

prior_summary