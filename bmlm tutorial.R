library(bmlm)
library(brms)
library(qgraph)

head(BLch9)
#MLM means isolating within and between person components of X,M,Y. BLch9 already
#has subject-mean deviated (within-person) components (XMY) but if you had to obtain them from the first three columns for each person use isolate

BLch9_2 <- isolate(BLch9,
        by = "id",
        value = c("fwkstrs", "fwkdis", "freldis"))
head(BLch9_2)
#the _cw ("centered within") now are within-person ("subject mean deviated"). 
#Notice they're the same as the XMY vars (again, we already had it there, this is just to show how to use the function)
#Notice also that X is the same for all of fwkstrs ==3 (that's just showing that X is tied to that column, m, and y respectively for the next two)
#_cb is used to label centered-between

fit <- mlm(d = BLch9_2,
           id = "id",
           x = "fwkstrs_cw", #IV
           m = "fwkdis_cw", #mediator
           y = "freldis_cw", #DV
           iter = 2000, 
           cores = 4)#default is use 4 chains, uses first HALF of iteration as warmup, thus 4chains 2000iters each 1000 warmup, we just specify here using STAN code that fits with MLM
#notice, no use of priors in this one yet
#you'd include that using STAN methods (need to look into this, use Adam's code as base)

mlm_summary(fit) #yay it works!

mlm_path_plot(fit, level = .95, text = T,
              digits = 2)
#text = T gives us the top left corner diagnostics
#me = average mediated effect
#c = total effect
#pme = proportion mediated effect
#cov(a,b) = covariance of the varying a and b parameters (i.e. paths a and b)

mlm_path_plot(fit, level = .95, text = T,
              digits = 2,
              random = F) #turns random effects (i.e. uninteresting, but present, errors) off. Notice the SDs are removed in this version

###Parameter plotting
#mlm_pars_plot(type = X), where X can be a hist, coef, or violin plot
mlm_pars_plot(fit, pars = "me") #default is hist
mlm_pars_plot(fit, type = "hist", pars = c("me", "c", "pme", "covab"))

mlm_pars_plot(fit, type = "coef", level = .99) #coefficient plot with point estimates
mlm_pars_plot(fit, type = "violin", color = "dodgerblue2")

mlm_pars_plot(fit, type = "hist", pars = c(
  "tau_cp", "Omega[1,2]", "Omega[1,3]",
  "Omega[2,1]", "tau_b", "Omega[2,3]",
  "Omega[3,1]", "Omega[3,2]", "tau_a"),
  nrow = 3, color = "skyblue4")
#the varying effects correlations (omegas) and standard devs (tau)
#histograms are mcmc samples of plausible paramter values from the corresponding posterior

mlm_pars_plot(fit, pars = c("tau_a", "tau_b", "tau_cp",
                            "corrab"), nrow=2, color="dodgerblue2")

#traceplots for convergence
library(bayesplot)
mcmc_trace(as.data.frame(fit), pars = c("a", "b", "cp", "corrab"))
#corrab is correlation of a and b (btw)

###Plotting fitted values
#use mlm_spaghetti_plot() for mapping predicted values
ps <- mlm_spaghetti_plot(mod = fit, d = BLch9_2,
                         x = "fwkstrs_cw",
                         m = "fwkdis_cw",
                         y = "freldis_cw")

#ps now has two plots: one for x->m, one for m->y
ps[1]
#population fitted values are 95% credible intervalsshaded interval around a thick line
#subject deviations are in lighter lines

ps[2]



######Trying for myself
My1 <- mlm(d = df,
           id = "id",
           x = "x", #IV
           m = "m", #mediator
           y = "y1", #DV
           iter = 2000, 
           cores = 4)
mlm_summary(My1)
mcmc_trace(as.data.frame(My1), pars = c("a", "b", "cp", "corrab"))

My2 <- mlm(d = df,
           id = "id",
           x = "x", #IV
           m = "m", #mediator
           y = "y2", #DV
           iter = 2000, 
           cores = 4)
My12 <- mlm(d = df,
            id = "id",
            x = "x", #IV
            m = "m", #mediator
            y = cbind("y1","y2"), #Doesn't like the cbind or mvbind in bmlm at this point
            iter = 2000, 
            cores = 4)
