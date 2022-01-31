library(bmlm)
library(brms)
library(qgraph)
library(tidyverse)
library(bayestestR)
options(mc.cores = parallel::detectCores())


df <- BLch9 %>% 
  filter(time < 3)
df <- subset(df, select = -c(fwkstrs, fwkdis, freldis))
#Just to have og model to compare to. Subsequent models should be mirroring these kind of effect sizes
fit <- mlm(d = BLch9,
           id = "id",
           x = "x", #IV
           m = "m", #mediator
           y = "y", #DV
           iter = 5000, 
           cores = 4)
mlm_path_plot(fit, level = .95, text = T,
              digits = 2)
#I need Wide format for simple regression BRMS, I need Long for multilevel stuff in BMLM, or nested BRMS stuff
#lets pretend we have 2-2-1 structure. Only have a varying y1,y2. X and M are same across both time points for each person
df_w <- reshape(df, idvar = "id", timevar = "time", direction = "wide")
df_w <- subset(df_w, select = -c(x.2,m.2))
df_l <- reshape(df_w, idvar = "id", 
                varying = c("y.1", "y.2"),
                v.name = "y",
                direction = "long")
df_l <- df_l[order(df_l$id),]
#Ok, have two equivalent, mock dataframes
#Both as 2-2-1: only y changes across time
#using numeric data

#Now to compare BMLM output to BRMS (both nested and unnested) output

#####BMLM Version. Numeric X
fit_bmlm <- mlm(d = df_l,
             id = "id",
             x = "x.1", #IV
             m = "m.1", #mediator
             y = "y", #DV
             iter = 2000, 
             cores = 4)

fit_bmlm_5k <- mlm(d = df_l,
           id = "id",
           x = "x.1", #IV
           m = "m.1", #mediator
           y = "y", #DV
           iter = 5000, 
           cores = 4)#default is use 4 chains, thus 4chains 5000iters each, we just specify here using STAN code that fits with MLM
mlm_path_plot(fit_bmlm, level = .95, text = T,
              digits = 2)
mlm_path_plot(fit_bmlm_5k, level = .95, text = T,
              digits = 2)
#text = T gives us the top left corner diagnostics
#me = average mediated effect
#c = total effect
#pme = proportion mediated effect
#cov(a,b) = covariance of the varying a and b parameters (i.e. paths a and b)


######BRMS Version
#need wide df here if using basic regression techniques.
#MIGHT be able to use long if trying to do time nesting (but that seems like a lot...)

###Single DV, Time as covariate: Numeric X
fit_brmsA <- brm(m.1 ~ x.1,
                 df_w,
                 warmup = 1000, iter = 2000, chains = 4)

fit_brmsB <- brm(y.2 ~ y.1 + x.1 + m.1,
    df_w,
    warmup = 1000, iter = 2000, chains = 4)
summary(fit_brmsA) #ok this is actually rather similar to the bmlm a path estimates
summary(fit_brmsB) #this is where brms piece-wise looks different from bmlm--x and m look very different

fit_brmsA_5k <- brm(m.1 ~ x.1,
                    df_w,
                    warmup = 1000, iter = 5000, chains = 4)

fit_brmsB_5k <- brm(y.2 ~ y.1 + x.1 + m.1,
                    df_w,
                    warmup = 1000, iter = 5000, chains = 4)
fit_brmsC_5k <- brm(y.2 ~ y.1 + x.1,
                    df_w,
                    warmup = 1000, iter = 5000, chains = 4) #to check the "total" C path
summary(fit_brmsA_5k)#ok, this still looks really consistent with the bmlm estimates
summary(fit_brmsB_5k)#Ok, looks a bit better for matching. Still different, and less precise too, i think largely because of that extra covariate and the call out of an intercept(?). The b path has a switched sign, which is concerning--both to my results and to the tutorial's known b path (~0.15)
#doing this brms would be alright if I knew what was up with that sign error on the b path

##Single DV, Time as Nested: Numeric X
fit_brmsB_MLM <- brm(y ~ x.1 + m.1 + time + (time| id),
                  df_l,
                  warmup = 1000, iter = 5000, chains = 4)
summary(fit_brmsB_MLM) #this actually kinda works tbh
#I only specify the time element to make sure the variance is properly aportioned to the right parameters. I don't care about the effect of time that much otherwise
#otherwise, the b and c effects actually look pretty good compared to the full bmlm
#if this version can do categorical x too, which may make it a winner
#let's make proper dummy variables, dividing up the x variable into 3 categories: 0,1,2
df_l <- df_l %>% 
  mutate(x.c = as.factor(as.numeric(cut(df_l$x.1, 3))-1)) #the -1 is just to make "012", rather than "123" categories

df_w <- df_w %>% 
  mutate(x.c = as.factor(as.numeric(cut(df_w$x.1, 3))-1))

contrasts(df_l$x.c)#check the contrast coding is good
model.matrix(~df_l$x.c)[, 2:3] #one way to isolate the dummy coding for a whole column

df_l <- df_l %>% 
  mutate(xD1 = model.matrix(~df_l$x.c)[, 2]) %>% 
  mutate(xD2 = model.matrix(~df_l$x.c)[, 3]) #inserting those dummy vars as explicit columns

#now, I want to see if the dummy coded brms matches with the dummy coded bmlm--at least for path a
####Categorical X, BRMS

fit_brmsA_5k_C <- brm(m.1 ~ x.c,
                    df_w,
                    warmup = 1000, iter = 5000, chains = 4) #a path with categorical x, this should look like the two bmlm dummy a paths
fit_brmsB_5k_C <- brm(y.2 ~ y.1 + x.c + m.1,
                    df_w,
                    warmup = 1000, iter = 5000, chains = 4) #trying time as covariate again, but with categorical X
fit_brmsB_MLM_C <- brm(y ~ x.c + m.1 + time + (time| id),
                     df_l,
                     warmup = 1000, iter = 5000, chains = 4) #should be more precise version of time as covariate model, but I want to know if it works w/categorical X

summary(fit_brmsA_5k_C)
summary(fit_brmsB_5k_C)
summary(fit_brmsB_MLM_C)
#A* ok, longer chains make this seem like less of a catastrophe, but brms and bmlm are still different enough to be suspicious. I favor the bmlm b/c of its greater precision and functionality (e.g. figures) 
#Option B) Use manual dummy coding in each model to force allowance with BMLM (e.g. 0 vs 1, KM vs other: PA vs other: Neutral vs Other--now need 3 models to isolate each effect I think?)
#Option B)* if using this, I will be limited to just ROPE comparisons I think? Might still be able to do 'KM vs other' model vs 'PA vs Other' model kind of stuff
#B* I think this works if using the 'analysis visualizing' (see pwpt in S1b thesis folder)
#Option C) Use the brm, piece-meal version. Though my initial testing here suggests they're nowhere near as precise, potentially problematically so given the small size of effects
#I'm suspicious of the signage and relative lack of precision in these estimates.


###mvbind-ed DV--yeah, not working in BMLM unfortunately, but not worth pursuing in brms if I don't like the brms process for this

#OK, so now to try with the manual dummy coding. Will need to re-jig the df again to make proper dummies. See 'analysis visualizing' pwpt


####Categorical X, BMLM
#This means I have to run the mediation twice--one for each dummy--
#to get the effects of each dummy var (each manip, all relative to 0 condition var which does not get a param estimate other than intercept)
#The mediator, b path, should be the same for each of the dummy models. The a path and c path should change though.
fitD1 <- mlm(d = df_l,
             id = "id",
             x = "xD1", #IV
             m = "m.1", #mediator
             y = "y", #DV
             iter = 5000, 
             cores = 4)# should be path a part 1, path c part 1
fitD2 <- mlm(d = df_l,
             id = "id",
             x = "xD2", #IV
             m = "m.1", #mediator
             y = "y", #DV
             iter = 5000, 
             cores = 4) #should be path b part 2, path c part 2
#b path for both bmlm dummy models should be the same. The only thing that should really differ is path a (what I mainly care about) and c
mlm_path_plot(fitD1, level = .95, text = T,
              digits = 2)
mlm_path_plot(fitD2, level = .95, text = T,
              digits = 2)
#hm alright, the b path is pretty close. Correct sign and magnitude though, CIs overlap a lot too. But not as neat as I would prefer tbh
#compared to the fit5k, b path is actually pretty decent considering it's all stochastic
#Using this method, I can only use ROPE, not model comparisons (I think) since I won't be able to tell what the full mediation
#model is, since the c path will be partitioned into each dummy effect's remaining direct effect on Y. Technically able to create the two c' effects as per the equation, but won't show as a clean simple SINGLE c' effect
#bmlm is way more precise than brms for mediation (kinda duh I guess--also in part from appropriately tying y.1,y.2 variance to each other)
#but bmlm needs manual dummy modeling to use multicategorical X. 
#i.e. each mediation has 2 dummy models. 2 mediations (cuz two Ms, PA and pKM) = 4 models for each DV
#6 DVs, means 24 models total if using bmlm
#re-running for consistency with different chain lengths, means ~30 models(?)
#a lot, but using Bayes, means model creation is a little less a problem

####mvbind M and/or Y?

#first need to make a correlated y var
#maybe try this with m too--you could treat both pKM and PA mediators similarly actually, both correlated dvs on the a path
df_l$y2 <- -1*(df_l$y + rnorm(length(df_l$y), 0.5,2))
df_l$m2.1 <- 1*(df_l$m.1 + rnorm(length(df_l$m.1), 0.5,2))
with(data = df_l, cor(y, y2)) #decent correlation. Negative (as we'd expect sat or frus) and a decent size r = 0.39
#cor(df_l$y, df_l$y2) same thing as above, just different format
with(data = df_l, cor(m.1, m2.1)) #decent correlation. r=0.5, we want to see positive corr here

#next step is to see how to do the step 1 and step 2 models (with and without correlation) using brms mvbind

options(mc.cores = parallel::detectCores())#I have 8 cores?? cool
set_rescor(FALSE)
fit_m <- brm(mvbind(m.1,m2.1) ~ x.c,
            df_l,
            warmup = 1000, iter = 5000)
summary(fit_m)

fit_y <- brm(mvbind(y,y2) ~ x.c + m.1 + m2.1 + time + (1+time|p|id),
             df_l,
             warmup = 1000, iter = 5000, chains = 4) #check what random correlations are being made here
summary(fit_y) #ok so the correlation goes for each one combination of intercept and time for each DV (6 of them)
#need to find a way to isolate the time and intercept pieces, b/c those make sense paired to themselves but not to each 

fit_y <- brm(mvbind(y,y2) ~ x.c + m.1 + m2.1 + time + (time+0|p|id) + (1|q|id),
             df_l,
             warmup = 1000, iter = 5000, chains = 4)#woo this works for the correlations!!!
#I would then do this for the other two pair of Ys (sat/frus)

#####Ok, I think this works. What does the wrapping do though? bf and calling brms later?
#lets compare the fit_m and fit_y I did, to what happens when we try to put them in a single "mediation" model

m_mod <- bf(
  mvbind(m.1,m2.1) ~ x.c
  )
y_mod <- bf(
  mvbind(y,y2) ~ x.c + m.1 + m2.1 + time + (time+0|p|id) + (1|q|id)
)
my_mod <- brm(m_mod + y_mod + set_rescor(FALSE),
              data = df_l,
              warmup = 1000, iter = 5000, chains = 4) #this doesn't work with mvbind I guess, can't tell the DVs

m_mod <- bf(m.1 ~ x.c)
y_mod <- bf(y ~ x.c + m.1 + m2.1 + time + (1 + time|id))
my_mod <- brm(m_mod + y_mod + set_rescor(FALSE),
              data = df_l,
              warmup = 1000, iter = 5000, chains = 4) #this seems happier. Doesn't like mvbind apparently when combining two together

summary(my_mod)#parameters look remarkably similar to the mvbind-ed fit_m, fit_y estimates--good to see overall

#####(A) Final step, try the main pieces again but with made up priors. Just to see how to properly get the prior notation working

#BRMS has priors for each population effect (i.e. fixed effect)
#default is flat, completely uninformed over the reals


##Broad strokes
#get_prior(fit_m) #doesn't work apparently
prior_summary(fit_m)
fit_m2 <- brm(mvbind(m.1,m2.1) ~ x.c,
             df_l,
             prior = c(
               set_prior("normal(.1,.1)", coef = "x.c1", resp = c("m1", "m21")),
               set_prior("normal(0,5)", class = "Intercept", resp = c("m1", "m21"))
             ),
             warmup = 1000, iter = 5000)

#ok notice, that in the og fit_m output, the effects have weird names b/c we have 2 dvs mvbinded.
#ok notice, b/c we have categorical x here, we have dummy vars with subsequent 'dummy effects'
#the priors are finnicky to specify, e.g. notice how brms removes "." in the resp var so now it thinks the resp var is m1 and m21 (not m2.1)
prior_summary(fit_m2)
summary(fit_m)
summary(fit_m2)
#Whew this works with mvbind! The priors noticeably are affecting the estimates of the correct parameters I specified
#And I can see the priors being correctly specified in the prior_summary: I can specify priors in mvbind then!


#fit_y <- brm(mvbind(y,y2) ~ x.c + m.1 + m2.1 + time + (time+0|p|id) + (1|q|id),
#             df_l,
#             warmup = 1000, iter = 5000, chains = 4)

##Finer comb


m_mod <- bf(m.1 ~ x.c)
y_mod <- bf(y ~ x.c + m.1 + m2.1 + time + (1 + time|id))
my_mod <- brm(m_mod + y_mod + set_rescor(FALSE),
              data = df_l,
              prior = c(
                set_prior("normal(.1,.1)", coef = "x.c1", resp = "y"),
                set_prior("normal(0,5)", class = "Intercept", resp = "y")
              ),
              warmup = 1000, iter = 5000)
prior_summary(my_mod) #ok prior setting works this way too it seems

###Does mediation package work?
m_mod <- bf(m.1 ~ x.c)
y_mod <- bf(y ~ x.c + m.1 + m2.1 + time + (1 + time|id))
my_mod <- brm(m_mod + y_mod + set_rescor(FALSE),
              data = df_l,
              warmup = 1000, iter = 5000)
mediation(my_mod) #eh, doesn't seem to get the multicategorical X part down, just shows
#treatment as x.c2, not x.c1 (the other category).
summary(my_mod)
#Also the reported indirect effect of x.c2 is almost exactly the same as the x.c2 effect in the model summary
#i.e. I don't think the mediation numbers are showing me anything notably different from what the model summary would tell me

###Ok, new comparison. (1|id) vs (1+time|id)
#considering it's just t1,t2, this may not actually matter since the random intercept is effectively capturing t1 and t2 is the dv anyway

fit1 <- brm(y ~ x.1 + m.1 + time + (time| id),
    df_l,
    warmup = 1000, iter = 8000, chains = 4)

fit2 <- brm(y ~ x.1 + m.1 + time + (1 | id),
            df_l,
            warmup = 1000, iter = 8000, chains = 4)
summary(fit1)
summary(fit2)#very slight difference, just on the m effect (to the .01 decimal change in CI and estimate)
#so overall, probs not a big difference with this kind of t1, t2 data I'd say
#for overidentification sake, may thik about just using random intercept
#trying with categorica X just in case
fit1 <- brm(y ~ x.c + m.1 + time + (time| id),
            df_l,
            warmup = 1000, iter = 8000, chains = 4)

fit2 <- brm(y ~ x.c + m.1 + time + (1 | id),
            df_l,
            warmup = 1000, iter = 8000, chains = 4)
summary(fit1)
summary(fit2)
#as expected, the a effects get way wider with less variation from the categorical x
#m effects stay about the same, both between these two and compared to the numeric x versions


###New thing from Adam--use t1 measures as random slopes? ie (1 + t1|ID)
#this should use df_w I think...

#something like this?
fit1 <- brm(y.2 ~ x.c + m.1 + (1+y.1|id),
             df_w,
             warmup = 1000, iter = 5000, chains = 4)
#compare to? I think this, this is the 'final' model that I would probs use, without the mvbind
fit2 <- brm(y ~ x.c + m.1 + time + (1|id),
             df_l,
             warmup = 1000, iter = 5000, chains = 4)

summary(fit1)
summary(fit2)
##Looks to me like the former has a type s error in the m path. I also think the latter model
##uses time/phase/whatever more consistently to the way i've seen longitudinal models get run in the past online and in course
##I trust that latter model construction better


##also, try out the difference between m ~ x, m ~x + (1|ID)
fit1 <- brm(mvbind(m.1,m2.1) ~ x.c,
             df_l,
             warmup = 1000, iter = 5000)

fit2 <- brm(mvbind(m.1,m2.1) ~ x.c + (1|id),
             df_l,
             warmup = 1000, iter = 5000)#having a lot of trouble with this sampling--taking 30mins per chain
summary(fit1)
summary(fit2)
##really hard time fitting the 2nd model here (1|id) with a  lot of sampling warnings
##I'd be careful with it. 


###then one last confirmation of whether the bf() + bf(), in brms() works with mvbind or not
###initial tests suggest brm isn't happy with mvbind within a bf
fit_m <- bf(mvbind(m.1, m2.1) ~ x.1)
fit_y <- bf(mvbind(y,y2) ~ x.c + m.1 + m2.1 + time + (time|id))
model_my <- brm(fit_m + fit_y, 
                 set_rescor(FALSE),
                 data = df_l,
                resp = c("m.1", "m2.1", "y", "y2"))
#yeah, this isn't super happy with using mvbind twice--it gets confused what is the response variable
#and I can't seem to find out how to specify the resp= part properly
#I don't think that's worth it