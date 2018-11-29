# Load and install necessary packages
requiredPackages <- c("igraph", "ggplot2", "Matrix", "ggthemes", "gridExtra")

for (pac in requiredPackages) {
    if(!require(pac,  character.only=TRUE)){
        install.packages(pac, repos="http://cran.rstudio.com")
        library(pac,  character.only=TRUE)
    } 
}
rm(pac)
rm(requiredPackages)


# Set WD and load data
wd = getwd()
if(grepl("nora", wd)) {
    setwd("~/Documents/18-19/CSN/LABS/05/src")
} else {
    # Set Piero Working directory
    setwd("~/Uni/csn/CSN-Lab-05/src")
}
rm(wd)

#################
####  Flags  ####
#################

APPLY_PROFILING = FALSE

LOAD_EXISTING_RUN_BA = TRUE
LOAD_EXISTING_RUN_BA_RAND = TRUE
LOAD_EXISTING_RUN_BA_NO_GROWTH = TRUE



TESTING = TRUE

#################
### Variables ###
#################

t.max = 10000
timestamps <- c(1, 10, 100, 1000) # 4 vertices to track

#################
### Sourcing  ###
#################
source("commonFunctions.R")
source("01.BA_original.R")
source("02.BA_growth_rand_att.R")
source("03.BA_no_growth.R")

m.0 = 2
if(TESTING){
   theoretic_ki = c()
   exact_scaled_therotic_ki = c()
   # approximated vertex degree over time (slide 4)
   approx_scaled_therotic_ki = rep(m.0*sqrt(t.max), t.max+1)
   m.0 = 2
   for (t in seq(t.max+1)) {
      # vertex degree over time (slide 4)
      theoretic_ki = append(theoretic_ki, (m.0*((t.max/t)**0.5)))
      # scaled vertex degree over time (slide 4)
      exact_scaled_therotic_ki = append(exact_scaled_therotic_ki, 
                                        sqrt(t)*theoretic_ki[t])
      
   }
}

#plot(seq(t.max+1), exact_scaled_therotic_ki)
#plot(seq(t.max+1), approx_scaled_therotic_ki)





#df.t.k = data.frame(table.BA$sequ, theoretic_ki)
# ggplot(data = table.BA) +
#    aes(x = table.BA$sequ, y = table.BA$ks.4) +
#    geom_point(color = '#ef562d') +
#    geom_hline(yintercept = df.t.k[1003,2]) +
#    theme_minimal()


# model selection ---------------------------------------------------------

# sequence of values
t = seq(t.max+1)

# models without intercept
model0 = nls(theoretic_ki~a*t,start=list(a=1))
model1 = nls(theoretic_ki~a*sqrt(t),start=list(a=1))
model2 = nls(theoretic_ki~a*(t)^b, start = list(a=0.1, b=0.1), control = nls.control(tol = 1e-1 ))
model3 = nls(theoretic_ki~a*exp(c*t), start = list(a=1,c=0.0001))
model4 = nls(theoretic_ki~a*(log(abs(t+d1))), start=list(a=0.1 ,d1=1), control = nls.control(tol = 0.7))

# models with intercept
model0i = nls(theoretic_ki~a*t+d,start=list(a=1, d=1))
model1i = nls(theoretic_ki~a*sqrt(t)+d,start=list(a=1, d=1))
model2i = nls(theoretic_ki~a*((t)^b)+d, start = list(a=0.1, b=0.1, d=0.1), algorithm = "port")
model3i = nls(theoretic_ki~d+a*exp(c*t), start = list(a=1,c=0.001, d=0.1),
            control = nls.control(maxiter = 100, tol=0.4))
model4i = nls(theoretic_ki~a*(log(abs(t+d1)))+d2, start=list(a=0.1 ,d1=1, d2=1), control = nls.control(tol = 0.7))

wt_inter = list(model0, model1, model2, model3, model4)
w_inter = list(model0i, model1i, model2i, model3i, model4i)

# get the sum of residuals squared of the models and the aic value for a model
get_RSS = function(x) return(sum(theoretic_ki-predict(x))^2)
get_AIC = function(rss, par) return(t.max*log(2*pi) + t.max*log(rss/t.max) + t.max + 2*(par + 1))
get_params = function(x) return(as.vector(summary(x)$parameters[,1]))

# parameters of each model
mod_params = c(1,1,2,2,2)
# get parameter estimate from models..
pm = lapply(wt_inter, get_params)
# replace missing parameters with NA (some models have 1 params others 2)..
pm = lapply(pm, `length<-`, max(lengths(pm)))
# transform it in matrix..
pm = do.call(rbind, pm)
# residual sum of squares
rss = sapply(wt_inter, get_RSS)
rss_i = sapply(w_inter, get_RSS)
# Akaike information criterion
aic_vect = mapply(get_AIC,rss, mod_params)
aic_vect_i = mapply(get_AIC, rss_i, mod_params+1) #+1 because all the models have now the intercept


# model selection growth + preferential attachment 
ms_gp = data.frame("RSS no inter"= rss, "AIC no inter"=aic_vect, "RSS inter"=rss_i, "AIC inter"= aic_vect_i)
rownames(ms_gp) = c("linear regr", "plaw 0.5", "plaw", "expo", "log")
ms_gp[which.min(ms_gp$AIC.no.inter),]
