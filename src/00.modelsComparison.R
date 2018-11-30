rm(list = ls())
# Load and install necessary packages
requiredPackages <- c("igraph", "ggplot2", "Matrix", "ggthemes", "gridExtra",
                      "minpack.lm")

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
# if(TESTING){
#    theoretic_ki = c()
#    exact_scaled_therotic_ki = c()
#    # approximated vertex degree over time (slide 4)
#    approx_scaled_therotic_ki = rep(m.0*sqrt(t.max), t.max+1)
#    m.0 = 2
#    for (t in seq(t.max+1)) {
#       # vertex degree over time (slide 4)
#       theoretic_ki = append(theoretic_ki, (m.0*((t.max/t)**0.5)))
#       # scaled vertex degree over time (slide 4)
#       exact_scaled_therotic_ki = append(exact_scaled_therotic_ki, 
#                                         sqrt(t)*theoretic_ki[t])
#       
#    }
# }


test = as.data.frame(cbind(table.BA[,-1],
                           table.BA.no.growth[,-1],
                           table.BA.Rand[,-1]))

# model selection ---------------------------------------------------------

# sequence of values
t = seq(t.max+1)
# get the sum of residuals squared of the models,the aic value and the parmeters found

# parameters of each model
mod_params = c(1, 1, 2, 2, 2, 2, 2, 3, 3, 3)

prova = function(i){
   
# models without intercept
model0 = nlsLM(i ~ a*t, start=list(a=1))
model1 = nlsLM(i ~ a*sqrt(t),start=list(a=1))
model2 = nlsLM(i ~ a*(t)^b, start = list(a=0.1, b=1))
model3 = nlsLM(i ~ a*exp(c*t), start = list(a=1,c=0.0001))
model4 = nlsLM(i ~ a*(log(abs(t+d1))), start=list(a=0.1 ,d1=1))

# models with intercept
model0i = nlsLM(i ~ a*t+d,start=list(a=1, d=1))
model1i = nlsLM(i ~ a*sqrt(t)+d,start=list(a=1, d=1))
model2i = nlsLM(i ~ a*((t)^b)+d, start = list(a=0.1, b=1, d=1), control=nls.lm.control(maxiter = 150))
model3i = nlsLM(i ~ d+a*exp(c*t), start = list(a=1,c=0.01, d=0.1))
model4i = nlsLM(i ~ a*(log(abs(t+d1)))+d2, start=list(a=10 ,d1=1000, d2=-100))

# the models computer are now inserted in a list to compact the code
models_list = list(model0, model1, model2, model3, model4,
                   model0i, model1i, model2i, model3i, model4i)

get_RSS = function(x) return(sum(i-predict(x))^2) #check
get_AIC = function(rss, par) return(t.max*log(2*pi) + t.max*log(rss/t.max) + t.max + 2*(par + 1))
get_params = function(x) return(as.vector(summary(x)$coefficients[,1]))

# get parameter estimate from models..
pm = lapply(models_list, get_params)
# replace missing parameters with NA (some models have 1 params others 2, others 3)..
pm = lapply(pm, `length<-`, max(lengths(pm)))
# transform it into matrix to append in the dataframe
pm = do.call(rbind, pm)

# residual sum of squares
rss = sapply(models_list, get_RSS)

# Akaike information criterion
aic_vect = mapply(get_AIC,rss, mod_params)

model_name = c("linear regr", "plaw 0.5", "plaw", "expo", "log",
                    "linear regr+i", "plaw 0.5+i", "plaw+i", "expo+i", "log+i")
# model selection growth + preferential attachment 
ms_gp = data.frame("Model" = model_name,
                   "RSS"= rss,      # residual square of sum no intercept
                   "AIC"=aic_vect,  # aic no intercept
                   "Param1"=pm[,1],          # param 1
                   "Param2"=pm[,2],          # param 2
                   "Param3"=pm[,3])          # param 
# add model names to df

# show the one that fit better
better_fit = ms_gp[which.min(ms_gp$AIC),]

}

df = do.call(rbind, apply(test,2, prova))

a = as.numeric(better_fit["Param1"])
b = as.numeric(better_fit["Param2"])

plot(table.BA$ks.2,main = "Fitting best model found", ylim = c(1,200), xlab = "t", ylab = "k(t)")
curve(0.0006748097*x+6.963243, col="red", add = TRUE, lwd = 2)
curve(2*sqrt(x), col="green", add = TRUE, lwd=2)
legend("bottomright", legend = c("Empirical k_i", "Best fitting"),
       lty = 1, col = c("black", "red"))




