library(ashr)
library(mashr)

#Simulate some data
set.seed(1)
simdata = simple_sims(500,5,1)

#Construct mashr object
data = set_mash_data(simdata$Bhat, simdata$Shat)

#Set up canonical covariance matrices
U.c = cov_canonical(data)  
print(names(U.c))

#Fit the model
m.c = mash(data, U.c, optmethod = "mixEM")

#Get pairwaise sharing
print(get_pairwise_sharing(m.c, factor = 0.5))

#Get mixture proportions
print(get_estimated_pi(m.c))


### Data-driven analysis ####
#Run condition-by-condition
data = set_mash_data(simdata$Bhat, simdata$Shat)
m.1by1 = mash_1by1(data)
strong = get_significant_results(m.1by1,0.05)

U.pca = cov_pca(data,5,strong)
print(names(U.pca))

U.ed = cov_ed(data, U.pca, strong)
m.ed = mash(data, U.ed)
print(get_loglik(m.ed),digits = 10)

print(get_pairwise_sharing(m.ed, factor = 0.5))

mash_plot_meta(m.c,get_significant_results(m.c)[1])

