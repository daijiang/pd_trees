# phylogenies ----
tree_list = list(purpose = multi2di(tree_marx),
                      apg = multi2di(tree_apg),
                      zanne = multi2di(tree_zanne),
                      otl = multi2di(tree_otl),
                      random = multi2di(tree_random))

# simulations ----
library(ape)
library(phylolm)
library(plyr)
library(dplyr)
library(mvMORPH)
library(tidyverse)

# trait simulation
nsim = 100
true_var = c(0.2, 0.75, 1.5)
true_alpha = c(0.005, 0.05, 0.5, 1)
true_model = c("BM", "OU")
n_row_bm = length(true_var) * nsim * length(tree_list) # var * nsim * ntree
n_row_ou = length(true_var) * nsim * length(tree_list) * length(true_alpha) # var * nsim * ntree * alpha
results_df = tibble(true_model = c(rep("BM", n_row_bm), rep("OU", n_row_ou)),
                        true_var = rep(rep(true_var, each = nsim * length(tree_list)),  1 + length(true_alpha)),
                        true_alpha = c(rep(0, n_row_bm), rep(true_alpha, each = n_row_ou/length(true_alpha))),
                        nitem = rep(rep(1:nsim, each = length(tree_list)), (n_row_bm + n_row_ou)/(nsim * length(tree_list))),
                        tree = rep(names(tree_list), (n_row_bm + n_row_ou)/length(tree_list)),
                        seed = rep(1:((n_row_bm + n_row_ou)/length(tree_list)), each = length(tree_list)))

phy_sig_1 = function(x, tree_list){
  set.seed(x$seed[1] + 1000)
  # simulate trait with tree_purpose
  if(x$true_model[1] == "BM"){
    trait_sim = rTraitCont(tree_list$purpose, model = "BM", sigma = sqrt(x$true_var[1]))
    trait_sim = mvMORPH::mvSIM(tree = tree_list$purpose, nsim = 1, model = "BM1", 
                               param = list(ntrait = 1, theta = 0, sigma = sqrt(x$true_var[1])))[,1]
    
  }
  
  if(x$true_model[1] == "OU"){
    trait_sim = rTraitCont(tree_list$purpose, model = "OU", sigma = sqrt(x$true_var[1]),
                           alpha = x$true_alpha[1], theta = 0)
  }
  
  plyr::adply(x, 1, function(xx){
    z0 <- try(phylolm(trait_sim ~ 1, phy = tree_list[[xx$tree]], model = "BM"))
    z1 <- try(phylolm(trait_sim ~ 1, phy = tree_list[[xx$tree]], model = "lambda", 
                      upper.bound = 1, lower.bound = 0))
    z2 <- try(phylolm(trait_sim ~ 1, phy = tree_list[[xx$tree]], model = "OUfixedRoot", 
                      upper.bound = 30, lower.bound = 0))
    data_frame(phySigLambda = if(inherits(z1, "try-error")) NA else z1$optpar, 
               phySigK = phytools::phylosig(tree_list[[xx$tree]], trait_sim, method = "K"),
               sigma2_est_bm = if(inherits(z0, "try-error")) NA else z0$sigma2, 
               sigma2_est_ou = if(inherits(z2, "try-error")) NA else z2$sigma2,
               aic_bm = if(inherits(z0, "try-error")) NA else z0$aic, 
               aic_ou = if(inherits(z2, "try-error")) NA else z2$aic, 
               alpha_est_ou = if(inherits(z2, "try-error")) NA else z2$optpar,
               loglik_em = if(inherits(z0, "try-error")) NA else z0$logLik, 
               loglik_ou = if(inherits(z2, "try-error")) NA else z2$logLik)
  })
}

results_sig = group_by(results_df, true_model, true_var, true_alpha, nitem) %>% 
  do(phy_sig_1(., tree_list = tree_list_marx))
