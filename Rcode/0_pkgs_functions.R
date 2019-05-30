library(ape)
library(picante)
library(brranching)
library(phylocomr)
if(!require(dli55)) devtools::install_github("daijiang/dli55")
# library(dli55)
library(parallel)
library(tidyverse)

# simulate community data
# I fixed species richness to the same for all sites, determined by 
# nspp = Ntip(tree) * sp_prop
sim_comm = function(nsite = 50, tree, lambda = 10, sp_prop = 0.1, sp_rich = NULL){
  nspp = Ntip(tree)
  sp_rich = ifelse(is.null(sp_rich), floor(nspp * sp_prop), sp_rich)
  comm_sim = matrix(0, nrow = nsite, ncol = nspp)
  for(i in 1:nsite){
    obs_idx = sample(x = nspp, size = sp_rich) # fix sp rich for all sites
    # print(obs_idx)
    obs_abund = rpois(n = length(obs_idx), lambda = lambda) # abundance data
    if(any(obs_abund == 0)) obs_abund[obs_abund == 0] = 1
    # print(obs_abund)
    comm_sim[i, obs_idx] = obs_abund
  }
  comm_sim = as.data.frame(comm_sim)
  names(comm_sim) = tolower(tree$tip.label)
  row.names(comm_sim) = paste("s", 1:nsite, sep = "_")
  dli55::rm_sp_noobs(comm_sim)
}

# to look at correlations among multiple phylo diversity measurements
cor_table = function(df, ind = "pd", tree_seq = "tree_pine"){
  # print(df$dat[1])
  indc = quo(ind)
  dat = select(df, tree, rowid, pds = !!indc) %>% spread("tree", "pds") %>% select(-rowid)
  x = names(dat)
  x = c(x[-which(x == "tree_random")], x[which(x == "tree_random")]) # put tree_random as the last one
  x = c(x[which(x == tree_seq)], x[-which(x == tree_seq)]) # put tree_seq as the first one
  
  xm = matrix(data = 0, nrow = length(x), ncol = length(x))
  dimnames(xm) = list(x, x)
  xm[lower.tri(xm)] = NA
  all_comb = reshape2::melt(xm, na.rm = TRUE) %>% select(-value) %>% 
    rename(tree1 = Var1, tree2 = Var2) %>% mutate(tree1 = as.character(tree1),
                                                  tree2 = as.character(tree2),
                                                  same_tree = tree1 == tree2,
                                                  ord = seq_along(tree1)) %>% 
    as_tibble()
  
  out1 = mutate(filter(all_comb, same_tree == FALSE), 
                subdat = map2(.x = tree1, .y = tree2, .f = ~ select(dat, .x, .y))) %>% 
    mutate(cor_models = map(subdat, ~try(cor.test(unlist(.[,1]), unlist(.[,2]), method = "spearman"))))
  # use LMM should be better than individual lms
           # lm_models = map(subdat, ~lm(as.formula(paste(names(.), collapse = " ~ ")), data = .)), # force the line go through the origin can be misleading if x and y are very different
           # lm_coefs = map(lm_models, broom::tidy),
           # r_sqr = map(lm_models, ~round(summary(.)$r.squared, 5))) %>% unnest(r_sqr)
  # diag
  out2 = filter(all_comb, same_tree == TRUE)
  bind_rows(out1, out2)
}

cor_table2 = function(df, index_list = list(pd = "pd", pd.prop = "pd.prop", mpd = "mpd",
                                            vpd = "vpd", mntd = "mntd", psv = "psv",
                                            pse = "pse"), ...){
  map(index_list, cor_table, df = df, ...) %>% 
    bind_rows(.id = "index") 
}

lmm_table = function(df, ind = "pd", beta = FALSE, tree_seq = "tree_pine", ...){
  indc = quo(ind)
  if(beta){
    dat = select(df, dat, tree, site1, site2, pds = !!indc) %>% spread("tree", "pds")
  } else {
    dat = select(df, dat, tree, rowid, pds = !!indc) %>% spread("tree", "pds") %>% select(-rowid)
  }
  
  all_comb = data_frame(tree1 = rep(tree_seq, 4), 
                        tree2 = c("tree_random", "tree_apg", "tree_zanne", "tree_otl"))
  
  if(beta){
    lmm = function(df, standardise = TRUE, origin = TRUE){
      if(standardise) df[,1:2] = scale(df[, 1:2])
      if(origin){
        xlm = lme4::lmer(as.formula(paste0(paste(names(df)[1:2], collapse = " ~ "), 
                                           " + (1|dat) - 1 + (1|dat:site1) + (1|dat:site2) + (0 + ", names(df)[2], "|dat)")), 
                         data = df)
      } else {
       xlm = lme4::lmer(as.formula(paste0(paste(names(df)[1:2], collapse = " ~ "), 
                                     " + (1|dat) + (1|dat:site1) + (1|dat:site2) + (0 + ", names(df)[2], "|dat)")), 
                   data = df)
      }
     xlm
    }
    out1 = mutate(all_comb, 
                  subdat = map2(.x = tree1, .y = tree2, .f = ~ select(dat, .x, .y, dat, site1, site2))) %>% 
      mutate(lmm_models = map(subdat, lmm, standardise = TRUE))
  } else {
    lmm = function(df, standardise = TRUE, origin = TRUE){
      if(standardise) df[,1:2] = scale(df[, 1:2])
      if(origin){
       xlm =  lme4::lmer(as.formula(paste0(paste(names(df)[1:2], collapse = " ~ "), 
                                     " + (1|dat) - 1 + (0 + ", names(df)[2], "|dat)")), 
                   data = df) # go through origin
      } else {
       xlm = lme4::lmer(as.formula(paste0(paste(names(df)[1:2], collapse = " ~ "), 
                                     " + (1|dat) + (0 + ", names(df)[2], "|dat)")), 
                   data = df)
      }
      xlm
    }
    out1 = mutate(all_comb, 
                  subdat = map2(.x = tree1, .y = tree2, .f = ~ select(dat, .x, .y, dat))) %>% 
      mutate(lmm_models = map(subdat, lmm, standardise = TRUE))
  }

  out1
}

lmm_table2 = function(df, index_list = list(pd = "pd", pd.prop = "pd.prop", mpd = "mpd",
                                            vpd = "vpd", mntd = "mntd", psv = "psv",
                                            pse = "pse"), 
                      beta = FALSE, ...){
  map(index_list, lmm_table, df = df, beta = beta, ...) %>% 
    bind_rows(.id = "index") 
}

get_lmm_coef = function(lmm, origin = TRUE){ # `origin` to be the same with lmm_table()!!
  df = broom::tidy(lmm)
  if(origin){
    df1 = select(df[1,], -group, -term)
    names(df1) = c("slope", "slope.se", "slope.t")
    df1$slope.re.sd = unlist(df[3, "estimate"])
    df1$resid.sd = unlist(df[4, "estimate"])
    
    ci = lme4::confint.merMod(lmm, parm = "beta_", method = "Wald")
    df1$slope_lower = ci[1]
    df1$slope_upper = ci[2]
    # df2 = rr2::R2(lmm)
    # df3 = as.data.frame(t(df2[, 2]))
    # names(df3) = df2$R2s
    # out = bind_cols(df1, df3)
    out = df1
  } else {
    df$term[1] = "Intercept"
    df1 = bind_cols(
      df[1:2, 1:2] %>% spread("term", "estimate") %>% setNames(c("intercept", "slope")),
      df[1:2, c(1, 3)] %>% spread("term", "std.error") %>% setNames(c("intercept.se", "slope.se")),
      df[1:2, c(1, 4)] %>% spread("term", "statistic") %>% setNames(c("intercept.t", "slope.t")),
      df[-c(1,2), 2, drop = F] %>% t() %>% as.data.frame() %>% 
        setNames(c("intercept.re.sd", "slope.re.sd", "resid.sd")))
    ci = lme4::confint.merMod(lmm, parm = "beta_", method = "Wald")[2,]
    df1$slope_lower = ci[1]
    df1$slope_upper = ci[2]
    df2 = rr2::R2(lmm)
    df3 = as.data.frame(t(df2[, 2]))
    names(df3) = df2$R2s
    out = bind_cols(df1, df3)
  }
  out
}

get_lmm_coef_beta = function(lmm, origin = TRUE){ # to be the same with lmm_table()!!
  df = broom::tidy(lmm)
  if(origin){
    df1 = select(df[1,], -group, -term)
    names(df1) = c("slope", "slope.se", "slope.t")
    df1$slope.re.sd = filter(df, group == "dat")$estimate
    df1$resid.sd = filter(df, group == "Residual")$estimate
    ci = lme4::confint.merMod(lmm, parm = "beta_", method = "Wald")
    df1$slope_lower = ci[1]
    df1$slope_upper = ci[2]
    # df2 = rr2::R2(lmm)
    # df3 = as.data.frame(t(df2[, 2]))
    # names(df3) = df2$R2s
    # out = bind_cols(df1, df3)
    out = df1
  } else {
    df$term[1] = "Intercept"
    df1 = bind_cols(
      df[1:2, 1:2] %>% spread("term", "estimate") %>% setNames(c("intercept", "slope")),
      df[1:2, c(1, 3)] %>% spread("term", "std.error") %>% setNames(c("intercept.se", "slope.se")),
      df[1:2, c(1, 4)] %>% spread("term", "statistic") %>% setNames(c("intercept.t", "slope.t")),
      df[-c(1,2), 2, drop = F] %>% t() %>% as.data.frame() %>% 
        setNames(c("icpt.site2.sd", "icpt.site1.sd", "slope.re.sd", "icpt.dat.sd", "resid.sd")))
    ci = lme4::confint.merMod(lmm, parm = "beta_", method = "Wald")[2,]
    df1$slope_lower = ci[1]
    df1$slope_upper = ci[2]
    df2 = rr2::R2(lmm, ce = F)
    df3 = as.data.frame(t(df2[, 2]))
    names(df3) = df2$R2s
    out = bind_cols(df1, df3)
  }
  out
}

# null.model = null.model, n.item = 0, abund.weight = TRUE, verbose = verbose
analyze_one = function(df, tree_list = c(tree_marx = tree_marx, 
                                         tree_zanne = tree_zanne,
                                         tree_apg = tree_apg,
                                         tree_otl = tree_otl,
                                         tree_random = tree_random),
                       alpha = TRUE, beta = TRUE, verbose = FALSE, ...){
  df_long = tibble::rownames_to_column(df, "site") %>% gather("sp", "freq", -site) %>% 
    filter(freq > 0) %>% arrange(site, sp) %>% select(site, freq, sp)
  
  # calculate PDs for each tree
  if(alpha){
    out_alpha = tree_list %>% 
      purrr::map(get_pd_alpha, samp_wide = df, samp_long = df_long,
          verbose = verbose, vpd = TRUE, ...) %>% 
      setNames(names(tree_list))
  }
  
  if(beta){
    out_beta = tree_list %>% 
      purrr::map(get_pd_beta, samp_wide = df, samp_long = df_long, 
          verbose = verbose, ...) %>% 
      setNames(names(tree_list))
  }
  if(alpha & beta){
    return(list(out_alpha = out_alpha, out_beta = out_beta))
  }
  
  if(alpha) return(list(out_alpha = out_alpha))
  if(beta) return(list(out_beta = out_beta))
}

analyze_interval = function(data = dat, begin = 1, n_to_analyze = 5, ncores = 5, 
                            treelist = c(tree_marx = tree_marx, 
                                         tree_zanne = tree_zanne,
                                         tree_apg = tree_apg,
                                         tree_otl = tree_otl,
                                         tree_random = tree_random),
                            ...){
  ends = begin + n_to_analyze - 1
  dat_test = data[begin:ends]
  mclapply(dat_test, function(x, ...){
    try(analyze_one(x, tree_list = treelist, alpha = TRUE, beta = TRUE, verbose = FALSE, ...))
  }, mc.cores = ncores)
}

#' calculate MPD and VPD (mean and variance of pairwise distance)
#' 
mvpd <- function(samp, dis, abundance.weighted=FALSE){
  N <- dim(samp)[1]
  mpd = numeric(N); vpd = numeric(N)
  for (i in 1:N) {
    # cat("row ", i)
    sppInSample <- names(samp[i, samp[i, ] > 0])
    if (length(sppInSample) > 1) {
      sample.dis <- dis[sppInSample, sppInSample]
      if (abundance.weighted) {
        sample.weights <- t(as.matrix(samp[i,sppInSample,drop=FALSE])) %*% as.matrix(samp[i,sppInSample,drop=FALSE])
        wm = weighted.mean(sample.dis,sample.weights)
        mpd[i] <- wm
        vpd[i] = sum(sample.weights[lower.tri(sample.weights)] * (sample.dis[lower.tri(sample.dis)] - wm)^2)
      }
      else {
        mpd[i] <- mean(sample.dis[lower.tri(sample.dis)])
        vpd[i] <- var(sample.dis[lower.tri(sample.dis)])
      }
    }
    else{
      mpd[i] = NA
      vpd[i] = NA
    }
  }
  data.frame(site = row.names(samp), mpd = mpd, vpd = vpd, stringsAsFactors = F)
}


