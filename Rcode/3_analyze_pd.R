index_list_alpha = list(#pd.root = "pd.root", 
  pd.uroot = "pd.uroot", pd.uroot.z = "pd.uroot.z",
  mpd = "mpd", mpd.z = "mpd.z", 
  mntd = "mntd", mntd.z = "mntd.z",
  psv = "psv", 
  # pse = "pse", 
  vpd = "vpd")

index_list_beta = list(unif = "unif",
                       mpd_beta = "mpd_beta",
                       mntd_beta  = "mntd_beta",
                       pcd_beta = "pcd_beta",
                       pcdc_beta = "pcdc_beta",
                       pcdp_beta = "pcdp_beta",
                       mpd_beta_z = "mpd_beta_z",
                       mntd_beta_z = "mntd_beta_z")
# dataset_name = "pine"
for(dataset_name in c("pine", "marx", "fl")){
  cat(dataset_name)
  tree_build = paste0("tree_", dataset_name)
  
  results_path = paste0("results/", dataset_name)
  if(!dir.exists(results_path)) dir.create(results_path)
  
  results = readRDS(paste0("data/output/", dataset_name, "/results_", dataset_name, ".rds"))
  
  # long format results ----
  res_list = map(results, ~map(., bind_rows, .id = "tree")) %>% transpose() %>% 
    map(bind_rows, .id = "dat")
  
  # alpha results
  # correlations
  # for each dataset, cal cor among 5 trees, for 5+ pd measures, take a while.
  if(!file.exists(paste0(results_path, "/results_alpha.rds"))){
    res_sum_alpha = res_list$out_alpha %>% 
      group_by(dat) %>%
      do(cor_table2(df = ., index_list = index_list_alpha, tree_seq = tree_build)) %>% 
      ungroup()
    unique(sapply(res_sum_alpha$cor_models, class))
    saveRDS(res_sum_alpha, file = paste0(results_path, "/results_alpha.rds"))
  } else {
    res_sum_alpha = readRDS(paste0(results_path, "/results_alpha.rds"))
  }
  
  if(!file.exists(paste0(results_path, "/cor_alpha.rds"))){
    cor_alpha1 = select(res_sum_alpha, dat, index, tree1, tree2, same_tree, ord, cor_models) %>% 
      filter(same_tree == FALSE) %>% 
      mutate(cor_coef = map(cor_models, broom::tidy)) %>% 
      select(-cor_models) %>% 
      unnest(cor_coef) %>% select(dat, index, tree1, tree2, ord, estimate, p.value) %>% 
      mutate(estimate = round(estimate, 6), p.value = round(p.value, 5)) %>% 
      group_by(index, tree1, tree2) %>% 
      summarise(ave_cor = mean(estimate, na.rm = T),
                ave_p = mean(p.value, na.rm = T),
                sig_p_prop = mean(p.value < 0.05, na.rm = T), # proportion of significant correlations
                median_cor = median(estimate, na.rm = T),
                median_p = median(p.value, na.rm = T)) %>% ungroup()
    cor_alpha = filter(res_sum_alpha, dat == "dat_1") %>%
      select(index, tree1, tree2, ord, same_tree) %>% 
      left_join(cor_alpha1, by = c("index", "tree1", "tree2")) %>% 
      mutate(ave_cor = ifelse(is.na(ave_cor), 1, ave_cor),
             median_cor = ifelse(is.na(median_cor), 1, median_cor)) %>%
      arrange(index, ord) %>% 
      mutate(tree1 = factor(tree1, levels = unique(tree2)), 
             tree2 = factor(tree2, levels = unique(tree2)),
             ave_cor = round(ave_cor, 2),
             dataset = dataset_name) 
    saveRDS(cor_alpha, file = paste0(results_path, "/cor_alpha.rds"))
  } else {
    cor_alpha = readRDS(paste0(results_path, "/cor_alpha.rds"))
  }
  
  # beta diversity ----
  res_list$out_beta
  # correlations
  if(!file.exists(paste0(results_path, "/results_beta.rds"))){
    res_sum_beta = mutate(filter(res_list$out_beta, site1 != "multi_sites"), 
                          rowid = paste0(site1, site2)) %>% 
      group_by(dat) %>% 
      do(cor_table2(df = ., index_list = index_list_beta, tree_seq = tree_build)) %>% 
      ungroup()
    saveRDS(res_sum_beta, file = paste0(results_path, "/results_beta.rds"))
  } else {
    res_sum_beta = readRDS(paste0(results_path, "/results_beta.rds"))
  }
  
  if(!file.exists(paste0(results_path, "/cor_beta.rds"))){
    cor_beta1 = select(res_sum_beta, dat, index, tree1, tree2, same_tree, ord, cor_models) %>% 
      filter(same_tree == FALSE) %>% 
      mutate(cor_coef = map(cor_models, broom::tidy)) %>% 
      select(-cor_models) %>% 
      unnest(cor_coef) %>% select(dat, index, tree1, tree2, ord, estimate, p.value) %>% 
      mutate(estimate = round(estimate, 6), p.value = round(p.value, 5)) %>% 
      group_by(index, tree1, tree2) %>% 
      summarise(ave_cor = mean(estimate, na.rm = T),
                ave_p = mean(p.value, na.rm = T),
                sig_p_prop = mean(p.value < 0.05, na.rm = T), # proportion of significant correlations
                median_cor = median(estimate, na.rm = T),
                median_p = median(p.value, na.rm = T)) %>% ungroup()
    cor_beta = filter(res_sum_beta, dat == "dat_1") %>%
      select(index, tree1, tree2, ord, same_tree) %>%
      left_join(cor_beta1, by = c("index", "tree1", "tree2")) %>%
      mutate(ave_cor = ifelse(is.na(ave_cor), 1, ave_cor),
             median_cor = ifelse(is.na(median_cor), 1, median_cor)) %>%
      arrange(index, ord) %>%
      mutate(tree1 = factor(tree1, levels = unique(tree2)),
             tree2 = factor(tree2, levels = unique(tree2)),
             ave_cor = round(ave_cor, 2),
             dataset = dataset_name)
    saveRDS(cor_beta, file = paste0(results_path, "/cor_beta.rds"))
  } else {
    cor_beta = readRDS(paste0(results_path, "/cor_beta.rds"))
  }
  
  # lmm alpha -------
  res_list$out_alpha 
  # use the first 40 datasets
  if(!file.exists(paste0(results_path, "/lmm_alpha_models.rds"))){
    res_alpha_lmm = lmm_table2(res_list$out_alpha,
                               index_list = index_list_alpha,
                               tree_seq = tree_build, REML = F)
    saveRDS(res_alpha_lmm, file = paste0(results_path, "/lmm_alpha_models.rds"))
  } else {
    res_alpha_lmm = readRDS(paste0(results_path, "/lmm_alpha_models.rds"))
  }

  if(!file.exists(paste0(results_path, "/lmm_alpha.rds"))){
    res_alpha_lmm2 = select(res_alpha_lmm, -subdat) %>% 
      mutate(lmm_coef = map(lmm_models, get_lmm_coef)) %>% 
      select(-lmm_models) %>% 
      unnest(lmm_coef)
    
    # lmm tables
    res_alpha_lmm3 = select(res_alpha_lmm2, index, tree1, tree2, # intercept, 
                            ends_with("sd"), starts_with("R2"), starts_with("slope")) 
    # res_alpha_lmm3[, 4:12] = round(res_alpha_lmm3[, 4:12], 4)
    res_alpha_lmm3$dataset = dataset_name
    saveRDS(res_alpha_lmm3, file = paste0(results_path, "/lmm_alpha.rds"))
  }
  
  # lmm beta ------
  res_list$out_beta 
  # too many data point
  set.seed(123456)
  # if(!file.exists(paste0(results_path, "/lmm_beta_models.rds"))){
    system.time(res_beta_lmm <- lmm_table2(filter(res_list$out_beta, 
                                      # dat %in% paste0("dat_", sample(1000, 100)),
                                      site1 != "multi_sites"
                                      ), 
                               index_list = index_list_beta[1:4],
                               beta = TRUE, REML = F, tree_seq = tree_build))
    # saveRDS(select(res_beta_lmm, -subdat), file = paste0(results_path, "/lmm_beta_models.rds"))
  # } else {
  #   res_beta_lmm = readRDS(paste0(results_path, "/lmm_beta_models.rds"))
  # }

  if(!file.exists(paste0(results_path, "/lmm_beta.rds"))){
    res_beta_lmm2 = filter(res_beta_lmm, index != "pcdc_beta") %>% 
      mutate(lmm_coef = map(lmm_models, get_lmm_coef_beta)) %>% 
      unnest(lmm_coef)
    # res_beta_lmm2$lmm_models[[3]] %>% rr2::R2.lr()
    
    # lmm tables
    res_beta_lmm3 = select(res_beta_lmm2, index, tree1, tree2, #intercept, 
                           starts_with("slope"), starts_with("R2")) 
    # res_beta_lmm3[, 4:8] = round(res_beta_lmm3[, 4:8], 4)
    res_beta_lmm3$dataset = dataset_name
    as.data.frame(res_beta_lmm3)
    
    saveRDS(res_beta_lmm3, file = paste0(results_path, "/lmm_beta.rds"))
  }
  
}
