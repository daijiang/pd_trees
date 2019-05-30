source("Rcode/0_pkgs_functions.R")

# phylogeny ----
download.file(url = "http://datadryad.org/bitstream/handle/10255/dryad.151979/phy.Alpes.taxized.taxref.tre?sequence=1", destfile = "data/Marx.tre")
tree_marx = ape::read.tree("data/Marx.tre") %>% collapse.singles()
tree_marx$tip.label = tolower(tree_marx$tip.label)
tree_marx$node.label = paste0("node", tree_marx$node.label) # all numbers! cause problem with phylocom
# phylocomr:::write_tree_(tree_marx, digits = 6) %>%
#   cat(file = "data/tree_marx.tre", sep = "\n")
sp = tree_marx$tip.label

sp_1 = phylomatic_names(taxa = sp, db = "apg")
sp_2 = grep(pattern = "^NA", x = sp_1, value = T, ignore.case = T)
sp_3 = gsub(pattern = "NA/.+/(.*)$", replacement = "\\1", sp_2) %>%
  phylomatic_names()
sp_1[grep(pattern = "^NA", x = sp_1, value = F, ignore.case = T)] = sp_3
any(grep(pattern = "^NA", x = sp_1, value = T, ignore.case = T))

# zanne 2014 -----
# get phylogeny from Phylomatic
tree_zanne = phylomatic(taxa = sp_1, taxnames = F, get = "POST", storedtree = "zanne2014")
tree_zanne$tip.label = tolower(tree_zanne$tip.label)

# apg III ----
tree_apg = phylomatic(taxa = sp_1, taxnames = F, get = "POST", storedtree = "R20120829")
class(tree_apg) = "phylo" # remove phylomatic as a class
write.tree(tree_apg, file = "data/tree_apg.tre")
# no branch length
tree_apg_bladj = ph_bladj(ages = "data/ages", phylo = "data/tree_apg.tre")
tree_apg_bladj = read.tree(text = tree_apg_bladj)
tree_apg_bladj$tip.label = tolower(tree_apg_bladj$tip.label)
tree_apg = tree_apg_bladj

# random phylogeny ----
max(phytools::nodeHeights(tree_marx))
max(phytools::nodeHeights(tree_zanne))
max(phytools::nodeHeights(tree_apg_bladj))
tree_random = ape::rcoal(n = Ntip(tree_marx), tip.label = tree_zanne$tip.label) %>%
  geiger::rescale("depth", max(phytools::nodeHeights(tree_marx))) # coalescent tree

# open tree of life ----
library(rotl)
resolved_sp = tnrs_match_names_2(sp, context_name = "Land plants")
resolved_sp2 = bind_rows(resolved_sp)
resolved_sp3 = setdiff(tolower(sp), resolved_sp2$search_string) %>%
  gsub(pattern = "_subsp.*$", replacement = "", .) %>%
  gsub(pattern = "_var.*$", replacement = "", .) %>%
  gsub(pattern = "ervilia_sylvatica", replacement = "vicia_sylvatica", .) %>%
  tnrs_match_names(context_name = "Land plants") %>% unique()
dim(resolved_sp2)
sum(is.na(resolved_sp2$unique_name))
resolved_sp4 = filter(resolved_sp2, !is.na(unique_name)) %>%
  bind_rows(resolved_sp3)
sum(is.na(resolved_sp4$ott_id))
resolved_sp4$search_string[duplicated(resolved_sp4$search_string)]
filter(resolved_sp4, search_string %in% resolved_sp4$search_string[duplicated(resolved_sp4$search_string)])
resolved_sp4 = unique(resolved_sp4)
taxon_map = structure(resolved_sp4$search_string, names = resolved_sp4$unique_name)
taxon_map[taxon_map == "vicia_sylvatica"] = "ervilia_sylvatica" # put it back

tree_otl = try(tol_induced_subtree(ott_ids = resolved_sp4$ott_id, label_format = "name"))
filter(resolved_sp4, ott_id == 1066887)
# the following ott id were not in the tree synthesis database
no_id = c(643726, 1066887, 762660, 103745, 464805, 160980, 798537, 1025132, 4730993, 8597)
tree_otl = tol_induced_subtree(ott_ids = filter(resolved_sp4, !ott_id %in% no_id)$ott_id, label_format = "name")
tree_otl$tip.label[tree_otl$tip.label == "Vicia_sylvatica"] = "Ervilia_sylvatica"

# get node ages from timetree.org
# timetree.org does not allow programmingly explore their website
itree_otl$info %>% filter(!is.na(median)) %>%
  select(node, estimated) %>%
  write.table(file = "data/ages", quote = F, sep = "\t", row.names = F, col.names = F)

# call phylomcom bladj
system("cd data && ~/Documents/phylocom-4.2/src/phylocom bladj -f tree_otl_nobrlen.tre > tree_otl_brlen.tre")
tree_otl = read.tree("data/tree_otl_brlen.tre")
tree_otl$tip.label = tolower(tree_otl$tip.label)

# prune all phylogenies to have the same set of species
length(intersect(tree_marx$tip.label, tree_otl$tip.label))
tree_marx = drop.tip(tree_marx, tree_marx$tip.label[!tree_marx$tip.label %in% tree_otl$tip.label])
tree_zanne = drop.tip(tree_zanne, tree_zanne$tip.label[!tree_zanne$tip.label %in% tree_otl$tip.label])
tree_apg = drop.tip(tree_apg, tree_apg$tip.label[!tree_apg$tip.label %in% tree_otl$tip.label])
tree_random = drop.tip(tree_random, tree_random$tip.label[!tree_random$tip.label %in% tree_otl$tip.label])
tree_otl = drop.tip(tree_otl, tree_otl$tip.label[!tree_otl$tip.label %in% tree_zanne$tip.label])

## simulate marx datasets ----
n_sim = 1000
dat_marx = map(1:n_sim, function(x){
  set.seed(x)
  sim_comm(nsite = 30, tree = tree_marx, lambda = 6, sp_rich = 50)
})
names(dat_marx) = paste0("dat_", 1:n_sim)

for(i in seq(1, n_sim, 5)){
  cat("analyze ", i, " to ", i+4, "\n")
  testr = analyze_interval(dat_marx, begin = i, n_to_analyze = 5, ncores = 5,
                           treelist = c(tree_marx = tree_marx, 
                                        tree_zanne = tree_zanne,
                                        tree_apg = tree_apg,
                                        tree_otl = tree_otl,
                                        tree_random = tree_random),
                           abund.weight = FALSE, null.model.phylomeasures = TRUE,
                           n.item = 999, null.model.phylocom = TRUE, verbose = F)
  results_marx[i:(i+4)] = testr
  saveRDS(results_marx, file = "data/output/marx/results_marx.rds")
  cat("analyzed ", i, " to ", i+4, "at ", toString(Sys.time()), "\n", file = "marx.txt", append = T)
}

# saveRDS(results_marx, file = "data/results_marx.rds")
# results_marx = readRDS("data/results_marx.rds")
# map_chr(results_marx, class)

# # some simulations did not work, re-run it
# results_marx = readRDS("data/output/marx/results_marx.rds")
# which(map_int(results_marx, length) != 2)
# results_marx$dat_196$out_alpha$tree_marx
# results_marx$dat_156$out_alpha$tree_marx
# x$out_alpha$tree_marx
# 
# w_null = res_list$out_alpha %>% filter(is.na(sr)) %>% pull(dat) %>% unique() 
# 
# saveRDS(results_marx, file = "data/output/marx/results_marx.rds")
# 
# for(i in seq(156, 195, 5)){
#   cat("analyze ", i, " to ", i+4, "\n" )
#   testr = analyze_interval(dat_marx, begin = i, n_to_analyze = 5, ncores = 5,
#                            treelist = c(tree_marx = tree_marx, 
#                                         tree_zanne = tree_zanne,
#                                         tree_apg = tree_apg,
#                                         tree_otl = tree_otl,
#                                         tree_random = tree_random),
#                            abund.weight = FALSE, null.model.phylomeasures = TRUE,
#                            n.item = 999, null.model.phylocom = TRUE, verbose = F)
#   results_marx[i:(i+4)] = testr
#   cat("analyzed ", i, " to ", i+4, "at ", toString(Sys.time()), "\n", file = "marx.txt", append = T)
# }
# 
# saveRDS(results_marx, file = "data/output/marx/results_marx.rds")
