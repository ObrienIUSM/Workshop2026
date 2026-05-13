library(FastPG)

in_path  <- "C:/Users/ananamat/Desktop/Obrien webinar"
out_path <- "C:/Users/ananamat/Desktop/Obrien webinar"
data_scaled <- readRDS(file.path(in_path, "scaled_umap_clusters_30um.rds"))

# Subset marker columns 
markers <- c("MPO","CD8","CD3","CD45RO","CD20","CD4","Ki67","CD45","CD11c","GATA3","Vimentin","FOXP3","SLC12A3","OPN","PLSCR1","Nestin","VCAM1","KIM1","AQP1","CD90","pMLKL","CD123","CD56","WT_1","pcJUN","MCP1","Cystatin3","LC3","Fibronectin","PROM1","CD31","NaKATPase","IGFBP7","Ecadherin","Podocalyxin","CollagenIV","Bcatenin","LRP2","CD68","aSMA","AQP2","CD206","HLA.DR","UMOD","Cytokeratin8","Biopsy","Object","X","Y","cluster","umap1","umap2", "n_neighbors", "cnt_0","cnt_1","cnt_2","cnt_3","cnt_4","cnt_5","cnt_6","cnt_7","cnt_8","cnt_9","cnt_10","cnt_11","cnt_12","cnt_13","cnt_0_ave0","cnt_1_ave0","cnt_2_ave0","cnt_3_ave0","cnt_4_ave0","cnt_5_ave0","cnt_6_ave0","cnt_7_ave0","cnt_8_ave0","cnt_9_ave0","cnt_10_ave0","cnt_11_ave0","cnt_12_ave0","cnt_13_ave0")

scaled.subset <- data_scaled[, markers] 

use_cols <- c("cnt_0_ave0","cnt_1_ave0","cnt_2_ave0","cnt_3_ave0","cnt_4_ave0","cnt_5_ave0","cnt_6_ave0","cnt_7_ave0","cnt_8_ave0","cnt_9_ave0","cnt_10_ave0","cnt_11_ave0","cnt_12_ave0","cnt_13_ave0")

data.mat <- as.matrix(scaled.subset[, use_cols])

# Run FastPG clustering
system.time({clusters <- fastCluster(data.mat, k = 500, num_threads = 20)})
scaled.subset$NE_cluster <- as.factor(clusters$communities)

saveRDS(scaled.subset, file.path(out_path, "scaled_umap_clusters_30um_NEclusters.rds"))

# Cell composition
library(dplyr)
library(tidyr)

df <- scaled.subset

prop_long <- df %>%
  group_by(NE_cluster, cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(NE_cluster) %>%
  mutate(total = sum(n), pct = 100 * n / total) %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = unique(cluster[order(as.integer(sub(".*_", "", cluster)))]))) %>%
  arrange(NE_cluster, cluster)

prop_wide <- prop_long %>%
  select(NE_cluster, cluster, pct) %>%
  pivot_wider(names_from = cluster, values_from = pct, values_fill = 0)
write.csv(prop_wide, file = file.path(out_path, "NEcluster_composition.csv"), row.names = FALSE)
