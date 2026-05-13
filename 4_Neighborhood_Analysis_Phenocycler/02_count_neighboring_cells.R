library(dbscan); library(data.table)

in_path  <- "C:/Users/ananamat/Desktop/Obrien webinar"
out_path <- "C:/Users/ananamat/Desktop/Obrien webinar"
data_scaled <- readRDS(file.path(in_path, "scaled_umap_clusters"))

# Function: neighborhood counts
compute_neighborhood_counts_by_sample <- function(data, radius, cluster_col = "cluster"
){
  dt <- as.data.table(copy(data))
  stopifnot(all(c("Object","Biopsy","X","Y", cluster_col) %in% colnames(dt)))
  dt[[cluster_col]] <- as.character(dt[[cluster_col]])
  clusters <- as.character(0:13)
  res_all <- list()
  for(s in unique(dt$Biopsy)){
    message("Processing Biopsy: ", s)
    sub <- dt[Biopsy == s]
    coords <- as.matrix(sub[, .(X, Y)])
    fr <- frNN(coords, eps = radius)
    res <- lapply(seq_len(nrow(sub)), function(i){
      idx <- fr$id[[i]]
      idx <- idx[idx != i]
      k <- length(idx)
      cnts <- if(k == 0){
        setNames(integer(length(clusters)), paste0("cnt_", clusters))
      } else {
        tab <- table(factor(sub[[cluster_col]][idx], levels = clusters))
        setNames(as.integer(tab), paste0("cnt_", clusters))}
      c(list(Object = sub$Object[i], Biopsy = s, X = sub$X[i], Y = sub$Y[i], n_neighbors = k),cnts)})
    res_all[[s]] <- rbindlist(lapply(res, as.data.table), fill = TRUE)}
  rbindlist(res_all, fill = TRUE)}

# Run
radius_pixel <- 58.8
neighbor_counts <- compute_neighborhood_counts_by_sample(data = data_scaled, radius = radius_pixel)

data_scaled_updated <- merge(data_scaled, neighbor_counts, by = c("Object", "Biopsy", "X", "Y"), all.x = TRUE, sort = FALSE)

# add mean-zero columns (cnt_L3_*_ave0) with tiny jitter
cnt_cols <- paste0("cnt_", 0:13)
nvec <- data_scaled_updated$n_neighbors
set.seed(1)
eps <- 1e-6
for (col in cnt_cols) {
  prop <- ifelse(nvec == 0, 0, data_scaled_updated[[col]] / nvec)
  prop <- prop + rnorm(length(prop), sd = eps)   
  data_scaled_updated[[paste0(col, "_ave0")]] <- prop - mean(prop)}

out_rds_file <- file.path(in_path, "scaled_umap_clusters_30um.rds")
fwrite(data_scaled_updated, file.path(out_path, "neighbor_counts.csv"))
saveRDS(data_scaled_updated, out_rds_file)


