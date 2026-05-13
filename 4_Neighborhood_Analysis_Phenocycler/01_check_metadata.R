# Check metadata
df <- readRDS("C:/Users/ananamat/Desktop/Obrien webinar/scaled_umap_clusters")

# Dimensions (rows, cols)
dim(df)
# Column names
colnames(df)
# Count each category
table(df$Biopsy)

