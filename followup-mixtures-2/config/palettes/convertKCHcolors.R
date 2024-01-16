library(R.matlab)
library(tidyverse)
library(foreach)

# Read in Matlab object with colors.
palette <- readMat("config/palettes/KCHcolors-gg.mat")

# Iterate through each level and extract the colors to format as a dataframe.
hexpalette <- foreach(level=seq(2,6), .combine="rbind") %do% {
  # Extract colors and taxa from Matlab object.
  colors <- as.data.frame(palette[["colors.kat"]][[level]])
  colnames(colors) <- c("R","G","B")
  colors$taxa <- unlist(palette[["OTUs.kat"]][[level]])
  # Convert colors from RGB to hex.
  colors <- colors %>%
    mutate(hex=rgb(R,G,B), level=level)
  colors
}

# Extract a shortened version of taxa names.
# Note that these shortened taxa names may not be unique,
# i.e. there are taxa of unknown genera in most orders.
hexpalette <- hexpalette %>%
  mutate(taxashort=gsub(".*__","",taxa))

# Export R-based color palette.
write.table(hexpalette %>% 
              dplyr::select(level, taxa, taxashort, hex, R, G, B),
            "config/palettes/KCHcolors-gg.txt", row.names=FALSE, quote=TRUE)
