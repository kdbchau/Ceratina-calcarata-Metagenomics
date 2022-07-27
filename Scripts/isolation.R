library(raster)
library(gdistance)
library(gstudio)
library(tidyverse)
library(ggpubr)
library(geosphere)
library(ecodist)

raster_landuse <- raster("../output_compressed.tif") #Land Use Tiff file 
e <- extent( c(-80.05,-78.7,43.43,44.1)) #Crop to study area 
raster_landuse_c <- crop(raster_landuse, e)

plot(raster_landuse_c)

# Open a png file
png("raster_landuseGTA.png")
# 2. Create a plot
plot(raster_landuse_c, xlab="Longitude", ylab="Latitude")
# Close the png file
dev.off()

classes <- table(getValues(raster_landuse_c)) # Tabulation of class counts
colors <- raster_landuse_c@legend@colortable  # Class colors


png("raster_landuseGTA_legend.png")
# 2. Create a plot
plot(raster_landuse_c,col=colors[2:10], colNA=colors[1], xlab="Longitude", ylab="Latitude")
# Close the png file
dev.off()

### These values will change depending on hypothesis
### Resistanc values indicated in Table 1 
tmp <- raster_landuse_c
tmp[ raster_landuse_c == 1 ] <- 0.9 #Water
tmp[ raster_landuse_c == 2 ] <- 0.1 #Trees
tmp[ raster_landuse_c == 3 ] <- 0.1 #Grass
tmp[ raster_landuse_c == 4 ] <- 0.9 #Flooded vegetation
tmp[ raster_landuse_c == 5 ] <- 0.1 #Crops
tmp[ raster_landuse_c == 6 ] <- 0.1 #Scrub/shrub
tmp[ raster_landuse_c == 7 ] <- 0.9 #Built Area
tmp[ raster_landuse_c == 8 ] <- 0.9 #Bare ground

cat_landuse <- ratify(tmp, count=TRUE)
rat <- levels(cat_landuse)[[1]]
levels(cat_landuse) <- rat

png("raster_landuse_cost.png")
plot(cat_landuse, col = terrain.colors(255), xlab="Longitude", ylab="Latitude")
dev.off()

tr <- transition( 1/cat_landuse, transitionFunction = mean, directions = 4 )
tr <- geoCorrection( tr, type="c", multpl=FALSE, scl=FALSE)
points <- read_tsv("../samples_map_162_sorted_alphabetically.txt") %>%
  select(Sample, Lat, Long)
sp_points <- SpatialPoints(points[,3:2])

png("samples_on_landuse_terrain.png")
plot(cat_landuse, col = terrain.colors(255), xlab="Longitude", ylab="Latitude")
plot(sp_points, add = TRUE)
dev.off()

### Example shortest path between two random points 
png("ShortestPathExample.png")
path.1 <- shortestPath( tr, sp_points[112], sp_points[83], output="SpatialLines")
plot(tmp , xlab="Longitude", ylab="Latitude", col = terrain.colors(255))
lines( path.1, col="red")
points(points[c(112,83),3:2],pch=16, col="red")
dev.off()

### Calculating distance using least-cost path 

eDist <- costDistance(tr, sp_points)
eDist <- as.matrix(eDist)
rownames(eDist) <- colnames(eDist) <- points$Sample
eDist[1:10,1:10]

geoDist <- distm(as.matrix(points[3:2]), fun = distGeo)

df <- data.frame( Landuse_Distance = eDist[ lower.tri(eDist)],
                  Geo_Distance = geoDist[ lower.tri(geoDist)])

df <- df[ !is.infinite(df$Landuse_Distance),]

write.table(df, file='distances_hyp1.tsv', quote=FALSE, sep='\t', col.names = NA)


#load relatedness file (comparison to self removed, 13041 pairwise comparisons, double check that the order is the same as in geo distance)
rel <- read_tsv("../LCPath_results/posthwe.relatedness.noself")
hyp1_df <- cbind(landuse_hyp1_df, CS_hyp1, rel)
MRM(hyp1_df$RELATEDNESS_AJK ~ hyp1_df$Landuse_Distance_Hyp1 + hyp1_df$Geo_Distance, nperm = 1000)

ggplot(hyp1_df,aes(x=Landuse_Distance_Hyp1,y=RELATEDNESS_AJK)) + geom_point(alpha = 0.3) + 
  stat_smooth(method=lm, formula = y ~ x) + xlab("Resistance Distance (Land use)") + ylab("Yang's Relatedness") + theme_classic() + stat_cor(method = "pearson", label.x = 0.8, label.y = 0.09)

### Export raster files
writeRaster(cat_landuse, file = "cat_landuse_hyp1.asc")
sites <- rasterize(sp_points, cat_landuse)
writeRaster(sites,filename="sites_hyp1.asc")   

### The above script used for other hypothesis, except resistance values changed accordingly

### Environmental distance script 

points <- read_tsv("../samples_map_162_sorted_alphabetically.txt") %>% select(Sample, Lat, Long)

envdata <- read_tsv("../env_variables_162_samples_standardized_temp_prec.txt")

envdatavar <- envdata %>% select(Annual_Temperature_Standardized, Annual_Precipitation_Standardized)
env_dist_canb <- dist(envdatavar, method = "canberra")
env_dist_canb_matrix <- as.matrix(env_dist_canb)

geoDist <- distm(as.matrix(points[3:2]), fun = distGeo)


env_df <- data.frame( Environmental_Distance = env_dist_canb_matrix[lower.tri(env_dist_canb_matrix)],
                  Geo_Distance = geoDist[ lower.tri(geoDist)])



