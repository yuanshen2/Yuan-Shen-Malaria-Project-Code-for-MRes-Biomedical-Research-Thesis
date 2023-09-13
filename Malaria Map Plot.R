library(sf)
#Plot a polygon
xcoords<- c(14,14,15,15,15)
ycoords <- c(4,5,5,4,4)
poly1 <- sp::Polygon(cbind(xcoords,ycoords))
firstPoly <- sp::Polygons(list(poly1), ID = "A")
firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
CRS_new<-CRS(SRS_string = "EPSG:4326")
proj4string(firstSpatialPoly) <- CRS_new
plot(firstSpatialPoly)

library(sf)

# Read the Africa shapefile
africa1 <- st_read(dsn = "C:/Users/Administrator/Documents/Africa Map Visualization/Africa.shp")

# Read the Africa+ shapefile
africa2 <- st_read(dsn = "C:/Users/Administrator/Documents/Africa Map Visualization/Africa+.shp")

# Identify columns unique to each dataset
cols_africa1 <- setdiff(names(africa1), names(africa2))
cols_africa2 <- setdiff(names(africa2), names(africa1))

# Add columns from africa1 to africa2 with NA values
for (col in cols_africa1) {
  africa2[[col]] <- NA
}

# Add columns from africa2 to africa1 with NA values
for (col in cols_africa2) {
  africa1[[col]] <- NA
}

# Ensure the column order is the same before using rbind
africa2 <- africa2[, names(africa1)]

# Use rbind to combine the datasets
africa_shp <- rbind(africa1, africa2)


# Plot DRC
library(sf)
DRC <- africa_shp[which(africa_shp$country_na == 'Democratic Republic of the Congo'),]
plot(st_geometry(DRC), col="red")
africa_shp$country_na

# Install and load required packages
install.packages(c("INLA", "cli", "rlang", "PrevMap", "ggplot2", "scales", "raster"))
library(INLA)
library(cli)
library(rlang)
library(PrevMap)
library(ggplot2)
library(scales)
library(raster)

# Load the .RData file
Malariadata <- inputs$data_all_wa_ea_drc
print(Malariadata)
Malariadata$latitude

#load required packages
library(rgeos)

# Convert the sf object to SpatialPolygons
africa_sp <- as(africa_shp, "Spatial")

# Simplify using gSimplify
africa_sp_simp <- gSimplify(africa_sp, tol = 0.05, topologyPreserve = TRUE)

# Convert back to sf object (optional)
africa_shp_simp <- st_as_sf(africa_sp_simp)

# Plot
plot(africa_shp_simp)


# Assuming sample_locations is loaded from the .RData file and is a data frame
# If it's not a data frame, you might need to http://127.0.0.1:25845/graphics/plot_zoom_png?width=1536&height=814adjust this section accordingly

# Create the data frame and SpatialPointsDataFrame
mydf <- structure(list(longitude = Malariadata$longitude, 
                       latitude = Malariadata$latitude), 
                  .Names = c("longitude", "latitude"), 
                  class = "data.frame", 
                  row.names = c(NA, nrow(Malariadata)))

xy <- mydf[,c(1,2)]
spdf <- SpatialPointsDataFrame(coords = xy, 
                               data = mydf, 
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Plot the sample locations on the map
plot(spdf, col = "red", add = TRUE)

# Convert BF (sf object) to SpatialPolygonsDataFrame (sp object)
DRC_sp <- as(DRC, "Spatial")

library(sf)
#Plot a polygon
xcoords<- c(14,14,15,15,15)
ycoords <- c(4,5,5,4,4)
poly1 <- sp::Polygon(cbind(xcoords,ycoords))
firstPoly <- sp::Polygons(list(poly1), ID = "A")
firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
CRS_new<-CRS(SRS_string = "EPSG:4326")
proj4string(firstSpatialPoly) <- CRS_new
plot(firstSpatialPoly)

library(sf)

# Read the Africa shapefile
africa1 <- st_read(dsn = "C:/Users/Administrator/Documents/Africa Map Visualization/Africa.shp")

# Read the Africa+ shapefile
africa2 <- st_read(dsn = "C:/Users/Administrator/Documents/Africa Map Visualization/Africa+.shp")

# Identify columns unique to each dataset
cols_africa1 <- setdiff(names(africa1), names(africa2))
cols_africa2 <- setdiff(names(africa2), names(africa1))

# Add columns from africa1 to africa2 with NA values
for (col in cols_africa1) {
  africa2[[col]] <- NA
}

# Add columns from africa2 to africa1 with NA values
for (col in cols_africa2) {
  africa1[[col]] <- NA
}

# Ensure the column order is the same before using rbind
africa2 <- africa2[, names(africa1)]

# Use rbind to combine the datasets
africa_shp <- rbind(africa1, africa2)


# Plot DRC
library(sf)
DRC <- africa_shp[which(africa_shp$country_na == 'Democratic Republic of the Congo'),]
plot(st_geometry(DRC), col="red")
africa_shp$country_na

# Install and load required packages
install.packages(c("INLA", "cli", "rlang", "PrevMap", "ggplot2", "scales", "raster"))
library(INLA)
library(cli)
library(rlang)
library(PrevMap)
library(ggplot2)
library(scales)
library(raster)

# Load the .RData file
Malariadata <- inputs$data_all_wa_ea_drc
print(Malariadata)
Malariadata$latitude

#load required packages
library(rgeos)

# Simplify the Africa shape using gSimplify
africa_shp_simp <- gSimplify(africa_shp, tol = 0.05)
plot(africa_shp_simp)

# Assuming sample_locations is loaded from the .RData file and is a data frame
# If it's not a data frame, you might need to http://127.0.0.1:25845/graphics/plot_zoom_png?width=1536&height=814adjust this section accordingly

# Create the data frame and SpatialPointsDataFrame
mydf <- structure(list(longitude = Malariadata$longitude, 
                       latitude = Malariadata$latitude), 
                  .Names = c("longitude", "latitude"), 
                  class = "data.frame", 
                  row.names = c(NA, nrow(Malariadata)))

xy <- mydf[,c(1,2)]
spdf <- SpatialPointsDataFrame(coords = xy, 
                               data = mydf, 
                               proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# Plot the sample locations on the map
plot(spdf, col = "red", add = TRUE)

# Convert DRC (sf object) to SpatialPolygonsDataFrame (sp object)
DRC_sp <- as(DRC, "Spatial")

# Set the CRS of BF_sp to match spdf
library(sp)
proj4string(DRC_sp) <- proj4string(spdf)
# Use over to find intersections
spdf_DRC <- over(spdf, DRC_sp)

# Filter points within DRC
DRC_points <- spdf[!is.na(spdf_DRC$country_na),]

# Plot DRC and the sample locations
plot(DRC_sp)
points(DRC_points, col="red", pch=20)

