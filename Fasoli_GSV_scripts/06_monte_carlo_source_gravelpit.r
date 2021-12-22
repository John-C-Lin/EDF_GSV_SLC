# Dec. 14th, 2021 by John C. Lin (John.Lin@utah.edu)
# Adapt Ben Fasoli's "06_monte_carlo_source.r" from landfill/CH4 case to Gravel Pit case

#!/usr/bin/env Rscript
rm(list = ls())
gc()

##############
#tracer <- "nox_ppb_ex"
tracer <- "pm25_ugm3_ex"
R.threshold <- 0.2   # threshold in correlation coefficient to be considered in "source-optiimized" flux calculation
regenerate.footTF <- FALSE   # whether or not to regenerate all footprints (TRUE), or load from previous run (FALSE)
##############
#JCL:  setwd("/uufs/chpc.utah.edu/common/home/lin-group8/btf/stilt-sims/google-street-view/analysis")
setwd("/uufs/chpc.utah.edu/common/home/u0791084/PROJECTS/Lin_EDF_GSV_SLC/Fasoli_GSV_scripts")
options(dplyr.summarise.inform = FALSE)

library(arrow); library(leaflet); library(mapview)
library(ncdf4); library(parallel); library(raster); library(tidyverse)

rasterOptions(chunksize = 1e11, maxmemory = 1e11)


get_footprint_filename_by_uuid <- function(uuid) {
    file.path(
        "/uufs/chpc.utah.edu/common/home/",
        "lin-group8/btf/stilt-sims/",
        "google-street-view/stilt/out/by-id",
        uuid, paste0(uuid, "_foot.nc")
    )
}


load_footprint_raster <- function(filename) {
    readAll(raster(filename))
}



# Define source and observation regions ----------------------------------------
# small source area over landfill
# source_bbox <- extent(-112.053, -112.035, 40.743, 40.749)
# medium source area containing landfill and surrounding area
# source_bbox <- extent(-112.07, -112.028, 40.742, 40.768)
# JCL:  subset of Gravel Pit area near Beck Street as a BOX
#source_bbox <- extent(-111.92, -111.9035, 40.794, 40.81)
# JCL:  subset of Gravel Pit area near Beck Street as a POLYGON
#x_coords <- c(-111.92, -111.9035, -111.9035, -111.92, -111.92)  # test with rectangle first
#y_coords <- c(40.794, 40.794, 40.81, 40.81, 40.794)             # test with rectangle first
#  Gravel Pit potential source coordinates as POLYGONS
x_coords <- c(-111.9151162369536, -111.9173277749282, -111.9221653805507, -111.9212393539757, -111.9148133173202, -111.9000109387472, -111.9086904007313, -111.9076404616799, -111.9111505165306, -111.9151162369536)
y_coords <- c(40.83139914238476, 40.81837534323808, 40.81261585646753, 40.80431024481059, 40.79126112721225, 40.79499345590276, 40.81794770917191, 40.82293997609, 40.83183646194003, 40.83139914238476) 

poly1 <- sp::Polygon(cbind(x_coords,y_coords))
firstPoly <- sp::Polygons(list(poly1), ID = "A")
source_bbox <- sp::SpatialPolygons(list(firstPoly))


# Line observations during case study day 2019-08-08
# observations_bbox <- extent(-112.075, -112.025, 40.741 - 1e-8, 40.741 + 1e-8)
# observations <- readRDS('data/receptors.rds')   %>%
#   filter(time > as.POSIXct('2019-08-08', tz = 'America/Denver'),
#          time < as.POSIXct('2019-08-09', tz = 'America/Denver')) %>%
#   filter(longitude >= observations_bbox@xmin,
#          longitude <= observations_bbox@xmax,
#          latitude >= observations_bbox@ymin,
#          latitude <= observations_bbox@ymax)


# Line observations for entire methane record
# observations_bbox <- extent(-112.09, -112.025, 40.741 - 0.005, 40.741 + 0.005)
# observations <- read_feather("gsv-data/by_receptor/by_receptor.feather") %>%
#     filter(
#         !is.na(ch4d_ppm_ex),
#         longitude >= observations_bbox@xmin,
#         longitude <= observations_bbox@xmax,
#         latitude >= observations_bbox@ymin,
#         latitude <= observations_bbox@ymax
#     )

# Area observations for entire methane record
# observations_bbox <- extent(-112.09, -112.01, 40.715, 40.775)
# observations <- read_feather("gsv-data/by_receptor/by_receptor.feather") %>%
#     filter(
#         !is.na(ch4d_ppm_ex),
#         longitude >= observations_bbox@xmin,
#         longitude <= observations_bbox@xmax,
#         latitude >= observations_bbox@ymin,
#         latitude <= observations_bbox@ymax
#     )

# JCL:  observations covering Gravel Pit area along Beck Street
observations_bbox <- extent(-111.926, -111.9035, 40.794, 40.83)
observations <- read_feather("gsv-data/by_receptor/by_receptor.feather") %>%
    filter(
        longitude >= observations_bbox@xmin,
        longitude <= observations_bbox@xmax,
        latitude >= observations_bbox@ymin,
        latitude <= observations_bbox@ymax,
    )
nox_ppb_ex <- observations$no_ppb_ex + observations$no2_ppb_ex
observations <- data.frame(observations, nox_ppb_ex)

filenames <- get_footprint_filename_by_uuid(observations$uuid)
if(regenerate.footTF){
  footprint_raster_list <- mclapply(filenames, load_footprint_raster, mc.cores = 500)
  footprint_rasters <- readAll(brick(footprint_raster_list))
  footprint_full_domain <- footprint_rasters
  saveRDS(footprint_full_domain,"footprint_full_domain_gravel.rds")
} else {
  print("reading in footprint_full_domain_gravel.rds from previous run")
  footprint_full_domain <- readRDS("footprint_full_domain_gravel.rds")
} # if(regenerate.footTF){

# filter footprint and observations by selected tracer
footprint_full_domain <- subset(footprint_full_domain,subset=which(!is.na(observations[,tracer])))
observations <- observations[!is.na(observations[,tracer]),]

observation_count <- nrow(observations)
measured <- observations[,tracer]

footprint_rasters <- crop(footprint_full_domain, source_bbox)
footprints <- as.array(footprint_rasters)
# JCL:
saveRDS(footprints,"footprints_gravel.rds")


# Randomized source field ------------------------------------------------------
random_iteration_count <- 1000
source_dim <- c(dim(footprints)[1:2], random_iteration_count)
sources <- array(runif(prod(source_dim)), dim = source_dim)


stilt_weighted_sources <- array(dim = c(dim(footprints), random_iteration_count))
predictions <- array(dim = c(observation_count, random_iteration_count))
r <- array(dim = random_iteration_count)

for (i in 1:random_iteration_count) {
    print(paste("Monte Carlo simulation",i,"out of",random_iteration_count,"realizations"))
    source <- sources[, , i]

    # JCL: footprint multiplied with randomized source 
    stilt_weighted_source <- footprints * c(source)
    # JCL: vector of predicted [tracer] across all receptors
    prediction <- apply(stilt_weighted_source, 3, sum)

    stilt_weighted_sources[, , , i] <- stilt_weighted_source
    predictions[, i] <- prediction
    # JCL: correlation between vector of measured vs predicted [tracer]
    r[i] <- cor(measured, prediction)
}

# total sensitivity of observations to each source location and value
# source_influence <- apply(stilt_weighted_sources, c(1, 2, 4), sum)
# source_r <- apply(source_influence, c(1, 2), function(x) cor(x, r))
# JCL: correlation between randomized emissions at a particular gridcell
#      and the vector of correlation coefficients over the # of random iterations
#      presumably, if a gridcell is a "hotspot", then there should be positive correlation
#      between that gridcell and the vector of correlation coefficients
source_r <- apply(sources, c(1, 2), function(x) cor(x, r))
saveRDS(source_r,"source_r.rds")

png(paste0("img/monte_carlo/source_r_gravel_",tracer,".png"), width = 2000, height = 2000, res = 300)
require("fields")
image.plot(source_r,main=tracer)
dev.off()

# Observation : Prediction correlation map -------------------------------------
image <- raster(source_r, template = footprint_rasters[[1]])
color_range <- c(-0.5, 0.5)
image[image < min(color_range)] <- min(color_range)
image[image > max(color_range)] <- max(color_range)

m <- leaflet() %>%
    # addProviderTiles('CartoDB.Positron') %>%
    addProviderTiles("Esri.WorldImagery") %>%
    addCircles(
        lng = observations$longitude,
        lat = observations$latitude,
        radius = 50,
        fillOpacity = 0.5,
        weight = 1,
        color = colorNumeric("Greys",
            domain = c(measured)
        )(measured),
        popup = as.character(measured)
    ) %>%
    # addRasterImage(log10(footprint_full_domain[[which.max(measured)]]),
    #                opacity = 0.75,
    #                colors = colorNumeric('Greys',
    #                                      domain = c(-4, 0),
    #                                      na.color = 'transparent')) %>%
    addRasterImage(image,
        opacity = 0.75,
        colors = colorNumeric("RdBu",
            domain = color_range,
            reverse = T,
            na.color = "transparent"
        )
    ) %>%
    addLegend(
        position = "bottomleft",
        pal = colorNumeric("RdBu",
            domain = color_range,
            reverse = T,
            na.color = "transparent"
        ),
        values = color_range
    )
m
# mapshot(m, file = "source-figures/correlation_map.png")

# Randomized source maps -------------------------------------------------------
# for (i in 1:4) {
#   source_raster <- raster(sources[,,i], template = footprint_rasters)
#   colors <- colorNumeric('Spectral',
#                          domain = c(0, 1),
#                          reverse = T)
#   m <- leaflet() %>%
#     addProviderTiles('Esri.WorldImagery') %>%
#     addRasterImage(source_raster,
#                    opacity = 0.5,
#                    colors = colors) %>%
#     addLegend(position = 'bottomleft',
#               pal = colors,
#               values = c(0, 1))
#   mapshot(m, file = paste0('source-figures/source_', i, '.png'))
# }
#
# for (i in 1:4) {
#   source_raster <- raster(source_influence[,,i], template = footprint_rasters)
#   colors <- colorNumeric('Spectral',
#                          domain = c(0, .1),
#                          reverse = T)
#   m <- leaflet() %>%
#     addProviderTiles('Esri.WorldImagery') %>%
#     addRasterImage(source_raster,
#                    opacity = 0.5,
#                    colors = colors) %>%
#     addLegend(position = 'bottomleft',
#               pal = colors,
#               values = c(0, .1))
#   mapshot(m, file = paste0('source-figures/source_influence_', i, '.png'))
# }


# Distribution of Observation : Prediction correlations ------------------------
png("img/monte_carlo/density_gravel.png", width = 2000, height = 2000, res = 300)
h <- hist(r, breaks = 25)
h$density <- h$counts / sum(h$counts) * 100
plot(h,
    freq = F,
    col = "steelblue",
    main = "Observation : Prediction correlation",
    xlab = "Pearson correlation coefficient (r)",
    ylab = "Frequency  (% total count)"
)
dev.off()


# Observation : Prediction scatterplot -----------------------------------------
predicted <- apply(footprints, 3, sum)
mod <- lm(measured ~ predicted)
df <- data.frame(
    predicted = predicted,
    measured = measured
)
ggplot(df, aes(x = predicted, y = measured)) +
    geom_point() +
    geom_abline(
        slope = mod$coefficients[2],
        intercept = mod$coefficients[1],
        color = "red"
    ) +
    theme_classic() +
    labs(
        x = bquote("Source Normalized Footprint  " ~
        frac(ppm, frac(mu * mol, m^2 ~ s))),
        y = tracer,
        title = "Gravel Pit observations vs STILT",
        subtitle = bquote("Flux" == .(round(mod$coefficients[2], 1)) ~ frac(mu * mol, m^2 ~ s) ~
        "    " ~ r^2 == .(round(summary(mod)$r.squared, 3)))
    )
ggsave("img/monte_carlo/source-binary-scatterplot_gravel.png", width = 6, height = 5)


# Threshold in correlation coefficient defines the source region
# source_mask <- source_r > quantile(source_r, 0.99)
# source_mask <- source_r > 0.05
source_mask <- source_r > R.threshold
binary_source <- array(0, dim = dim(source_mask))
binary_source[source_mask] <- 1

predicted <- apply(footprints * c(binary_source), 3, sum)
# filter out times when footprint is 0
#sel <- predicted > 0
#predicted <- predicted[sel]
#measured <- measured[sel]
mod <- lm(measured ~ predicted)
summary(mod)

df <- data.frame(
    predicted = predicted,
    measured = measured
)
ggplot(df, aes(x = predicted, y = measured)) +
    geom_point() +
    geom_abline(
        slope = mod$coefficients[2],
        intercept = mod$coefficients[1],
        color = "red"
    ) +
    theme_classic() +
    labs(
        x = bquote("Source Normalized Footprint  " ~
        frac(ppm, frac(mu * mol, m^2 ~ s))),
        y = tracer,
        title = "Gravel Pit observations vs STILT",
        caption = paste("For source gridcells with R >",R.threshold),
        subtitle = bquote("Flux" == .(round(mod$coefficients[2], 1)) ~ frac(mu * mol, m^2 ~ s) ~
        "    " ~ r^2 == .(round(summary(mod)$r.squared, 3)))
    )
ggsave("img/monte_carlo/source-optimized-scatterplot_gravel.png", width = 6, height = 5)


# source_raster <- raster(binary_source, template = footprint_rasters)
# values(source_raster) <- 1
# m <- leaflet() %>%
#     addProviderTiles("Esri.WorldImagery") %>%
#     addRasterImage(source_raster,
#         opacity = 0.75,
#         color = "red"
#     )
# m
# # mapshot(m, file = "img/monte_carlo/source_binary_map.png")

# source_raster <- raster(binary_source, template = footprint_rasters)
# source_raster[source_raster == 0] <- NA
# m <- leaflet() %>%
#     addProviderTiles("Esri.WorldImagery") %>%
#     addRasterImage(source_raster,
#         opacity = 0.75,
#         colors = "red"
#     )
# m
# mapshot(m, file = "img/monte_carlo/source_optimized_map.png")