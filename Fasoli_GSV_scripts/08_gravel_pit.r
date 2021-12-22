#!/usr/bin/env Rscript

#JCL setwd("/uufs/chpc.utah.edu/common/home/lin-group8/btf/stilt-sims/google-street-view/analysis")
setwd("./")

options(dplyr.summarise.inform = FALSE)

library(arrow)
library(glue)
library(knitr)
library(leaflet)
library(parallel)
library(raster)
library(tidyverse)


local_time <- function(time_string) {
    as.POSIXct(time_string, tz = "America/Denver")
}

seq_002 <- function(x, y) {
    seq(x, y, by = 0.002)
}


by_receptor <- read_feather("gsv-data/by_receptor/by_receptor.feather")
str(by_receptor)

bbox <- extent(-111.935, -111.901, 40.791, 40.841)


by_receptor <- by_receptor %>%
    mutate(
        ch4d_ppm_ex = ifelse(ch4d_ppm_n > 5, ch4d_ppm_ex, NA),
        co2d_ppm_ex = ifelse(co2d_ppm_n > 5, co2d_ppm_ex, NA),
        bc_ngm3_ex = ifelse(bc_ngm3_n > 5, bc_ngm3_ex, NA),
        co_ppb_ex = ifelse(co_ppb_n > 5, co_ppb_ex, NA),
        pm25_ugm3_ex = ifelse(pm25_ugm3_n > 5, pm25_ugm3_ex, NA)
    ) %>%
    filter(
        longitude >= bbox@xmin,
        longitude <= bbox@xmax,
        latitude >= bbox@ymin,
        latitude <= bbox@ymax
    )


simulation_ids <- by_receptor$uuid
filenames <- paste0(simulation_ids, "_foot.nc")

# Load and cache footprint rasters
filepaths <- file.path(
    "/uufs/chpc.utah.edu/common/home/lin-group8/btf/stilt-sims/google-street-view/stilt",
    "out", "footprints", filenames
)
rasters <- filepaths[file.exists(filepaths)] %>%
    mclapply(function(f) sum(brick(f)), mc.cores = 100) %>%
    brick() %>%
    readAll() %>%
    setZ(simulation_ids[file.exists(filepaths)], name = "job_id")
saveRDS(rasters, "data/gravelpit_footprint_cache.rds")
rasters <- readRDS("data/gravelpit_footprint_cache.rds")


source_locations <- bind_rows(list(
    expand.grid(
        longitude = seq_002(-111.913, -111.911),
        latitude = seq_002(40.819, 40.831)
    ),
    expand.grid(
        longitude = seq_002(-111.915, -111.911),
        latitude = seq_002(40.815, 40.817)
    ),
    expand.grid(
        longitude = seq_002(-111.917, -111.911),
        latitude = seq_002(40.809, 40.813)
    ),
    expand.grid(
        longitude = seq_002(-111.915, -111.907),
        latitude = 40.807
    ),
    expand.grid(
        longitude = seq_002(-111.913, -111.907),
        latitude = 40.805
    ),
    expand.grid(
        longitude = seq_002(-111.911, -111.907),
        latitude = 40.803
    ),
    data.frame(
        longitude = -111.909,
        latitude = 40.801
    )
))
write.csv(source_locations, "gravelpit_source_locations.csv", row.names = F)


source_cell_index <- cellFromXY(rasters, source_locations[c("longitude", "latitude")])
source_footprint_totals <- data.frame(
    job_id = getZ(rasters),
    gravelpit_footprint = colSums(rasters[source_cell_index])
)

#JCL(211119): bug--can't find "job_id" column in by_receptor
colnames(by_receptor)[colnames(by_receptor)=="uuid"] <- "job_id"

gridded_subdomain_with_footprint <- by_receptor %>%
    inner_join(source_footprint_totals, by = "job_id")
write.csv(gridded_subdomain_with_footprint,
    "gridded_subdomain_with_footprint.csv",
    row.names = F
)


# with(gridded_subdomain_with_footprint,
#      cor(pm25_ugm3_ex, gravelpit_footprint, use = 'pairwise.complete.obs'))
#
# mod <- lm(pm25_ugm3_ex ~ gravelpit_footprint, data = gridded_subdomain_with_footprint)
# summary(mod)

# gridded_subdomain_with_footprint <- gridded_subdomain_with_footprint %>%
#   mutate(pm25_co2 = (pm25_ugm3_ex + 1) / (co2d_ppm_ex + 1))# %>%
# # filter(!is.na(pm25_co2), pm25_co2 < 1e8, pm25_co2 > 1e-8)
# with(gridded_subdomain_with_footprint,
#      cor(pm25_co2, gravelpit_footprint, use = 'pairwise.complete.obs'))
#
# mod <- lm(pm25_co2 ~ gravelpit_footprint, data = gridded_subdomain_with_footprint)
# summary(mod)


# bbox <- extent(-111.935, -111.901, 40.791, 40.841)

# gravel_pit_receptors <- expand.grid(
#   longitude = seq(bbox@xmin, bbox@xmax, by = 0.002),
#   latitude = seq(bbox@ymin, bbox@ymax, by = 0.002)
# )

gravel_pit_receptors <- unique(by_receptor[c("longitude", "latitude")])

write.csv(gravel_pit_receptors, "gravel_pit_receptor_locations_002.csv", row.names = F)


gridded_subdomain_with_footprint_case_studies <- gridded_subdomain_with_footprint %>%
    filter(((time >= local_time("2019-07-25")) & (time <= local_time("2019-07-26"))) |
        ((time >= local_time("2019-08-08")) & (time <= local_time("2019-08-09"))) |
        ((time >= local_time("2019-08-22")) & (time <= local_time("2019-08-23"))) |
        ((time >= local_time("2019-12-10")) & (time <= local_time("2019-12-11"))) |
        ((time >= local_time("2020-01-29")) & (time <= local_time("2020-01-30")))) # %>%
# filter(gravelpit_footprint > 1e-2)

mod <- lm(pm25_ugm3_ex ~ gravelpit_footprint,
    data = gridded_subdomain_with_footprint_case_studies
)
summary(mod)


case_study_dates <- c(
    local_time("2019-07-25"),
    local_time("2019-08-08"),
    local_time("2019-08-22"),
    local_time("2019-12-10"),
    local_time("2020-01-29")
)

i <- 2
# i <- i + 1
date <- case_study_dates[i]
print(date)
gridded_subdomain_with_footprint %>%
    filter(
        time >= date,
        time <= date + 86400
    ) %>%
    filter(pm25_ugm3_n > 5) %>%
    with(cor(pm25_ugm3_ex, gravelpit_footprint))
# ggplot(aes(x = longitude, y = latitude, color = pm25_ugm3_ex)) +
# geom_point() +
# scale_color_viridis_c() +
# theme_classic() +
# labs(x = NULL,
#      y = NULL,
#      color = NULL)

date <- local_time("2019-08-08")
drive_pass_locations <- bind_rows(list(
    data.frame(longitude = -111.905, latitude = 40.795),
    data.frame(longitude = -111.907, latitude = 40.797),
    data.frame(longitude = -111.909, latitude = 40.799),
    data.frame(longitude = -111.911, latitude = 40.799),
    data.frame(longitude = -111.911, latitude = 40.801),
    data.frame(longitude = -111.913, latitude = 40.801),
    data.frame(longitude = -111.913, latitude = 40.803),
    data.frame(longitude = -111.915, latitude = 40.803),
    data.frame(longitude = -111.917, latitude = 40.805),
    data.frame(longitude = -111.917, latitude = 40.807),
    data.frame(longitude = -111.919, latitude = 40.807),
    data.frame(longitude = -111.919, latitude = 40.809),
    data.frame(longitude = -111.919, latitude = 40.811),
    data.frame(longitude = -111.919, latitude = 40.813),
    data.frame(longitude = -111.919, latitude = 40.815),
    data.frame(longitude = -111.917, latitude = 40.817),
    data.frame(longitude = -111.915, latitude = 40.821),
    data.frame(longitude = -111.915, latitude = 40.823),
    data.frame(longitude = -111.915, latitude = 40.825),
    data.frame(longitude = -111.915, latitude = 40.827),
    data.frame(longitude = -111.915, latitude = 40.829),
    data.frame(longitude = -111.915, latitude = 40.831),
    data.frame(longitude = -111.915, latitude = 40.833),
    data.frame(longitude = -111.913, latitude = 40.837),
    data.frame(longitude = -111.911, latitude = 40.839)
))

data_on_day <- gridded_subdomain_with_footprint %>%
    filter(
        time >= date,
        time <= date + 86400
    ) %>%
    filter(pm25_ugm3_n > 5)
# inner_join(drive_pass_locations, by = c('longitude', 'latitude'))

data_on_day %>%
    select(latitude, longitude, time, pm25_ugm3_ex) %>%
    write.csv("gravelpit_aug8_pm25.csv", row.names = F)

mod <- lm(pm25_ugm3_ex ~ gravelpit_footprint, data = data_on_day)
summary(mod)
with(
    data_on_day,
    plot(gravelpit_footprint, pm25_ugm3_ex)
)
abline(mod, col = "red")


max_measurement_index <- which.max(data_on_day$pm25_ugm3_ex)
job_id <- data_on_day$job_id[max_measurement_index]

simulation_ids <- getZ(rasters)
job_index <- which(simulation_ids == job_id)


latitude <- data_on_day$latitude[i]
longitude <- data_on_day$longitude[i]

r <- rasters[[job_index]]
# r <- sum(brick(paste0('footprints/gravelpit/', job_id, '.nc')))

# leaflet() %>%
#   addProviderTiles('CartoDB.Positron') %>%
#   addRasterImage(r, opacity = 0.5) %>%
#   addMarkers(lng = longitude, lat = latitude)

r %>%
    as.data.frame(xy = T) %>%
    setNames(c("longitude", "latitude", "footprint")) %>%
    filter(footprint > 1e-8) %>%
    mutate(footprint = format(footprint, scientific = FALSE)) %>%
    write.csv("gravelpit_aug8_footprint_example.csv", row.names = F, quote = F)

with(
    gridded_subdomain_with_footprint_case_studies,
    plot(gravelpit_footprint, pm25_ugm3_ex)
)
abline(mod, col = "red")

quantile(gridded_subdomain_with_footprint_case_studies$pm25_ugm3_ex, na.rm = T)