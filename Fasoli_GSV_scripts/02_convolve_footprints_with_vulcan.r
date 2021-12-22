#!/usr/bin/env Rscript

library(arrow)
library(dplyr)
library(ncdf4)
library(parallel)
library(raster)


rasterOptions(chunksize = 1e9, maxmemory = 1e11)


sector_names <- c(
    "airport", "cement", "cmv", "commercial", "elec_prod", "industrial",
    "nonroad", "onroad", "rail", "residential", "total"
)


get_footprint_filename_by_uuid <- function(uuid) {
    file.path(
        "/uufs/chpc.utah.edu/common/home/",
        "lin-group8/btf/stilt-sims/",
        "google-street-view/stilt/out/footprints_01",
        uuid, paste0(uuid, "_foot.nc")
    )
}


get_vulcan_sector <- function(sector) {
    filepath <- file.path(
        "/uufs/chpc.utah.edu/common/home/",
        "lin-group9/inventories/vulcan/v3/",
        "data/rds/gsv/2015"
    )
    filename <- file.path(filepath, paste0(sector, ".rds"))
    readRDS(filename)
}


get_vulcan_sectors <- function() {
    vulcan <- mclapply(sector_names,
        get_vulcan_sector,
        mc.cores = length(sector_names)
    )
    names(vulcan) <- sector_names
    vulcan
}


get_enhancement_by_uuid <- function(uuid, vulcan, vulcan_time) {
    footprint_filename <- get_footprint_filename_by_uuid(uuid)

    if (!file.exists(footprint_filename)) {
        message("Missing ", footprint_filename)
        return(NA)
    }

    footprint_raster <- raster(footprint_filename)

    footprint_time <- getZ(footprint_raster)
    if (nchar(footprint_time) < 11) {
        footprint_time <- paste(footprint_time, "00")
    }
    footprint_time <- as.POSIXct(footprint_time, tz = "UTC", format = "%Y-%m-%d %H")

    time_index <- match(footprint_time, vulcan_time)
    mclapply(
        vulcan,
        function(sector) {
            product <- resample(footprint_raster, sector) * subset(sector, time_index)
            sum(values(product), na.rm = T)
        },
        mc.cores = length(vulcan)
    )
}


is_leap_year <- function(year) {
    is_leap <- (year %% 4) == 0
    is_leap <- is_leap & ((year %% 100 != 0))
    is_leap <- is_leap | ((year %% 400) == 0)
    is_leap
}


scale_vulcan <- function(vulcan, start, stop) {
    vulcan_time <- getZ(vulcan[[1]])

    start <- as.POSIXct(trunc(start, units = "hours"))
    stop <- as.POSIXct(trunc(stop, units = "hours")) + 3600

    start_year <- as.numeric(strftime(start, tz = "UTC", "%Y"))
    stop_year <- as.numeric(strftime(stop, tz = "UTC", "%Y"))

    start_day_index <- as.numeric(strftime(start, tz = "UTC", "%j"))
    stop_day_index <- as.numeric(strftime(stop, tz = "UTC", "%j"))
    stop_day_index <- ifelse(is_leap_year(stop_year), stop_day_index - 1, stop_day_index)

    vulcan_start_day_of_week <- as.numeric(strftime(
        min(vulcan_time),
        tz = "UTC", format = "%u"
    )) # thu, 4

    layer_order <- c()
    layer_time <- as.POSIXct(character())
    for (year in start_year:stop_year) {
        start_day_of_week <- as.numeric(strftime(
            as.POSIXct(paste0(start_year, "-01-01")),
            format = "%u"
        ))

        shift_days <- vulcan_start_day_of_week - start_day_of_week
        shift_order <- c((shift_days + 1):nlayers(vulcan[[1]]), 1:shift_days)

        if (is_leap_year(year)) {
            leap_day_of_week <- as.numeric(strftime(
                as.POSIXct(paste0(year, "-02-29")),
                format = "%u"
            ))

            if (leap_day_of_week %in% c(1, 6)) {
                # Look forward if Feb 29 falls on Mon or Sat
                repeat_index <- 60 # Feb 29 = Mar 1
            } else if (leap_day_of_week %in% c(2:5, 7)) {
                # Look backward if Feb 29 falls Tue-Fri or Sun
                repeat_index <- 59 # Feb 29 = Feb 28
            }

            shift_order <- c(shift_order[1:59], shift_order[59:length(shift_order)])
        }

        shift_time <- as.POSIXct(paste0(year, "-01-01"), tz = "UTC") + (0:(length(shift_order) - 1)) * 3600

        layer_order <- c(layer_order, shift_order)
        layer_time <- c(layer_time, shift_time)
    }

    layer_index <- (start_day_index * 24 - 23):(length(layer_order) - ((365 - stop_day_index) * 24 - 23))
    layer_order <- layer_order[layer_index]
    layer_time <- layer_time[layer_index]

    mclapply(
        vulcan,
        function(sector) {
            sector <- subset(sector, layer_order)
            setZ(sector, layer_time)
        },
        mc.cores = length(vulcan)
    )
}


receptors <- read_feather("data/receptors/receptors.feather")

start <- min(receptors$time)
stop <- max(receptors$time)


# vulcan <- get_vulcan_sectors()
# vulcan <- scale_vulcan(vulcan, start, stop)
# vulcan_time <- getZ(vulcan[[1]])
# enhancement_list <- mclapply(receptors$uuid,
#                              get_enhancement_by_uuid,
#                              vulcan = vulcan,
#                              vulcan_time = vulcan_time,
#                              mc.cores = 104)
# saveRDS(enhancement_list, 'enhancement_list.rds')





enhancement_list <- readRDS("data/enhancement_list.rds")

is_error <- sapply(enhancement_list, function(enhancement_sectors) {
    any(
        sapply(enhancement_sectors, function(enhancement_sector) {
            class(enhancement_sector) == "try-error"
        })
    ) || any(is.na(enhancement_sectors))
})

error_index <- which(is_error)

# Spot check errored simulations
# uuid <- receptors$uuid[error_index[1]]
# get_enhancement_by_uuid(uuid, vulcan = vulcan, vulcan_time = vulcan_time)


na_sectors <- rep(NA, length(sector_names))
names(na_sectors) <- sector_names
list(na_sectors)
enhancement_list[is_error] <- list(na_sectors)


enhancement_df <- bind_rows(enhancement_list)
with_enhancements <- cbind(receptors, enhancement_df)
str(with_enhancements)
saveRDS(with_enhancements, "data/receptors/receptors_with_sectors.rds")