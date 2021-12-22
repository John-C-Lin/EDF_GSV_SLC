library(arrow)
library(tidyverse)
library(rgdal)

by_polygon_drive_day <- readRDS("data/receptors/by_polygon_drive_day.rds")
receptors_with_sectors <- readRDS("data/receptors/receptors_with_sectors.rds")


by_receptor <- receptors_with_sectors %>%
    rename(
        airport_ppm_ex = airport,
        cement_ppm_ex = cement,
        cmv_ppm_ex = cmv,
        commercial_ppm_ex = commercial,
        elec_prod_ppm_ex = elec_prod,
        industrial_ppm_ex = industrial,
        nonroad_ppm_ex = nonroad,
        onroad_ppm_ex = onroad,
        rail_ppm_ex = rail,
        residential_ppm_ex = residential,
        total_ppm_ex = total
    )

yyyy_mm <- strftime(by_receptor$time, "%Y-%m")
for (id in unique(yyyy_mm)) {
    mask <- yyyy_mm == id
    filename <- file.path("gsv-data/by_receptor/by_month", paste0(id, ".csv"))
    dir.create(dirname(filename), showWarnings = FALSE)``
    write_csv(by_receptor[mask, ], filename)
}
write_feather(by_receptor, "gsv-data/by_receptor/by_receptor.feather", compression = "zstd", compression_level = 9)



socio <- read_csv("data/socioeconomic/edf_polygons_v20191124.dat") %>%
    setNames(tolower(colnames(.))) %>%
    rename(
        polygon_id = polygon_name,
        population = pop,
        population_density = pop_density,
        race_minority_pct = race_pct_minority
    ) %>%
    mutate(
        race_white_pct = 100 * race_white / race_total,
        race_hispanic_pct = 100 * race_hispanic / race_total,
        race_black_pct = 100 * race_black / race_total,
        race_am_indian_pct = 100 * race_am_indian / race_total,
        race_asian_pct = 100 * race_asian / race_total,
        race_pac_isl_pct = 100 * race_pac_isl / race_total,
        race_other_pct = 100 * race_other / race_total,
        race_multi_pct = 100 * race_multiracial / race_total,
    ) %>%
    # The calc_pct_minority here is the same as the provided race_pct_minority
    # mutate(
    #     minority = race_hispanic + race_black + race_am_indian + race_asian + race_pac_isl + race_other + race_multiracial,
    #     calc_pct_minority = 100 * (minority / race_total)
    # ) %>%
    filter(race_total > 0) %>%
    select(
        polygon_id,
        area_sq_km,
        population,
        population_density,
        income,
        race_white_pct,
        race_hispanic_pct,
        race_black_pct,
        race_am_indian_pct,
        race_asian_pct,
        race_pac_isl_pct,
        race_other_pct,
        race_multi_pct,
        race_minority_pct
    )

by_polygon_drive_day <- by_polygon_drive_day %>%
    full_join(socio, by = "polygon_id")
write_feather(by_polygon_drive_day, "gsv-data/by_polygon_drive_day/by_polygon_drive_day.feather")
write_csv(by_polygon_drive_day, "gsv-data/by_polygon_drive_day/by_polygon_drive_day.csv")


by_polygon <- by_polygon_drive_day %>%
    group_by(polygon_id) %>%
    summarize(across(ends_with("_ex"), median, na.rm = T))

by_polygon <- full_join(by_polygon, socio, by = "polygon_id")
write_feather(by_polygon, "gsv-data/by_polygon/by_polygon.feather")
write_csv(by_polygon, "gsv-data/by_polygon/by_polygon.csv")


polygons <- readOGR("data/polygons/Priority_1_v2019-11-24.kml", stringsAsFactors = FALSE)
colnames(polygons@data) <- c("polygon_id", "polygon_description")
polygons@data <- full_join(polygons@data, by_polygon, by = "polygon_id")
polygons <- polygons[!is.na(polygons$total_ppm_ex), ]
writeOGR(polygons, "gsv-data/by_polygon/by_polygon.geojson", layer = "by-polygon", driver = "GeoJSON", overwrite = TRUE)