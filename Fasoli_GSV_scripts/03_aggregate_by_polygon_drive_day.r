#!/usr/bin/env Rscript

library(plotly)
library(tidyverse)

df <- readRDS("data/receptors/receptors_with_sectors.rds")
df <- df[!is.na(df$polygon_id), ]

# Cast to local time since aggregations are performed by daily drive (midnight - midnight)
attributes(df$time)$tzone <- "America/Denver"


minutes_in_polygon <- df %>%
    mutate(time = as.POSIXct(trunc(time, units = "mins"))) %>%
    (function(group) unique(group[c("time", "system_id", "polygon_id")])) %>%
    mutate(time = as.POSIXct(trunc(time, units = "days"))) %>%
    group_by(time, system_id, polygon_id) %>%
    summarize(minutes_in_polygon = n()) %>%
    ungroup()


# by_drive_pass <- df %>%
#     mutate(time = as.POSIXct(trunc(time, units = "days"))) %>%
#     group_by(time, system_id, polygon_id) %>%
#     filter(n() >= 60) %>%
by_drive_pass <- df %>%
    mutate(time = as.POSIXct(trunc(time, units = "days"))) %>%
    full_join(minutes_in_polygon, by = c("time", "system_id", "polygon_id")) %>%
    filter(minutes_in_polygon >= 60) %>%
    group_by(time, system_id, polygon_id) %>%
    summarize(
        ch4d_ppm_ex = sum(ch4d_ppm_ex * ch4d_ppm_n, na.rm = T) / sum(ch4d_ppm_n),
        co2d_ppm_ex = sum(co2d_ppm_ex * co2d_ppm_n, na.rm = T) / sum(co2d_ppm_n),
        bc_ngm3_ex = sum(bc_ngm3_ex * bc_ngm3_n, na.rm = T) / sum(bc_ngm3_n),
        co_ppb_ex = sum(co_ppb_ex * co_ppb_n, na.rm = T) / sum(co_ppb_n),
        no_ppb_ex = sum(no_ppb_ex * no_ppb_n, na.rm = T) / sum(no_ppb_n),
        no2_ppb_ex = sum(no2_ppb_ex * no2_ppb_n, na.rm = T) / sum(no2_ppb_n),
        pm25_ugm3_ex = sum(pm25_ugm3_ex * pm25_ugm3_n, na.rm = T) / sum(pm25_ugm3_n),
        airport_ppm_ex = sum(airport * n, na.rm = T) / sum(n),
        cement_ppm_ex = sum(cement * n, na.rm = T) / sum(n),
        cmv_ppm_ex = sum(cmv * n, na.rm = T) / sum(n),
        commercial_ppm_ex = sum(commercial * n, na.rm = T) / sum(n),
        elec_prod_ppm_ex = sum(elec_prod * n, na.rm = T) / sum(n),
        industrial_ppm_ex = sum(industrial * n, na.rm = T) / sum(n),
        nonroad_ppm_ex = sum(nonroad * n, na.rm = T) / sum(n),
        onroad_ppm_ex = sum(onroad * n, na.rm = T) / sum(n),
        rail_ppm_ex = sum(rail * n, na.rm = T) / sum(n),
        residential_ppm_ex = sum(residential * n, na.rm = T) / sum(n),
        total_ppm_ex = sum(total * n, na.rm = T) / sum(n)
    ) %>%
    ungroup()


# all cement and cmv values == 0 so remove
by_drive_pass <- select(by_drive_pass, -cement_ppm_ex, -cmv_ppm_ex)
saveRDS(by_drive_pass, "data/receptors/by_polygon_drive_day.rds")

# cor.mtest <- function(mat, ...) {
#     mat <- as.matrix(mat)
#     n <- ncol(mat)
#     p.mat<- matrix(NA, n, n)
#     diag(p.mat) <- 0
#     for (i in 1:(n - 1)) {
#         for (j in (i + 1):n) {
#             tmp <- cor.test(mat[, i], mat[, j], ...)
#             p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
#         }
#     }
#   colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
#   p.mat
# }

# observations <- select(by_drive_pass, -system_id, -polygon_id, -time)
# correlation_matrix <- cor(observations, use = 'pairwise.complete.obs')
# correlation_matrix_p <- cor.mtest(correlation_matrix)

# png('correlation_matrix_p.png', width = 1500, height = 1500, res = 144)
# corrplot(
#     correlation_matrix,
#     type = 'upper',
#     method = 'color',
#     addCoef.col = 'black',
#     tl.col = 'black',
#     p.mat = correlation_matrix_p,
#     sig.level = 0.05,
# )
# dev.off()