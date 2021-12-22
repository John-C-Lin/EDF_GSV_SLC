# library(corrplot)
# library(Hmisc)
# library(tidyverse)


# by_polygon <- readRDS('data/export/by_polygon.rds') %>%
#     mutate(nox_ppb_ex = no_ppb_ex + no2_ppb_ex) %>%
#     select(
#         co2d_ppm_ex,
#         pm25_ugm3_ex,
#         nox_ppb_ex,
#         no_ppb_ex,
#         no2_ppb_ex,
#         bc_ngm3_ex,
#         ch4d_ppm_ex,
#         everything()
#     )
# correlation <- select(by_polygon, -polygon_id) %>%
#     as.matrix() %>%
#     rcorr()
# r <- correlation$r
# p <- correlation$P
# png('img/correlation_matrix/by_polygon.png', width = 3000, height = 3000, res = 144)
# corrplot(
#     r,
#     type = 'upper',
#     method = 'color',
#     addCoef.col = 'black',
#     tl.col = 'black',
#     p.mat = p,
#     sig.level = 0.05,
# )
# dev.off()



# correlation <- by_polygon %>%
#     select(
#         co2d_ppm_ex,
#         pm25_ugm3_ex,
#         nox_ppb_ex,
#         no_ppb_ex,
#         no2_ppb_ex,
#         bc_ngm3_ex,
#         ch4d_ppm_ex,
#         vulcan_total_ppm_ex = total_ppm_ex,
#         income,
#         race_pct_white,
#         race_pct_minority
#     ) %>%
#     as.matrix() %>%
#     rcorr()
# r <- correlation$r
# p <- correlation$P
# png('img/correlation_matrix/by_polygon_subset.png', width = 1000, height = 1000, res = 144)
# corrplot(
#     r,
#     type = 'upper',
#     method = 'color',
#     addCoef.col = 'black',
#     tl.col = 'black',
#     p.mat = p,
#     sig.level = 0.05,
# )
# dev.off()



# by_polygon_drive_day <- readRDS('data/export/by_polygon_drive_day.rds') %>%
#     mutate(nox_ppb_ex = no_ppb_ex + no2_ppb_ex) %>%
#     select(
#         co2d_ppm_ex,
#         pm25_ugm3_ex,
#         nox_ppb_ex,
#         no_ppb_ex,
#         no2_ppb_ex,
#         bc_ngm3_ex,
#         ch4d_ppm_ex,
#         everything()
#     )
# correlation <- select(by_polygon_drive_day, -polygon_id, -system_id, -time) %>%
#     as.matrix() %>%
#     rcorr()
# r <- correlation$r
# p <- correlation$P
# png('img/correlation_matrix/by_polygon_drive_day.png', width = 3000, height = 3000, res = 144)
# corrplot(
#     r,
#     type = 'upper',
#     method = 'color',
#     addCoef.col = 'black',
#     tl.col = 'black',
#     p.mat = p,
#     sig.level = 0.05,
# )
# dev.off()



# correlation <- by_polygon_drive_day %>%
#     select(
#         co2d_ppm_ex,
#         pm25_ugm3_ex,
#         nox_ppb_ex,
#         no_ppb_ex,
#         no2_ppb_ex,
#         bc_ngm3_ex,
#         ch4d_ppm_ex,
#         vulcan_total_ppm_ex = total_ppm_ex,
#         income,
#         race_pct_white,
#         race_pct_minority
#     ) %>%
#     as.matrix() %>%
#     rcorr()
# r <- correlation$r
# p <- correlation$P
# png('img/correlation_matrix/by_polygon_drive_day_subset.png', width = 1000, height = 1000, res = 144)
# corrplot(
#     r,
#     type = 'upper',
#     method = 'color',
#     addCoef.col = 'black',
#     tl.col = 'black',
#     p.mat = p,
#     sig.level = 0.05,
# )
# dev.off()


library(arrow)
library(corrplot)
library(Hmisc)
library(tidyverse)


by_polygon <- read_feather("gsv-data/by_polygon/by_polygon.feather") %>%
    mutate(nox_ppb_ex = no_ppb_ex + no2_ppb_ex) %>%
    select(
        co2d_ppm_ex,
        pm25_ugm3_ex,
        nox_ppb_ex,
        no_ppb_ex,
        no2_ppb_ex,
        bc_ngm3_ex,
        ch4d_ppm_ex,
        everything()
    )

correlation <- select(by_polygon, -polygon_id) %>%
    as.matrix() %>%
    rcorr()
r <- correlation$r
p <- correlation$P
png("img/correlation_matrix/by_polygon.png", width = 3000, height = 3000, res = 144)
corrplot(
    r,
    type = "upper",
    method = "color",
    addCoef.col = "black",
    tl.col = "black",
    p.mat = p,
    sig.level = 0.05,
)
dev.off()



correlation <- by_polygon %>%
    select(
        co2d_ppm_ex,
        pm25_ugm3_ex,
        nox_ppb_ex,
        no_ppb_ex,
        no2_ppb_ex,
        bc_ngm3_ex,
        ch4d_ppm_ex,
        vulcan_total_ppm_ex = total_ppm_ex,
        income,
        race_white_pct,
        race_minority_pct
    ) %>%
    as.matrix() %>%
    rcorr()
r <- correlation$r
p <- correlation$P
png("img/correlation_matrix/by_polygon_subset.png", width = 1000, height = 1000, res = 144)
try(corrplot(
    r,
    type = "upper",
    method = "color",
    addCoef.col = "black",
    tl.col = "black",
    p.mat = p,
    sig.level = 0.05,
))
dev.off()



by_polygon_drive_day <- read_feather("gsv-data/by_polygon_drive_day/by_polygon_drive_day.feather") %>%
    mutate(nox_ppb_ex = no_ppb_ex + no2_ppb_ex) %>%
    select(
        co2d_ppm_ex,
        pm25_ugm3_ex,
        nox_ppb_ex,
        no_ppb_ex,
        no2_ppb_ex,
        bc_ngm3_ex,
        ch4d_ppm_ex,
        everything()
    )
correlation <- select(by_polygon_drive_day, -polygon_id, -system_id, -time) %>%
    as.matrix() %>%
    rcorr()
r <- correlation$r
p <- correlation$P
png("img/correlation_matrix/by_polygon_drive_day.png", width = 3000, height = 3000, res = 144)
corrplot(
    r,
    type = "upper",
    method = "color",
    addCoef.col = "black",
    tl.col = "black",
    p.mat = p,
    sig.level = 0.05,
)
dev.off()



correlation <- by_polygon_drive_day %>%
    select(
        co2d_ppm_ex,
        pm25_ugm3_ex,
        nox_ppb_ex,
        no_ppb_ex,
        no2_ppb_ex,
        bc_ngm3_ex,
        ch4d_ppm_ex,
        vulcan_total_ppm_ex = total_ppm_ex,
        income,
        race_white_pct,
        race_minority_pct
    ) %>%
    as.matrix() %>%
    rcorr()
r <- correlation$r
p <- correlation$P
png("img/correlation_matrix/by_polygon_drive_day_subset.png", width = 1000, height = 1000, res = 144)
try(corrplot(
    r,
    type = "upper",
    method = "color",
    addCoef.col = "black",
    tl.col = "black",
    p.mat = p,
    sig.level = 0.05,
))
dev.off()
