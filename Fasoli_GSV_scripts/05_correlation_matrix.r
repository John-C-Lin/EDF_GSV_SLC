library(arrow); library(corrplot)
library(Hmisc); library(tidyverse)
library(ggplot2)


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
col.jcl <- colorRampPalette(rev(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061")))
corrplot(
    r,
    type = "upper",
    method = "color",
    addCoef.col = "black",
    tl.col = "black",
    p.mat = p,
    sig.level = 0.05,
    col = col.jcl(200)
)
dev.off()


correlation <- by_polygon %>%
    select(
        pm25_ugm3_ex,
        nox_ppb_ex,
        bc_ngm3_ex,
        population_density,
        income,
        race_white_pct,
        race_black_pct,
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
    col = col.jcl(200),
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
    col = col.jcl(200),
)
dev.off()


correlation <- by_polygon_drive_day %>%
    select(
        pm25_ugm3_ex,
        nox_ppb_ex,
        bc_ngm3_ex,
        population_density,
        income,
        race_white_pct,
        race_black_pct,
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
    col = col.jcl(200),
))
dev.off()

# multiple regression of "by_polygon_drive_day" to understand the factors related to pollution
dat.all <- by_polygon_drive_day %>%
    select(
        pm25_ugm3_ex,
        nox_ppb_ex,
        bc_ngm3_ex,
        population_density,
        income,
        race_white_pct,
        race_black_pct,
        race_minority_pct
    )  %>% as.data.frame()
lm.pm25 <- lm(pm25_ugm3_ex ~ population_density + income + race_white_pct + race_black_pct + race_minority_pct,data=dat.all)
slm.pm25 <- step(lm.pm25)   # stepwise selection of regression model by AIC (see MASS 3rd Ed, pg. 186-188)
summary(slm.pm25)
lm.nox <- lm(nox_ppb_ex ~ population_density + income + race_white_pct + race_black_pct + race_minority_pct,data=dat.all)
slm.nox <- step(lm.nox)   # stepwise selection of regression model by AIC (see MASS 3rd Ed, pg. 186-188)
summary(slm.nox)
lm.bc <- lm(bc_ngm3_ex ~ population_density + income + race_white_pct + race_black_pct + race_minority_pct,data=dat.all)
slm.bc <- step(lm.bc)   # stepwise selection of regression model by AIC (see MASS 3rd Ed, pg. 186-188)
summary(slm.bc)


# multiple regression of "by_polygon" to understand the factors related to pollution
dat.all <- by_polygon %>%
    select(
        pm25_ugm3_ex,
        nox_ppb_ex,
        bc_ngm3_ex,
        population_density,
        income,
        race_white_pct,
        race_black_pct,
        race_minority_pct
    )  %>% as.data.frame()
lm.pm25 <- lm(pm25_ugm3_ex ~ population_density + income + race_white_pct + race_black_pct + race_minority_pct,data=dat.all)
slm.pm25 <- step(lm.pm25)   # stepwise selection of regression model by AIC (see MASS 3rd Ed, pg. 186-188)
summary(slm.pm25)
lm.nox <- lm(nox_ppb_ex ~ population_density + income + race_white_pct + race_black_pct + race_minority_pct,data=dat.all)
slm.nox <- step(lm.nox)   # stepwise selection of regression model by AIC (see MASS 3rd Ed, pg. 186-188)
summary(slm.nox)
lm.bc <- lm(bc_ngm3_ex ~ population_density + income + race_white_pct + race_black_pct + race_minority_pct,data=dat.all)
slm.bc <- step(lm.bc)   # stepwise selection of regression model by AIC (see MASS 3rd Ed, pg. 186-188)
summary(slm.bc)

# understand the socioeconomic data in the different polygons
dat <- by_polygon_drive_day
unique(dat$polygon_id)
sel <- dat$polygon_id == "Ivy_League"  # each polygon has constant values over drive/days
dat[sel,"population_density"]
# print out polygons, ordered by diff vars (large to small)
data <- tapply(dat$population_density,dat$polygon_id,unique)
print(rev(sort(data)))
data <- data.frame(polygon_id = names(data), population_density = data)
ggplot(data, aes(x=polygon_id, y=population_density, fill=polygon_id)) + geom_bar(stat = "identity",show.legend=FALSE) + coord_flip()
ggsave("popdens_by_polygon.png",path="img/correlation_matrix")

data <- tapply(dat$income,dat$polygon_id,unique)
print(rev(sort(data)))
data <- data.frame(polygon_id = names(data), income = data)
ggplot(data, aes(x=polygon_id, y=income, fill=polygon_id)) + geom_bar(stat = "identity",show.legend=FALSE) + coord_flip()
ggsave("income_by_polygon.png",path="img/correlation_matrix")

data <- tapply(dat$race_white_pct,dat$polygon_id,unique)
print(rev(sort(data)))
data <- data.frame(polygon_id = names(data), race_white_pct = data)
ggplot(data, aes(x=polygon_id, y=race_white_pct, fill=polygon_id)) + geom_bar(stat = "identity",show.legend=FALSE) + coord_flip()
ggsave("race_white_pct.png",path="img/correlation_matrix")

print(rev(sort(tapply(dat$race_white_pct,dat$polygon_id,unique))))
print(rev(sort(tapply(dat$race_minority_pct,dat$polygon_id,unique))))
print(rev(sort(tapply(dat$race_black_pct,dat$polygon_id,unique))))
data <- data.frame(polygon=names(pm25),pm25_ugm3_ex = pm25)
ggplot(data, aes(x=polygon, y=pm25_ugm3_ex)) + geom_bar(stat = "identity") + coord_flip()

# understand the pollution excess data in the different polygons
pm25 <- rev(sort(tapply(dat$pm25_ugm3_ex,dat$polygon_id,median,na.rm=TRUE)))
nox <- rev(sort(tapply(dat$nox_ppb_ex,dat$polygon_id,median,na.rm=TRUE)))
bc <- rev(sort(tapply(dat$bc_ngm3_ex,dat$polygon_id,median,na.rm=TRUE)))



ggplot(by_polygon_drive_day, aes(x=polygon_id, y=pm25_ugm3_ex, fill=polygon_id)) + 
    geom_boxplot() + theme(legend.position="none") + xlab("") + coord_flip()
ggsave("pm25_by_polygon.png",path="img/correlation_matrix")

ggplot(by_polygon_drive_day, aes(x=polygon_id, y=nox_ppb_ex, fill=polygon_id)) + 
    geom_boxplot() + theme(legend.position="none") + xlab("") + coord_flip()
ggsave("nox_by_polygon.png",path="img/correlation_matrix")

ggplot(by_polygon_drive_day, aes(x=polygon_id, y=bc_ngm3_ex, fill=polygon_id)) + 
    geom_boxplot() + theme(legend.position="none") + xlab("") + coord_flip()
ggsave("bc_by_polygon.png",path="img/correlation_matrix")