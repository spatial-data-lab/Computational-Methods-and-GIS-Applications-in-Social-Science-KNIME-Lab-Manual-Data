# Library##################################################
library(np)
library(raster)
library(rgdal)
library(spatialEco)
library(rstudioapi)
library(sf)
# Read Data
df<-knime.in
# Set save path
pathspace<-knime.flow.in[["knime.workspace"]]
stkdetif <- paste(pathspace,"/CMGISV3/Chp03/data/resulttif.tif",sep="")

# Define input/output parameters ################
varx <- "X"
vary <- "Y"
time <- "time"
test_start <-  "2011-11-01" # week format 2011-10-27
# Grid Size
gridsize_t <- 7 # 7
gridsize_xy <- 250 # 100 meters
# Bandwidth
bwx <- 359.6448 # 359.6448
bwy <- 702.1862 # 702.1862
bwt <- 230 # 230
probs <- 0.01

### Extract Data s################

df <- df[c(varx, vary, time)]
names(df) <- c("X", "Y", "OFFENSE_DATE")
print(paste("variabels ", varx, vary, time, "renamed as X,Y,OFFENSE_DATE"))

x1 <- df[, "X"]
y1 <- df[, "Y"]
# set time field format
df$OFFENSE_DATE <- as.Date(df$OFFENSE_DATE)
df$T <- as.integer(df$OFFENSE_DATE - min(df$OFFENSE_DATE) + 1)
df_sub <- df[df$OFFENSE_DATE < test_start, ]

x <- df_sub[, "X"]
y <- df_sub[, "Y"]
t <- df_sub[, "T"]

# filter test time range
test_end <- as.Date(test_start) + gridsize_t - 1
df_sub1 <- df[
df$OFFENSE_DATE >= test_start & df$OFFENSE_DATE <= test_end,
]
print(summary(df_sub1))

# Get test data during specific time span
x_test <- df_sub1[, "X"]
y_test <- df_sub1[, "Y"]


x_seq <- seq(floor(min(x1) - 2 * gridsize_xy), ceiling(max(x1) + 2 * gridsize_xy), by = gridsize_xy)
y_seq <- seq(floor(min(y1) - 2 * gridsize_xy), ceiling(max(y1) + 2 * gridsize_xy), by = gridsize_xy)
t_seq <- seq(ceiling(max(t)) + 1, ceiling(max(t)) + gridsize_t, length = gridsize_t)

print(paste(
	"Build a grid with", length(x_seq),
	"rows, ", length(y_seq), "cols, and ",
	length(x_seq) * length(y_seq), "cells in total"
	))

data_eval <- expand.grid(x = x_seq, y = y_seq, t = t_seq)

# Set bandwidth
bw <- npudensbw(~ x + y + t,
	bws = c(bwx, bwy, bwt),
	bwmethod = "cv.ml",
	bandwidth.compute = FALSE,
	ckertype = "epanechnikov"
	)
summary(bw)
fhat <- npudens(bws = bw)
result <- predict(fhat, newdata = data_eval)
df_result <- data.frame(data_eval, d = result)

# define normalization
norm_minmax <- function(x) {(x - min(x)) / (max(x) - min(x))}

# apply normalization
df_result$d <- norm_minmax(df_result$d) * 100
# aggregate mean value for all cells- long table
df_agg <- aggregate(
	list(d = df_result$d),
	by = list(x = df_result$x, y = df_result$y),
	mean
)

crs_obj <- CRS("+init=epsg:26915")
# save it as raster
ras_obj <- rasterFromXYZ(df_agg)
ras_obj <- rasterFromXYZ(df_agg)
proj4string(ras_obj) <- crs_obj
rf <- writeRaster(ras_obj, filename=stkdetif, format="GTiff", overwrite=TRUE)

### Significance test
##################################################

# significance test
sig_test <- function(ras_kde, ras_ref, probs, runs, sig_value) {
	qt <- quantile(ras_kde, probs = probs)
	hot_cells <- Which(ras_kde >= qt, cell = TRUE)
	df_hotcells <- data.frame(cellid = hot_cells, count = 0)
	for (i in 1:runs) {
 		hot_cells_simu <- sampleRandom(
   		ras_ref,
   		size = length(hot_cells), na.rm = TRUE, cells = TRUE
		 )
		simucellids <- hot_cells_simu[, 1]
		df_hotcells$count[df_hotcells$cellid %in% simucellids] <- df_hotcells$count[df_hotcells$cellid %in% simucellids] + 1
 		print(paste(i, " out of ", runs, " simulations"))}
	sigcellids <- df_hotcells$cellid[df_hotcells$count >= sig_value]
	ras_sighot <- ras_kde
	ras_sighot[] <- NA
	ras_sighot[sigcellids] <- ras_kde[sigcellids]
	return(list("cellids" = sigcellids, "sig_raster" = ras_sighot))
}

eva_metrics <- function(x_test, y_test, ras_sighot, probs) {
	df_xy <- data.frame(x = x_test, y = y_test)
	list_hitrate <- c()
	list_arearate <- c()
	list_pai <- c()
	list_nni <- c()
	list_nni_pvalue <- c()

	for (i in seq(0.99, probs, -0.005)) {
 		qt <- quantile(ras_sighot, probs = i)
 		hot_cells <- Which(ras_sighot >= qt, cell = TRUE)
 		ras_sighot_i <- ras_sighot
 		ras_sighot_i[] <- NA
 		ras_sighot_i[hot_cells] <- ras_sighot[hot_cells]
 		ras_ct <- rasterize(df_xy, ras_sighot_i, fun = "count")
 		list_ct <- as.vector(ras_ct[hot_cells])
 		hit_ct <- sum(list_ct, na.rm = TRUE)
 		hit_rate <- hit_ct / length(x_test)
 		list_hitrate <- c(list_hitrate, hit_rate)
 		area_sigcells <- length(hot_cells) * xres(ras_sighot) * yres(ras_sighot)
 		area_allcells <- nrow(ras_sighot) * ncol(ras_sighot) * xres(ras_sighot) * yres(ras_sighot)
 		ar_rate <- area_sigcells / area_allcells
 		list_arearate <- c(list_arearate, ar_rate)
 		list_pai <- c(list_pai, hit_rate / ar_rate)
 		sp <- rasterToPoints(ras_sighot_i, spatial = TRUE)
 		# new code: SpatialPointsDataFrame to POINT
 		sp <- st_as_sf(sp)
 		list_nni <- c(list_nni, nni(sp, "hull")$NNI)
 		list_nni_pvalue <- c(list_nni_pvalue, nni(sp, "hull")$p)
	}
	result <- data.frame(
		arearate = list_arearate,
		hitrate = list_hitrate, pai = list_pai, nni = list_nni,
		nni_p = list_nni_pvalue
	)
	print("done!")
	return(result)
}


sig_result <- sig_test(
	ras_obj, ras_obj,
	probs = probs, runs = 100, sig_value = 90
)
#print(sig_result$sig_raster)

eva_result <- eva_metrics(
	x_test, y_test, sig_result$sig_raster,
	probs = probs
)
