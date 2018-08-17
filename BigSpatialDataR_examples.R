## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ------------------------------------------------------------------------
list.of.packages <- c("plyr", "parallel", "GSIF", "ranger", "raster", 
                      "rgdal", "rgrass7", "snowfall", "lidR", "knitr", "tmap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies = TRUE)

## ------------------------------------------------------------------------
library(rgdal)
library(raster)
GDALinfo("./data/Indonesia_ESA_lcc_300m_2000.tif")

## ------------------------------------------------------------------------
make_LC_tiles <- function(i, tile.tbl, 
                          out.path="./tiled",
                          lc1="./data/Indonesia_ESA_lcc_300m_2000.tif",
                          lc2="./data/Indonesia_ESA_lcc_300m_2015.tif", 
                          leg.lcc){
  out.tif = paste0(out.path, "/T_", tile.tbl[i,"ID"], ".tif")
  if(!file.exists(out.tif)){
    m <- readGDAL(lc1, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]),
                 region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]),
                 output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]),
                 silent = TRUE)
    m@data[,2] <- readGDAL(lc2, offset=unlist(tile.tbl[i,c("offset.y","offset.x")]),
                  region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]),
                  output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]), 
                  silent = TRUE)$band1
    names(m) <- c("LC2000","LC2015")
    m <- as(m, "SpatialPixelsDataFrame")
    ## Focus only on pixels that show land cover change
    sel <- !m$LC2000==m$LC2015
    if(sum(sel)>0){
      m <- m[sel,]
      m$v <- paste(m$LC2000, m$LC2015, sep="_")
      m$i <- plyr::join(data.frame(NAME=m$v), leg.lcc, type="left")$Value
      writeGDAL(m["i"], out.tif, type="Int16", 
                options="COMPRESS=DEFLATE", mvFlag=-32768)
    }
  }
}

## ------------------------------------------------------------------------
leg <- read.csv("./data/ESA_landcover_legend.csv")
str(leg)
comb.leg <- expand.grid(leg$Value, leg$Value)
comb.leg$lcc <- paste(comb.leg$Var1, comb.leg$Var2, sep="_")

## ------------------------------------------------------------------------
leg.lcc <- data.frame(Value=1:nrow(comb.leg), NAME=comb.leg$lcc)
head(leg.lcc)

## ------------------------------------------------------------------------
library(raster)
library(GSIF)
## check whether the maps match perfectly to the same grid:
x <- raster::stack(paste0("./data/Indonesia_ESA_lcc_300m_", c(2000, 2015), ".tif"))
## OK!
obj <- GDALinfo("./data/Indonesia_ESA_lcc_300m_2000.tif")
## tile to 50km blocks:
tile.lst <- getSpatialTiles(obj, block.x=.5, return.SpatialPolygons=TRUE)
tile.tbl <- getSpatialTiles(obj, block.x=.5, return.SpatialPolygons=FALSE)
tile.tbl$ID <- as.character(1:nrow(tile.tbl))
head(tile.tbl)

## ----plot-tiles, echo=TRUE, fig.width=5, fig.cap="Tiling system based on the 50 km by 50 km tiles."----
te <- as.vector(extent(x))
library(tmap)
data("World")
tm_shape(World, xlim=te[c(1,2)], ylim=te[c(3,4)], projection="longlat") +
  tm_polygons() +
  tm_shape(as(tile.lst, "SpatialLines")) + tm_lines()

## ----plot-tile-lcc, echo=TRUE, fig.width=5, fig.cap="Single tile loaded into memory and plotted."----
## plot tile number 124:
i = 124
m <- readGDAL("./data/Indonesia_ESA_lcc_300m_2000.tif", 
              offset=unlist(tile.tbl[i,c("offset.y","offset.x")]),
              region.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]),
              output.dim=unlist(tile.tbl[i,c("region.dim.y","region.dim.x")]))
plot(raster(m), legend=FALSE, col=rgb(leg$R/255, leg$G/255, leg$B/255))

## ------------------------------------------------------------------------
library(snowfall)
sfInit(parallel=TRUE, cpus=parallel::detectCores())
sfExport("make_LC_tiles", "tile.tbl", "leg.lcc")
sfLibrary(rgdal)
sfLibrary(plyr)
out.lst <- sfClusterApplyLB(1:nrow(tile.tbl), 
                function(x){ make_LC_tiles(x, tile.tbl, leg.lcc=leg.lcc) })
sfStop()
## takes few seconds depending on the number of cores

## ----htop-8-cores, echo=FALSE, fig.cap="Fully parallelized computing using 8 cores. Displayed using htop software.", out.width="100%"----
knitr::include_graphics("./tex/htop_8_cores.png")

## ------------------------------------------------------------------------
t.lst <- list.files("./tiled", pattern=".tif", full.names=TRUE)
str(t.lst)

## ------------------------------------------------------------------------
out.tmp <- "./data/t_list.txt"
vrt.tmp <- "./data/indonesia.vrt"
cat(t.lst, sep="\n", file=out.tmp)
system(paste0('gdalbuildvrt -input_file_list ', out.tmp, ' ', vrt.tmp))
system(paste0('gdalwarp ', vrt.tmp, 
             ' \"./data/Indonesia_ESA_lcc_300m_change.tif\" ', 
             '-ot \"Int16\" -dstnodata \"-32767\" -co \"BIGTIFF=YES\" ',  
             '-multi -wm 2000 -co \"COMPRESS=DEFLATE\" -overwrite ',
             '-r \"near\" -wo \"NUM_THREADS=ALL_CPUS\"'))

## ------------------------------------------------------------------------
raster("./data/Indonesia_ESA_lcc_300m_change.tif")

## ----qgis-indonesia, echo=FALSE, fig.cap="Land cover class changes (2000 to 2015) for Kalimantan.", out.width="100%"----
knitr::include_graphics("./tex/Indonesia_ESA_lcc_300m_change_preview.jpg")

## ---- eval=FALSE---------------------------------------------------------
## laz2grid <- function(file, res=30, out.dir, prj){
##  out.file = paste0(out.dir, gsub(".laz", ".tif", basename(file)))
##  if(!file.exists(out.file)){
##    require(lidR); require(raster); require(sp)
##    lidar = lidR::readLAS(file)
##    dtm = lidR::grid_terrain(lidar, res=res, method = "knnidw")
##    dtm.r = as.raster(dtm)
##    raster::projection(dtm.r) <- paste(lidar@crs)
##    dtm.s = raster::projectRaster(dtm.r, crs=prj, res=30)
##    raster::writeRaster(dtm.s, filename=out.file, options="COMPRESS=DEFLATE", overwrite=TRUE)
##  }
## }

## ------------------------------------------------------------------------
pnt.ices <- read.csv("./data/icesat_1533635398523.csv")
pnt.ices <- pnt.ices[,1:3]
names(pnt.ices) <- c("x","y","mDLSM")
pnt.ices$type <- "ICESAT"
coordinates(pnt.ices) <- ~ x+y
proj4string(pnt.ices) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
pnt.lid <- sampleRandom(raster("./data/Boulder_LiDAR_30m.tif"), size=4e4, sp=TRUE)
names(pnt.lid) <- "mDLSM"
pnt.lid$type <- "LiDAR"
pnt.ices <- spTransform(pnt.ices, pnt.lid@proj4string)
pnt.all <- rbind(pnt.ices, pnt.lid)
## filter out values outside the range:
pnt.all <- pnt.all[pnt.all$mDLSM<5000&pnt.all$mDLSM>1000,]
## 48,649 points

## ------------------------------------------------------------------------
library(parallel)
library(raster)
cov.lst = c("Boulder_AW3D30s_30m_v1802.tif", 
            "Boulder_NED_30m.tif", "Boulder_Landsat_red_30m_2014.tif",
            "Boulder_Landsat_NIR_30m_2014.tif",
            "Boulder_HH_30m.tif", "Boulder_HV_30m.tif")
xB <- raster::stack(paste0("./data/", cov.lst))
## all covariates can be stacked!
ov.cov <- parallel::mclapply(cov.lst, function(i){ 
                  raster::extract(raster(paste0("./data/", i)), pnt.all) }, 
                  mc.cores = detectCores())
ov.cov <- data.frame(ov.cov)
names(ov.cov) = cov.lst

## ------------------------------------------------------------------------
rm.DLSM <- cbind(as.data.frame(pnt.all), ov.cov)
rm.DLSM$NDVI <- (rm.DLSM[,8] - rm.DLSM[,7]) / (rm.DLSM[,8] + rm.DLSM[,7]) * 100
str(rm.DLSM)

## ------------------------------------------------------------------------
library(ranger)
fm.DLSM <- mDLSM ~ Boulder_AW3D30s_30m_v1802.tif + Boulder_NED_30m.tif +
  Boulder_HH_30m.tif + Boulder_HV_30m.tif + NDVI
sel.na <- complete.cases(rm.DLSM[,all.vars(fm.DLSM)])
summary(sel.na)
m.DLSM <- ranger(fm.DLSM, rm.DLSM[sel.na,], num.trees=85, importance="impurity")
m.DLSM

## ------------------------------------------------------------------------
xl1.P <- as.list(ranger::importance(m.DLSM))
print(t(data.frame(xl1.P[order(unlist(xl1.P), decreasing=TRUE)])))

## ------------------------------------------------------------------------
objB <- GDALinfo("./data/Boulder_AW3D30s_30m_v1802.tif")
library(GSIF)
## tile to 10km blocks:
tileB.lst <- getSpatialTiles(objB, block.x=1e4, return.SpatialPolygons=TRUE)
tileB.tbl <- getSpatialTiles(objB, block.x=1e4, return.SpatialPolygons=FALSE)
tileB.tbl$ID <- as.character(1:nrow(tileB.tbl))

## ------------------------------------------------------------------------
predict_mDLSM <- function(m.DLSM, i, tileB.tbl, cov.lst, in.path="./data/", out.path="./tiledB/"){
  out.tif = paste0(out.path, "T_", tileB.tbl[i,"ID"], ".tif")
  if(!file.exists(out.tif)){
    covs.files <- paste0(in.path, cov.lst)
    newdata <- readGDAL(covs.files[1], offset=unlist(tileB.tbl[i,c("offset.y","offset.x")]),
                 region.dim=unlist(tileB.tbl[i,c("region.dim.y","region.dim.x")]),
                 output.dim=unlist(tileB.tbl[i,c("region.dim.y","region.dim.x")]),
                 silent = TRUE)
    for(j in 2:length(cov.lst)){
      newdata@data[,j] <- readGDAL(covs.files[j],
                 offset=unlist(tileB.tbl[i,c("offset.y","offset.x")]),
                 region.dim=unlist(tileB.tbl[i,c("region.dim.y","region.dim.x")]),
                 output.dim=unlist(tileB.tbl[i,c("region.dim.y","region.dim.x")]),
                 silent = TRUE)$band1
    }
    names(newdata) <- basename(covs.files)
    newdata$NDVI <- (newdata@data[,4] - newdata@data[,3]) / 
        (newdata@data[,4] + newdata@data[,3]) * 100
    newdata$NDVI <- ifelse(is.na(newdata$NDVI), 0, newdata$NDVI)
    sel.pr <- complete.cases(newdata@data)
    out <- predict(m.DLSM, newdata@data[sel.pr,])
    g <- as(newdata[1], "SpatialPixelsDataFrame")
    g[sel.pr,"m"] = out$predictions
    writeGDAL(g["m"], out.tif, type="Int16", 
                options="COMPRESS=DEFLATE", mvFlag=-32768)
  }
}

## ------------------------------------------------------------------------
library(snowfall)
sfInit(parallel=TRUE, cpus=parallel::detectCores())
sfExport("predict_mDLSM", "m.DLSM", "tileB.tbl", "cov.lst")
sfLibrary(rgdal)
sfLibrary(ranger)
outB.lst <- sfClusterApplyLB(1:nrow(tileB.tbl), 
                function(x){ predict_mDLSM(m.DLSM, x, tileB.tbl, cov.lst) })
sfStop()
## takes few minutes

## ------------------------------------------------------------------------
tB.lst <- list.files("./tiledB", pattern=".tif", full.names=TRUE)
outB.tmp <- "./data/b_list.txt"
vrtB.tmp <- "./data/boulder.vrt"
cat(tB.lst, sep="\n", file=outB.tmp)
system(paste0('gdalbuildvrt -input_file_list ', outB.tmp, ' ', vrtB.tmp))
system(paste0('gdalwarp ', vrtB.tmp, 
             ' \"./data/Boulder_mDLSM_30m.tif\" ', 
             '-ot \"Int16\" -dstnodata \"-32767\" -co \"BIGTIFF=YES\" ',  
             '-multi -wm 2000 -co \"COMPRESS=DEFLATE\" -overwrite ',
             '-r \"near\" -wo \"NUM_THREADS=ALL_CPUS\"'))

## ----comparison-shading, echo=FALSE, fig.cap="Comparison of the original AW3D30 vs the predicted Land Surface Model (mDLSM). Fitting a model to predict land surface model seems to solve the problem of artifacts / striping effects, but then it can introduce local artifacts in areas of higher vegetation and under-represented by training points.", out.width="100%"----
knitr::include_graphics("./tex/comprison_shading.png")

## ---- eval=FALSE---------------------------------------------------------
## tmp = tempfile(fileext = ".sdat")
## system(paste0('gdal_translate ./data/Boulder_mDLSM_30m.tif -of \"SAGA\" ', tmp))
## #system(paste0('saga_cmd ta_hydrology 15 -DEM=\"',
## #               gsub(".sdat", ".sgrd", tmp), '\" -TWI=\"',
## #               gsub(".sdat", "_twi.sgrd", tmp), '\"') )
## system(paste0('saga_cmd ta_morphometry 0 -ELEVATION=\"',
##                         gsub(".sdat", ".sgrd", tmp), '\" -SLOPE=\"',
##                         gsub(".sdat", "_slope.sgrd", tmp),
##                         '\" -C_PROF=\"',
##                         gsub(".sdat", "_cprof.sgrd", tmp), '\"') )
## system(paste0('gdal_translate ', gsub(".sdat", "_slope.sdat", tmp),
##               ' ./data/Boulder_mDLSM_slope_30m.tif -scale ',
##               ' -ot "Byte" -co \"COMPRESS=DEFLATE\"'))

## ---- eval=FALSE---------------------------------------------------------
## library(rgrass7)
## rname <- "./data/Boulder_mDLSM_30m.tif"
## # Set GRASS environment and database location
## loc <- initGRASS("/usr/lib/grass74", home="/data/tmp/",
##                  gisDbase="GRASS_TEMP", override=TRUE)
## execGRASS("r.in.gdal", flags="o", parameters=list(input=rname, output="mDLSM"))
## execGRASS("g.region", parameters=list(raster="mDLSM"))
## execGRASS("r.geomorphon", parameters=list(elevation="mDLSM", forms="mDLSMg"))
## #plot(readRAST("mDLSMg"))
## execGRASS("r.out.gdal", parameters=list(input="mDLSMg",
##           output="./data/Boulder_mDLSMg_30m.tif",
##           type="Byte", createopt="COMPRESS=DEFLATE"))
## ## clean-up
## unlink("./GRASS_TEMP", recursive = TRUE)
## unset.GIS_LOCK()
## unlink_.gislock()
## remove_GISRC()

