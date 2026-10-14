source("00_config.R", encoding = "UTF-8")
paths<-readRDS(file.path(OUT,"paths.rds")); dems<-rast(file.path(DIRS['rasters'],"Gdynia_DEM_2011_2019_2023_harmonized.tif"))
# Morphometry from harmonized 2023 DEM
z<-dems[[3]]
slope<-terrain(z,"slope",unit="degrees",neighbors=8)
tri<-terrain(z,"TRI",neighbors=8)
rough<-terrain(z,"roughness",neighbors=8)
# multiscale bathymetric position: 11 m and 51 m neighbourhoods
bpi<-function(x,n){x-focal(x,w=matrix(1,n,n),fun="mean",na.rm=TRUE)}
bpi11<-bpi(z,11); bpi51<-bpi(z,51)
m<-c(z,slope,tri,rough,bpi11,bpi51); names(m)<-c("DEM_2023","slope_deg","TRI","roughness","BPI_11m","BPI_51m")
writeRaster(m,file.path(DIRS['rasters'],"Gdynia_2023_morphometry.tif"),overwrite=TRUE,gdal=c("COMPRESS=DEFLATE","TILED=YES"))
# Align acoustic layers for interpretation only; do not treat RGB SSS as quantitative intensity.
template<-z
read_align<-function(f,method="bilinear") resample(project(rast(f),crs(template),method=method),template,method=method)
bbs<-read_align(paths$bbs2023); names(bbs)<-paste0("BBS_2023_",seq_len(nlyr(bbs)))
sss11<-read_align(paths$sss2011); names(sss11)<-paste0("SSS_2011_RGB",seq_len(nlyr(sss11)))
sss19<-read_align(paths$sss2019); names(sss19)<-paste0("SSS_2019_RGB",seq_len(nlyr(sss19)))
writeRaster(bbs,file.path(DIRS['rasters'],"BBS_2023_aligned_interpretive.tif"),overwrite=TRUE)
writeRaster(sss11,file.path(DIRS['rasters'],"SSS_2011_aligned_interpretive.tif"),overwrite=TRUE)
writeRaster(sss19,file.path(DIRS['rasters'],"SSS_2019_aligned_interpretive.tif"),overwrite=TRUE)
# Zonal morphometry inside/outside if boundary exists
if(file.exists(file.path(DIRS['vectors'],"dumpsite.gpkg"))){
 d<-vect(file.path(DIRS['vectors'],"dumpsite.gpkg")); inside<-global(mask(m,d),c("mean","sd","min","max"),na.rm=TRUE); fwrite(as.data.table(inside,keep.rownames="metric"),file.path(DIRS['tables'],"morphometry_inside_dumpsite.csv"))
}
