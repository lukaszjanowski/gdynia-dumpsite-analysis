source("00_config.R", encoding = "UTF-8")
dems<-rast(file.path(DIRS['rasters'],"Gdynia_DEM_2011_2019_2023_harmonized.tif"))
dump_f<-file.path(DIRS['vectors'],"dumpsite.gpkg"); dumpsite<-if(file.exists(dump_f)) vect(dump_f) else NULL
cell_area<-prod(res(dems))
summary_one<-function(dod,pair,lod,zone="full"){
 v<-values(dod,mat=FALSE); v<-v[is.finite(v)]; vt<-v[abs(v)>=lod]
 data.table(pair=pair,zone=zone,lod_m=lod,total_valid_area_m2=length(v)*cell_area,
 detectable_area_m2=length(vt)*cell_area,detectable_pct=100*length(vt)/length(v),
 raising_area_m2=sum(vt>0)*cell_area,lowering_area_m2=sum(vt<0)*cell_area,
 raising_volume_m3=sum(vt[vt>0])*cell_area,lowering_volume_m3=-sum(vt[vt<0])*cell_area,
 net_volume_m3=sum(vt)*cell_area,total_change_volume_m3=sum(abs(vt))*cell_area,
 mean_raising_m=if(any(vt>0)) mean(vt[vt>0]) else NA_real_,mean_lowering_m=if(any(vt<0)) mean(vt[vt<0]) else NA_real_,
 q02=quantile(v,.02),q50=median(v),q98=quantile(v,.98),min=min(v),max=max(v))
}
process_pair<-function(i,j,pair){
 dod<-dems[[j]]-dems[[i]]; names(dod)<-paste0("DoD_",pair)
 writeRaster(dod,file.path(DIRS['rasters'],paste0(names(dod),"_raw.tif")),overwrite=TRUE,gdal=c("COMPRESS=DEFLATE","TILED=YES"))
 out<-rbindlist(lapply(LOD_VALUES,function(l){
   a<-summary_one(dod,pair,l,"full_common_coverage")
   if(!is.null(dumpsite)){
    inside<-mask(dod,dumpsite); outside<-mask(dod,dumpsite,inverse=TRUE)
    rbind(a,summary_one(inside,pair,l,"inside_dumpsite"),summary_one(outside,pair,l,"outside_dumpsite"))
   } else a
 }))
 # primary classified raster: -1 lowering, 0 stable, 1 raising
 cls<-ifel(abs(dod)<LOD_PRIMARY,0,ifel(dod>0,1,-1)); names(cls)<-paste0("class_",pair,"_LoD",LOD_PRIMARY)
 writeRaster(cls,file.path(DIRS['rasters'],paste0(names(cls),".tif")),overwrite=TRUE,datatype="INT1S",gdal=c("COMPRESS=DEFLATE","TILED=YES"))
 out
}
res<-rbind(process_pair(1,2,"2011_2019"),process_pair(2,3,"2019_2023"),process_pair(1,3,"2011_2023"))
fwrite(res,file.path(DIRS['tables'],"dod_volume_summary_all_lod.csv"))
# official comparison for the two main intervals

# Optional comparison with unverified administrative records skipped.
message(
  "Optional disposal-record comparison skipped: ",
  "no verified administrative input was supplied."
)
