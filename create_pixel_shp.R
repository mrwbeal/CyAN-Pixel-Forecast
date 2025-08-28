library(sf)
library(feather)
library(dplyr)
library(ggplot2)

conus_lakes <- read_sf('data/OLCI_resolvable_lakes_2022_09_08/OLCI_resolvable_lakes_2022_09_08.shp')



f = list.files("/work/HAB4CAST/max_beal/CyANPixelForecast/data/cyan_pixel_conus_lakes",full.names = T)
df = read_feather(f[1])

getpixelshp =function(file){
  df = read_feather(file)
  cells = df %>% distinct(x,y,cell,COMID)
  return(cells)
}


pixelshp = list()
for (file in seq_along(f)) {
  print(file)
  pixelshp[[file]] = getpixelshp(f[file])
    
}

cells = pixelshp %>% bind_rows()


sf_object <- st_as_sf(cells, coords = c("x", "y"), crs = st_crs(conus_lakes))

st_write(sf_object,"/work/HAB4CAST/max_beal/CyANPixelForecast/data/cells_shapefile.shp")

# Test plots
# pixel = sf_object %>% filter(COMID==10038590)
# 
# lake = conus_lakes %>% filter(COMID==10038590)
# 
# ggplot() + geom_sf(data=lake) + geom_sf(data=pixel)
# 
# 
# 
# 
# temptest = mean_atemp%>% filter(COMID==10038590)
# pixel = conus_cells %>% filter(COMID==10038590)
# 
# temp_plot = merge(conus_cells,temptest)
# 
# temp_plot = st_buffer(temp_plot,100) #Buffer cell points by 1 m, exactextract only accepts polygons
# 
# exrast = crop(atemp_spat,temp_plot)
# 
# rastplt = as.data.frame(exrast,xy=T)
# 
# ggplot() + geom_tile(data=rastplt,aes(x=x,y=y,fill=PRISM_tmean_stable_4kmD2_20170101_bil)) + geom_sf(data=lake) + geom_sf(data=temp_plot, aes(fill=mean_atemp)) + scale_fill_viridis_c()
# 



