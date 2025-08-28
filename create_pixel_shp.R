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



