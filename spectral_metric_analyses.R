# refactor 1

# load packages -----------------------------------------------------------
library(raster)
library(tidyverse) #dplyr masks raster::extract + raster::select
library(vegan) # for calculate_sv
library(geometry) # for calculate_chv
library(sf)

# combining spectral band images into one tif file  -------------------------------------------------------

# file directory
mosaics_dir <- "data/2024_mosaics"

# folder list | recursive = won't pick folders within folders
folders <- list.dirs(mosaics_dir, full.names = FALSE, recursive = FALSE)
folders <- folders[folders != "point_clouds"]

# define order of bands
desired_band_order <- c("blue", "green", "red", "red_edge", "nir")

## NOTE: spectral band image tif file names must be named after their band (e.g., blue, nir, etc), 
#  otherwise change 'desired_band_order' to match file names
#  should be combined in wavelength order, esp for biodivmapR processes (i.e. as above)


# read in fishnet shapefile from fishnets dir, only want geometry col
subplots <- read_sf('data/fishnets/cons_fishnet.shp') %>%
  select('geometry')

# apply subplot ids
subplots$subplot_id <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep="_")))

# loop thru each folder 
for (folder in folders) {
  # create path
  folder_path <- file.path(mosaics_dir, folder)
  
  # list of tif files | \\. represents . (dots need to be escaped w \, \ need to be escaped with  \). $means at end of file name/string
  tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)
  
  # load as raster
  rasters <- lapply(tif_files, raster)
  
  # extract band names from file names
  band_names <- tools::file_path_sans_ext(basename(tif_files))
  
  # stack rasters and assign band names
  combined_image <- stack(rasters)
  names(combined_image) <- band_names
  
  # reorder the bands based on the desired band order
  combined_image <- combined_image[[match(desired_band_order, band_names)]]
  
  # create output file
  output_filename <- file.path("data_out/combined_rasters", paste0(folder, "_combined_image.tif"))
  writeRaster(combined_image, filename = output_filename, format = "GTiff", options="INTERLEAVE=BAND", overwrite = TRUE)
  
  plot(combined_image)
}

# spectral metric functions by subplot ------------------------------------

# read in fishnet shapefile from fishnets dir, only want geometry col
subplots <- read_sf('data/fishnets/cons_fishnet.shp') %>%
  select('geometry')

# apply subplot ids
subplots$subplot_id <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep="_")))

## check that subplots ids have been applied correctly
# calc centroids for easy vis
centroids <- st_centroid(subplots) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  mutate(subplot_id = subplots$subplot_id)

# map fishnet, make sure ids applied correctly
ggplot() +
  geom_sf(data = subplots) +
  geom_text(data = centroids, aes(x = X, y = Y, label = subplot_id), color = 'red') +
  labs(title = 'polygons with subplot id', x = 'lon', y = 'lat') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# load raster
raster_data <- stack('data_out/combined_rasters/cons_reflectance_combined_image.tif')

# correct names
names(raster_data) <- c("blue", "green", "red", "red_edge", "nir")

for (i in 1:nrow(subplots)) {
  subplot <- subplots[i, ]
  subplot_id <- subplot$subplot_id
  
  # convert subplot to spatial object
  subplot_sp <- as(subplot, "Spatial")
  
  # crop and mask the raster using the current subplot
  cropped_raster <- crop(raster_data, subplot_sp)
  masked_raster <- mask(cropped_raster, subplot_sp)
  
  # extract pixel values
  pixel_values <- as.data.frame(getValues(masked_raster))
  
  # add subplot id to pixel values df
  pixel_values$subplot_id <- subplot_id
  
  # add to list
  pixel_values_list[[i]] <- pixel_values
}

# combine all pixel values into single df
all_pixel_values <- bind_rows(pixel_values_list) %>%
  na.omit() #omits 15012324 NAs?

## CV, CHV, SV metric functions appropriated from https://github.com/ALCrofts/CABO_SVH_Forest_Sites/tree/v1.0

# cv function
calculate_cv <- function(pixel_values_df, # Dataframe with spectra reflectance values 
                    subplots, # What you want to calculate spectral diversity for, here it's plots.
                    wavelengths # Cols where spectral reflectance values are
){
cv <- pixel_values_df %>%
  select(c({{subplots}}, {{wavelengths}})) %>%
  group_by({{subplots}}) %>%
  summarise_all(~sd(.)/abs(mean(.))) %>%
  rowwise({{subplots}}) %>%
  summarise(CV = sum(c_across(cols = everything()), na.rm = T) / (ncol(.) - sum(is.na(c_across(everything())))))

return(cv)
}

cv_cons <- calculate_cv(all_pixel_values, subplot_id, c('blue':'nir'))

# sv function
calculate_sv <- function(pixel_values_df, # Dataframe with spectra reflectance values 
                         subplots, # What you want to calculate spectral diversity for, here it's plots.
                         wavelengths # Cols where spectral reflectance values are
){
  
  # a) Determine the number of spectral observations per plot
  spectral_points <- pixel_values_df %>% 
    group_by({{subplots}}) %>%
    summarise(points = n())
  
sv <- pixel_values_df %>%
  select(c({{wavelengths}}, {{subplots}})) %>%
  group_by({{subplots}}) %>%
  summarise_all(~sum((.x - mean(.x))^2)) %>%
  rowwise({{subplots}}) %>%
  summarise(SS = sum(c_across(cols = everything()))) %>%
  left_join(spectral_points) %>%
  summarise(SV = SS / (points - 1)) # divide by the number of spectral observations, here it can vary by plot

return(sv)
}

sv_cons <- calculate_sv(all_pixel_values, subplot_id, c('blue':'nir'))

# chv function
calculate_chv <- function(pixel_values_df,
                          subplots,
                          wavelengths
){
  
  subplot_of_interest <- deparse(substitute(subplots))
  
  PCA <- pixel_values_df %>%
    select(c({{wavelengths}})) %>%
    rda(scale = F) 
  
  CHV <- data.frame(PCA$CA$u) %>%
    select(c(1:3)) %>%
    cbind({{pixel_values_df}}) %>%
    select(c('PC1', 'PC2', 'PC3', {{subplots}})) %>%
    group_split({{subplots}}) %>%
    set_names(map(., ~unique(.[[subplot_of_interest]]))) %>%
    map(~convhulln(.x[-4], option = 'FA')) %>%
    map_dbl('vol') %>%
    as.data.frame() %>%
    rownames_to_column('plot_field_id') %>%
    rename(CHV = '.')
  
  rm(subplot_of_interest) 
  
  return(CHV)
}

chv_cons <- calculate_chv(all_pixel_values, subplot_id, c('blue':'nir'))


# now for ALL images and fishnets ------------------------------------------------------

raster_files <- list.files('data_out/combined_rasters', pattern = '_reflectance_combined_image.tif$', full.names = T)
subplot_files <- list.files('data/fishnets', pattern = '_fishnet.shp$', full.names = T)

all_pixel_values_list <- list()

for (raster_file in raster_files) {
  
  # identify the string that represents the site name
  identifier <- str_extract(basename(raster_file), "^[^_]+")
  
  #choose the corresponding subplot file
  subplot_file <- subplot_files[grep(paste0('^', identifier), basename(subplot_files))]

  # read in subplot file and select geometries
  subplots <- read_sf(subplot_file) %>% 
    select('geometry')
  
  # apply subplot ids
  subplots$subplot_id <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep="_")))
  
  # read in raster file
  raster_data <- stack(raster_file)
  
  # apply names - should be saved in wavelength order as per sect 1 of this script
  names(raster_data) <- c('blue', 'green', 'red', 'red_edge', 'nir')
  
  # create empty list
  pixel_values_list <- list()
  
  for (i in 1:nrow(subplots)){
    
    # select the i-th subplot and its id
    subplot <- subplots[i, ]
    subplot_id <- subplot$subplot_id
    
    # convert to spatial object
    subplot_sp <- as(subplot, "Spatial")
    
    # crop and mask raster using current subplot 
    cropped_raster <- crop(raster_data, subplot_sp)
    masked_raster <- mask(cropped_raster, subplot_sp)
    
    # extract pixel values
    pixel_values  <- as.data.frame(getValues(masked_raster))
    
    # add subplot id to pixel values df
    pixel_values$subplot_id <- subplot_id
    
    #add to list
    pixel_values_list[[i]] <- pixel_values
    
  }
  # combined all pixel values into one df for current raster
  all_pixel_values <- bind_rows(pixel_values_list) %>%
    na.omit()
  
  # add to overall list with all raster data pixel values 
  all_pixel_values_list[[identifier]] <- all_pixel_values

  }

final_pixel_values <- bind_rows(all_pixel_values_list, .id = 'identifier')



# apply spectral metrics to pixel value df --------------------------------

## CV, CHV, SV metric functions appropriated from https://github.com/ALCrofts/CABO_SVH_Forest_Sites/tree/v1.0
# cv function
calculate_cv <- function(pixel_values_df, subplots, wavelengths) {
  cv <- pixel_values_df %>%
    select(c({{subplots}}, {{wavelengths}})) %>%
    group_by({{subplots}}) %>%
    summarise_all(~sd(.)/abs(mean(.))) %>%
    rowwise({{subplots}}) %>%
    summarise(CV = sum(c_across(cols = everything()), na.rm = T) / (ncol(.) - sum(is.na(c_across(everything())))))
  
  return(cv)
}

# sv function
calculate_sv <- function(pixel_values_df, subplots, wavelengths) {
  spectral_points <- pixel_values_df %>%
    group_by({{subplots}}) %>%
    summarise(points = n())
  
  sv <- pixel_values_df %>%
    select(c({{wavelengths}}, {{subplots}})) %>%
    group_by({{subplots}}) %>%
    summarise_all(~sum((.x - mean(.x))^2)) %>%
    rowwise({{subplots}}) %>%
    summarise(SS = sum(c_across(cols = everything()))) %>%
    left_join(spectral_points) %>%
    summarise(SV = SS / (points - 1))
  
  return(sv)
}

# chv function
calculate_chv <- function(pixel_values_df, subplots, wavelengths) {
  subplot_of_interest <- deparse(substitute(subplots))
  
  PCA <- pixel_values_df %>%
    select(c({{wavelengths}})) %>%
    rda(scale = F)
  
  CHV <- data.frame(PCA$CA$u) %>%
    select(c(1:3)) %>%
    cbind({{pixel_values_df}}) %>%
    select(c('PC1', 'PC2', 'PC3', {{subplots}})) %>%
    group_split({{subplots}}) %>%
    set_names(map(., ~unique(.[[subplot_of_interest]]))) %>%
    map(~convhulln(.x[-4], option = 'FA')) %>%
    map_dbl('vol') %>%
    as.data.frame() %>%
    rownames_to_column('plot_field_id') %>%
    rename(CHV = '.')
  
  rm(subplot_of_interest) 
  
  return(CHV)
}

# create empty list to store results
results <- list()

# loop through each site (represented as 'identifier' from file name)
for (identifier in unique(final_pixel_values$identifier)) {
  
  # !! --> gets the values of current identifier ('unquotes' it)
  site_pixel_values <- final_pixel_values %>% filter(identifier == !!identifier)
  
  cv <- calculate_cv(site_pixel_values, subplot_id, c('blue', 'green', 'red', 'red_edge', 'nir'))
  sv <- calculate_sv(site_pixel_values, subplot_id, c('blue', 'green', 'red', 'red_edge', 'nir'))
  chv <- calculate_chv(site_pixel_values, subplot_id, c('blue', 'green', 'red', 'red_edge', 'nir'))
  
  results[[identifier]] <- list(CV = cv, SV = sv, CHV = chv)
}

# combine into dfs
combined_cv <- bind_rows(lapply(results, function(x) x$CV), .id = 'identifier')
combined_sv <- bind_rows(lapply(results, function(x) x$SV), .id = 'identifier')
combined_chv <- bind_rows(lapply(results, function(x) x$CHV), .id = 'identifier')

