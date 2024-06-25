# refactor 1


# load packages -----------------------------------------------------------
library(raster)
library(tidyverse) #dplyr masks raster::extract + raster::select
library(vegan) # for calculate_sv
library(geometry) # for calculate_chv

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

subplot_dimensions <- read.csv('data/fishnets/cons_fishnet.csv') 

# generate subplot IDs
subplot_ids <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep="_")))
subplot_dimensions <- subplot_dimensions %>% mutate(subplot_id = subplot_ids)

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

raster_data <- stack('data_out/combined_rasters/cons_reflectance_combined_image.tif')
names(raster_data) <- c("blue", "green", "red", "red_edge", "nir")

for (i in 1:nrow(subplots)) {
  subplot <- subplots[i, ]
  subplot_id <- subplot$subplot_id
  
  # Convert the current subplot to a Spatial object
  subplot_sp <- as(subplot, "Spatial")
  
  # Crop and mask the raster using the current subplot
  cropped_raster <- crop(raster_data, subplot_sp)
  masked_raster <- mask(cropped_raster, subplot_sp)
  
  # Extract pixel values
  pixel_values <- as.data.frame(getValues(masked_raster))
  
  # Add subplot ID to pixel values data frame
  pixel_values$subplot_id <- subplot_id
  
  # Add to the list
  pixel_values_list[[i]] <- pixel_values
}

# Combine all pixel values into a single data frame
all_pixel_values <- bind_rows(pixel_values_list) %>%
  na.omit() #omits 15012324 NAs?


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


