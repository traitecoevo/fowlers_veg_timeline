# refactor 1

# load packages -----------------------------------------------------------
library(raster)
library(tidyverse) # dplyr masks raster::extract + raster::select
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

# extract pixel values for all rasters by subplot  ------------------------------------------------------

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
    rownames_to_column('subplot_id') %>%
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

combined_metrics <- combined_cv %>%
  left_join(combined_sv, by = c("identifier", "subplot_id")) %>%
  left_join(combined_chv, by = c("identifier", "subplot_id"))

# field observations for every plot ---------------------------------------

# read survey data
survey_data <- read.csv('data/ausplots_march_24.csv')

# get unique site names
unique_sites <- unique(survey_data$site_unique) 
unique_sites <- unique_sites[unique_sites != ""]

# list to store results for all lists
all_site_results <- list()

# list to store community matrices- this is a temporary step, to check that community matrices are correct :)
community_matrices <- list()

# Loop through each unique site
for (site in unique_sites) {
  # Filter data for the current site
  site_survey_data <- survey_data %>%
    filter(site_unique == site)
  
  # Extract only direction of the transect (no numbers)
  site_survey_data$transect_direction <- gsub('[[:digit:]]+', '', site_survey_data$transect)
  
  # Extract only number of the transect (no direction)
  site_survey_data$transect_number <- as.numeric(gsub(".*?([0-9]+).*", "\\1", site_survey_data$transect))
  
  # Create variable for fixed transect direction (to order them all transects in the same direction)
  site_survey_data$transect_direction2 <- NA 
  
  # Create variable for fixed point number (inverse in some cases as if they had been collected in the same direction)
  site_survey_data$point_number2 <- NA 
  
  # Create XY empty variables for plot XY coordinates
  site_survey_data$X_plot <- NA 
  site_survey_data$Y_plot <- NA
  
  # For loop to homogenize transects and numbers. It converts all W-E to E-W and all N-S to S-N 
  for (i in 1:nrow(site_survey_data)){
    if (site_survey_data[i, "transect_direction"] == "W-E") {
      site_survey_data[i, "point_number2"] <- 100 - site_survey_data[i, "point_number"] # If transect E-W, transect fixed is W-E and inverse numbers
      site_survey_data[i, "transect_direction2"] <- "E-W"
    }
    if (site_survey_data[i, "transect_direction"] == "E-W") {
      site_survey_data[i, "point_number2"] <- site_survey_data[i, "point_number"] # If transect W-E, all stays the same
      site_survey_data[i, "transect_direction2"] <- "E-W"
    }
    if (site_survey_data[i, "transect_direction"] == "S-N") {
      site_survey_data[i, "point_number2"] <- site_survey_data[i, "point_number"] # If transect N-S, all stays the same
      site_survey_data[i, "transect_direction2"] <- "S-N"
    }
    if (site_survey_data[i, "transect_direction"] == "N-S") {
      site_survey_data[i, "point_number2"] <- 100 - site_survey_data[i, "point_number"] # If transect S-N, transect fixed is N-S and inverse numbers
      site_survey_data[i, "transect_direction2"] <- "S-N"
    }
  }
  
  # For loop to assign plotXY coordinates to each point intercept
  for (i in 1:nrow(site_survey_data)){
    if (site_survey_data[i, "transect_direction2"] == "E-W") {
      if (site_survey_data[i, "transect_number"] == 1){
        site_survey_data[i, "Y_plot"] <- 10
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 2){
        site_survey_data[i, "Y_plot"] <- 30
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 3){
        site_survey_data[i, "Y_plot"] <- 50
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 4){
        site_survey_data[i, "Y_plot"] <- 70
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 5){
        site_survey_data[i, "Y_plot"] <- 90
        site_survey_data[i, "X_plot"] <- site_survey_data[i, "point_number2"]
      }
    }
    if (site_survey_data[i, "transect_direction2"] == "S-N") {
      if (site_survey_data[i, "transect_number"] == 1){
        site_survey_data[i, "X_plot"] <- 10
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 2){
        site_survey_data[i, "X_plot"] <- 30
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 3){
        site_survey_data[i, "X_plot"] <- 50
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 4){
        site_survey_data[i, "X_plot"] <- 70
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
      if (site_survey_data[i, "transect_number"] == 5){
        site_survey_data[i, "X_plot"] <- 90
        site_survey_data[i, "Y_plot"] <- site_survey_data[i, "point_number2"]
      }
    }
  }
  
  # subplot rows and columns - +1 ensures 0 point values fall into correct subplot, 
  # pmin ensures 100 point values falls in correct subplot given +1
  site_survey_data$subplot_row <- pmin(ceiling((site_survey_data$Y_plot + 1) / 20), 5)
  site_survey_data$subplot_col <- pmin(ceiling((site_survey_data$X_plot + 1) / 20), 5)
  
  # single ID for subplot row and column
  site_survey_data$subplot_id <- paste(site_survey_data$subplot_row, site_survey_data$subplot_col, sep = "_")
  
  subplot_diversity <- site_survey_data %>%
    drop_na(standardised_name) %>%
    group_by(subplot_id) %>%
    summarise(species_richness = n_distinct(standardised_name))
  
  community_matrix <- site_survey_data %>%
    drop_na(standardised_name) %>%
    count(subplot_id, standardised_name) %>%
    spread(standardised_name, n, fill = 0)
  
  # store the community matrix in the list - this is a temp step to check!!!
  community_matrices[[site]] <- community_matrix
  
   #Remove unwanted column if it exists -- WHY DOES THIS COLUMN EXIST~!!!!>???>
  if ("V1" %in% colnames(community_matrix)) {
    community_matrix <- community_matrix %>% 
      dplyr::select(-V1)
  }
  
  # calculate diversity indices
  shannon_diversity <- diversity(community_matrix[, -1], index = "shannon")
  simpson_diversity <- diversity(community_matrix[, -1], index = "simpson")
  
  subplot_diversity <- subplot_diversity %>%
    mutate(shannon_diversity = shannon_diversity,
      simpson_diversity = simpson_diversity,
      site = site)
  
  # store  result for  current site
  all_site_results[[site]] <- subplot_diversity
}

# combine into one df
final_results <- bind_rows(all_site_results, .id = "site")

# add identifier based on site name so that this can be joined with spectral results 
# in refactor two you could remove this step by just calling the sites by their ausplot name...
# would require changing the file names of OG rasters
final_results <- final_results %>%
  mutate(identifier = recode(site,
                       "NSABHC009" = "emu",
                       "NSABHC010" = "emgrazed",
                       "NSABHC011" = "sandstone",
                       "NSABHC012" = "cons"))

head(final_results)


# erghhhhhhhhhhhhhhhhhhhh ... need to still clean up the survey data and think about stuff like... 
# dead shrub stuff, things IDd to genus level etc etc....
#community_matrices[["NSABHC012"]] |> View()

# combined field and spectral metrics -------------------------------------

combined_diversity_values <- left_join(final_results, combined_metrics, by = c('identifier', 'subplot_id'))

# pivot longer 
long_data <- combined_diversity_values %>%
  pivot_longer(cols = c(species_richness, shannon_diversity, simpson_diversity),
               names_to = "taxonomic_metric",
               values_to = "taxonomic_value") %>%
  pivot_longer(cols = c(CV, SV, CHV),
               names_to = "spectral_metric",
               values_to = "spectral_value")


# plot, showing spectral ~ taxonomic relos
ggplot(long_data, aes(x = spectral_value, y = taxonomic_value)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(taxonomic_metric ~ spectral_metric, scales = "free") +
  labs(
    title = "comparison of taxonomic + spectral diversity metrics ",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# now plot seperating plot location (emgrazed removed as don't have this raster yet lol oops)
long_data %>%
  filter(identifier != 'emgrazed') %>% 
  ggplot(aes(x = spectral_value, y = taxonomic_value)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE) +
  facet_grid(taxonomic_metric ~ spectral_metric + identifier, scales = "free") +
  labs(
    title = "comparison of taxonomic + spectral diversity metrics by site",
    x = "spectral diversity value",
    y = "taxonomic diversity value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

