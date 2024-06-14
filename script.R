
# load packages -----------------------------------------------------------

library(sf)
library(ggplot2)
library(galah)
library(tidyverse)
library(austraits)
library(ggthemes)
library(here)
library(devtools)
install_github("ternaustralia/ausplotsR", build_vignettes = TRUE, dependencies = TRUE)
library(ausplotsR)
library(raster)
library(RStoolbox)
library(terra)

# ausplot data ------------------------------------------------------------

#chose fowlers plots based on THIS map of plots https://www.tern.org.au/news-gap-filling/

plots_oi <- c('NSABHC0009', 'NSABHC0010', 'NSABHC0011', 'NSABHC0012', 'NSABHC0028', 'NSABHC0029')

veg <- get_ausplots(plots_oi, veg.vouchers = T, veg.PI = T)

vegPI <- veg$veg.PI

allveg <- veg$veg.vouch

#species richness for plots @ each survey 
veg_count_by_survey <- allveg %>%
  group_by(site_unique) %>%
  summarise(count = n()) 

#all species found across ausplot surveys, and surveys in which they were found
veg_in_plots <- allveg %>% 
  group_by(herbarium_determination) %>% 
  summarize(site_unique = paste(unique(site_unique), collapse = ", ")) %>%
  arrange(herbarium_determination)

#load in 2022 data
vegPI22 <- read.csv("data/FG_2022_visit/Fowlers Gap_PI_GrowthForm_data-ALL.csv") %>%
  subset(substring(visit_start_date, 8, 9) == '22')

#load in 2024 data
vegPI24 <- read.csv("data/ausplot_survey_march24.csv")



#this is for mapping the ausplot surveys on a cartesian plane - it corrects for transects being surveyed in diff directions
# Extract only direction of the transect (no numbers)
vegPI$transect_direction <- gsub('[[:digit:]]+', '', vegPI$transect)

# Extract only number of the transect (no direction)
vegPI$transect_number <- as.numeric(gsub(".*?([0-9]+).*", "\\1", vegPI$transect))

# Create variable for fixed transect direction (to order them all transects in the same direction)
vegPI$transect_direction2 <- NA 

# Create variable for fixed point number (inverse in some cases as if they had been collected in the same direction)
vegPI$point_number2 <- NA 

# Create XY empty variables for plot XY coordinates
vegPI$X_plot <- NA 
vegPI$Y_plot <- NA

# For loop to homogenize transects and numbers. It converts all E-W to W-E and all S-N to N-S
for (i in 1:nrow(vegPI)){
  if (vegPI[i, "transect_direction"] == "E-W") {
    vegPI[i, "point_number2"] <- 101 - vegPI[i, "point_number"] # If transect E-W, transect fixed is W-E and inverse numbers
    vegPI[i, "transect_direction2"] <- "W-E"
  }
  if (vegPI[i, "transect_direction"] == "W-E") {
    vegPI[i, "point_number2"] <- vegPI[i, "point_number"] # If transect W-E, all stays the same
    vegPI[i, "transect_direction2"] <- "W-E"
  }
  if (vegPI[i, "transect_direction"] == "N-S") {
    vegPI[i, "point_number2"] <- vegPI[i, "point_number"] # If transect N-S, all stays the same
    vegPI[i, "transect_direction2"] <- "N-S"
  }
  if (vegPI[i, "transect_direction"] == "S-N") {
    vegPI[i, "point_number2"] <- 101 - vegPI[i, "point_number"] # If transect S-N, transect fixed is N-S and inverse numbers
    vegPI[i, "transect_direction2"] <- "N-S"
  }
}

# For loop to assign plotXY coordinates to each point intercept
for (i in 1:nrow(vegPI)){
  if (vegPI[i, "transect_direction2"] == "W-E") {
    if (vegPI[i, "transect_number"] == 1){
      vegPI[i, "Y_plot"] <- 10
      vegPI[i, "X_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 2){
      vegPI[i, "Y_plot"] <- 30
      vegPI[i, "X_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 3){
      vegPI[i, "Y_plot"] <- 50
      vegPI[i, "X_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 4){
      vegPI[i, "Y_plot"] <- 70
      vegPI[i, "X_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 5){
      vegPI[i, "Y_plot"] <- 90
      vegPI[i, "X_plot"] <- vegPI[i, "point_number2"]
    }
  }
  if (vegPI[i, "transect_direction2"] == "N-S") {
    if (vegPI[i, "transect_number"] == 1){
      vegPI[i, "X_plot"] <- 10
      vegPI[i, "Y_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 2){
      vegPI[i, "X_plot"] <- 30
      vegPI[i, "Y_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 3){
      vegPI[i, "X_plot"] <- 50
      vegPI[i, "Y_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 4){
      vegPI[i, "X_plot"] <- 70
      vegPI[i, "Y_plot"] <- vegPI[i, "point_number2"]
    }
    if (vegPI[i, "transect_number"] == 5){
      vegPI[i, "X_plot"] <- 90
      vegPI[i, "Y_plot"] <- vegPI[i, "point_number2"]
    }
  }
}

#mapping 2012 emu pen survey from ausplots by family
vegPI %>% 
  filter(site_location_visit_id == 53604) %>% #emu pen visit
  drop_na(herbarium_determination) %>%
  ggplot(aes(x = X_plot, y = Y_plot, shape = family, color = family)) +
  geom_point() +
  scale_shape_manual(values = c(1:19)) +
  theme_classic() +
  labs(y = 'West (m)', x = 'South (m)')


# get ALA data and apply threat status ------------------------------------

galah_config(email = 'adelegemmell@hotmail.com')

FG_area <- st_read('data/unsw-fowlers.kml')

fowlers_veg <- galah_call() |>
  galah_identify("plantae") |>
  galah_geolocate(FG_area) |>
  galah_select(genus, family, group = 'basic') |>
  atlas_occurrences()

fowlers_veg$year <- as.numeric(substr(fowlers_veg$eventDate, 1, 4))

#view occurrence of ALA obs by year, filled by family
fowlers_veg %>%
  filter(year >= 1940) %>%
  ggplot() +
  geom_bar(aes(x = year, fill = family)) +
  theme_classic() +
  theme(legend.position = 'none')

# downloaded threatened status of NSW/SA plants from https://www.environment.gov.au/sprat-public/action/report 
threatflora <- read.csv('data/EPBC Threatened Flora plus SA.csv')

#create count df with counts of # records in ALA data
veg_count <- fowlers_veg %>%
  group_by(scientificName) %>%
  summarise(count = n()) 

#create count df with counts of # records in ausplot data
ausplot_counts <- vegPI %>%
  group_by(herbarium_determination) %>%
  summarise(count = n())

# find the matching conservation status for a given species name and status column (doesn't account for synonyms etc)
#i want to run threatflora, and veg_count through APCalign so that synonyms etc are accounted for
findConservationStatus <- function(species_name, status_column) {
  matching_index <- which(species_name %in% threatflora$scientificName |
                            species_name %in% threatflora$EPBC_Name |
                            species_name %in% threatflora$NSW_Name |
                            species_name %in% threatflora$IUCN_Name |
                            species_name %in% threatflora$SA_Name)
  
  if (length(matching_index) > 0) {
    # If there is a match, return the corresponding conservation status from the specified column
    return(threatflora[[status_column]][matching_index[1]])
  } else {
    # If no match is found, return NA or any default value
    return(NA)
  }
}

# create  new column for EPBC, NSW, IUCN, SA threat status
veg_count <- veg_count %>%
  mutate(EPBC_Status = mapply(findConservationStatus, scientificName, "EPBC_Status"),
         NSW_Status = mapply(findConservationStatus, scientificName, "NSW_Status"),
         IUCN_Status = mapply(findConservationStatus, scientificName, "IUCN_Status"),
         SA_Status = mapply(findConservationStatus, scientificName, "SA_Status"),
         SA = mapply(findConservationStatus, scientificName, "SA"),
         NSW = mapply(findConservationStatus, scientificName, "NSW"))


# austraits - extract annual traits ---------------------------------------

library('austraits')
austraits <- load_austraits(version = "5.0.0", path = "intro/downloads")

#extract life history traits
annual_perennial_traits <- austraits %>% extract_trait('life_history')

#down with tibbles
annual_perennial_traits <- left_join(annual_perennial_traits[["traits"]], annual_perennial_traits[["taxa"]], by = "taxon_name")

#unique taxa with life history traits
length(unique(annual_perennial_traits$taxon_name))
#30474

#different life history values
unique(annual_perennial_traits$value)
# [1] "perennial"                                       "biennial perennial"       # [3] "annual"                                          "annual perennial"         # [5] "biennial"                                        "annual biennial"          # [7] "annual short_lived_perennial"                    "ephemeral"                # [9] "short_lived_perennial"                           "annual ephemeral"        #[11] "annual biennial perennial"                       "biennial short_lived_perennial"                 
#[13] "annual ephemeral perennial"                      "annual biennial short_lived_perennial"          
#[15] "perennial short_lived_perennial"                 "biennial perennial short_lived_perennial"       
#[17] "annual perennial short_lived_perennial"          "ephemeral short_lived_perennial"                
#[19] "ephemeral perennial"                             "annual biennial perennial short_lived_perennial"

annual_trait <- unique(annual_perennial_traits$value[grepl("annual", annual_perennial_traits$value)])

fowlers_life_form <- annual_perennial_traits %>%
  filter(taxon_name %in% fowlers_veg$scientificName |
           taxon_name %in% ausplot_counts$herbarium_determination)

#this finds the life form traits listed for each species at fowlers (if available) and selects the life form that is 
#specified most often as the 'life form'. Do we want to do this like this, or take any record of 'annual'?
most_frequent_lifeform <- fowlers_life_form %>%
  group_by(taxon_name) %>%
  summarise(most_common_lifeform = names(sort(table(value), decreasing = TRUE))[1])

fowlers_veg <- fowlers_veg %>%
  left_join(most_frequent_lifeform, by = c("scientificName" = "taxon_name"))

fowlers_veg <- fowlers_veg %>%
  mutate(simp_lf = case_when(
    grepl('annual', most_common_lifeform, ignore.case = TRUE) ~ 'annual',
    grepl('perennial', most_common_lifeform, ignore.case = TRUE) ~ 'perennial',
    TRUE ~ most_common_lifeform
  ))

vegPI <- vegPI %>%
  left_join(most_frequent_lifeform, by = c("herbarium_determination" = "taxon_name"))


vegPI <- vegPI %>%
  mutate(simp_lf = case_when(
    grepl('annual', most_common_lifeform, ignore.case = TRUE) ~ 'annual',
    grepl('perennial', most_common_lifeform, ignore.case = TRUE) ~ 'perennial',
    TRUE ~ most_common_lifeform
  ))

#save csvs for ausplots/ALA + lifeform
#write.csv(vegPI, 'data_out/ausplots_species_with_lifeform')

#write.csv(fowlers_veg, 'data_out/ALA_FG_species_with_lifeform')


# plotting the ausplot species spatially ----------------------------------

#the emu pen site falls within:
#POLYGON((141.69736861 -31.08023472,141.69750889 -31.07934056,141.69855806 -31.07946139,141.69843861 -31.08032806,141.69736861 -31.08023472))

#convert plot based on SW coordinate

# Install and load necessary packages
install.packages("sp")
install.packages("rgdal")
library(sp)
library(rgdal)

# Define the geographic coordinates (latitude and longitude)
latitude <- -31.08023472  # Replace with the actual latitude of your site
longitude <- 141.69736861  # Replace with the actual longitude of your site

# Create a spatial points object
coordinates <- SpatialPoints(matrix(c(longitude, latitude), ncol = 2), 
                             proj4string = CRS("+proj=longlat +datum=WGS84"))

# Define the coordinate reference system (CRS) for Western NSW, Australia
crs <- CRS("+proj=utm +zone=54 +south +datum=WGS84")

# Transform the coordinates to eastings and northings
transformed_coordinates <- spTransform(coordinates, crs)

# Extract eastings and northings
eastings <- coordinates(transformed_coordinates)[1]
northings <- coordinates(transformed_coordinates)[2]

#new df for emu pen only data
emu_pen_veg <- vegPI %>%
  filter(site_location_name == 'NSABHC0009')

#adding SW corner eastings and northings to X and Y values - although this is not a legit fix, as the polygon isnt due north
emu_pen_veg$X_plot <- emu_pen_veg$X_plot + eastings
emu_pen_veg$Y_plot <- emu_pen_veg$Y_plot + northings

# Create a SpatialPointsDataFrame
coordinates <- cbind(emu_pen_veg$X_plot, emu_pen_veg$Y_plot)
crs <- CRS("+proj=utm +zone=54 +south +datum=WGS84")  # Replace with the appropriate CRS for your data
spatial_points <- SpatialPointsDataFrame(coordinates, emu_pen_veg, proj4string = crs)

# Transform to WGS84 (latitude and longitude)
wgs84_crs <- CRS("+proj=longlat +datum=WGS84")
spatial_points_wgs84 <- spTransform(spatial_points, wgs84_crs)

# Extract latitude and longitude
emu_pen_veg$Longitude <- coordinates(spatial_points_wgs84)[, 1]
emu_pen_veg$Latitude <- coordinates(spatial_points_wgs84)[, 2]

library(sf)

FG_WG <- read.csv('data/FG_WG_locations.csv')

coordinates <- cbind(FG_WG$Northing, FG_WG$Easting)
crs <- CRS("+proj=utm +zone=54 +south +datum=WGS84")  # Replace with the appropriate CRS for your data
spatial_points <- SpatialPointsDataFrame(coordinates, FG_WG, proj4string = crs)

# Transform to WGS84 (latitude and longitude)
wgs84_crs <- CRS("+proj=longlat +datum=WGS84")
spatial_points_wgs84 <- spTransform(spatial_points, wgs84_crs)

# Extract latitude and longitude
FG_WG$Longitude <- coordinates(spatial_points_wgs84)[, 1]
FG_WG$Latitude <- coordinates(spatial_points_wgs84)[, 2]

FG_WG$Rainfall_Gauge <- "Rainfall Gauge"

ggplot() +
  geom_sf(data = FG_area, fill = NA) +
  geom_point(data = fowlers_veg, alpha = 0.6, color = 'grey50', 
              aes(x = decimalLongitude, y = decimalLatitude, fill = "ALA obs")) +
  geom_point(data = FG_WG, size = 2, color = 'red',
             aes(x = Longitude, y = Latitude, shape = Rainfall_Gauge)) +
  geom_point(data = veg$site.info, size = 6, alpha = 0.6, shape = 15,
             aes(x = longitude, y = latitude, color = site_location_name)) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = 'white'),  # Set the panel background color to white
plot.background = element_rect(fill = "white", color = 'white'),   # Set the overall plot background color to white
legend.background = element_rect(fill = "white", color = 'white')) +
  scale_color_hue(labels = c('emu pen', 'emu grazed', 'sandstone', 'conservation', 'connors (mallee)', 'freislich (red gum)')) +
  scale_fill_manual(labels = c("ALA observations"), values = c("grey50")) +
  guides(fill = guide_legend(title = NULL),
         shape = guide_legend(title = NULL)) +
  labs(color = "ausplot site locations",
       x = "",
       y = "")

ggsave('maps_graphs/fowlers_obs_map.png')

unique(fowlers_veg$dataResourceName)

ggplot() +
  geom_sf(data = FG_area, fill = NA) +
  geom_point(data = fowlers_veg, alpha = 0.6, 
             aes(x = decimalLongitude, y = decimalLatitude, fill = simp_lf)) +
  geom_point(data = veg$site.info, size = 6, alpha = 0.6,
             aes(x = longitude, y = latitude, color = site_location_name)) +
  geom_point(data = FG_WG, 
             aes(x = Longitude, y = Latitude)) +
  scale_fill_manual(values = c("red", "orange", "green")) +
  theme_minimal() 

vegPI %>% filter(site_location_visit_id == 53607) |>
  View()



# comparing rainfall to site visit dates ----------------------------------

ausplot_site_dates <- read.csv('data/FG_2022_visit/Fowlers_gap_summary_site_data.csv')

FG_daily_rainfall <- read.csv('data/fowlers_rainfall.csv')

FG_daily_rainfall$Short.Date <- as.Date(paste(FG_daily_rainfall$Year, 
                                               FG_daily_rainfall$Month, 
                                               FG_daily_rainfall$Day, 
                                               sep = "-"), 
                                         format = "%Y-%m-%d")


#pivot longer
ausplot_site_dates_long <- ausplot_site_dates %>% 
  pivot_longer(cols = c(Established.Date, Revisit.Date.1, Revisit.Date.2, Revisit.Date.3),
               names_to = "Visit_Number",
               values_to = "Visit_Date") %>%
  filter(Visit_Date != "") %>%  # Remove rows with empty values in Visit_Date
  arrange(Plot.Name, Visit_Date) %>%
  mutate(Visit_Number = case_when(
    str_detect(Visit_Number, "Established") ~ 1,
    str_detect(Visit_Number, "Revisit.Date.1") ~ 2,
    str_detect(Visit_Number, "Revisit.Date.2") ~ 3,
    str_detect(Visit_Number, "Revisit.Date.3") ~ 4
  ))

#as.Date the date
ausplot_site_dates_long$Visit_Date <- as.Date(ausplot_site_dates_long$Visit_Date, format = "%d/%m/%Y")


ggplot() +
  geom_point(data = ausplot_site_dates_long, aes(x = Visit_Date, y = 100, color = Plot.Name), size = 3, position = position_jitter(height = 25)) +
  geom_line(data = FG_daily_rainfall[FG_daily_rainfall$Year >= 2012, ],
            aes(x = Short.Date, y = Rainfall)) + 
  scale_x_date(date_labels = "%Y-%m-%d") +
  scale_x_date(date_labels = "%Y-%m-%d", limits = c(as.Date("2012-01-01"), NA), date_breaks = "6 months") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

library(dplyr)
library(purrr)

# Assuming ausplot_site_dates_long is your main data frame
# Assuming FG_daily_rainfall is your rainfall data frame

# Function to calculate rainfall sum for previous n months
calculate_previous_rainfall <- function(visit_date, rainfall_data, n_months) {
  end_date <- visit_date
  start_date <- end_date %m-% months(n_months)
  
  subset_data <- rainfall_data %>%
    filter(Short.Date >= start_date, Short.Date <= end_date)
  
  sum_rainfall <- sum(subset_data$Rainfall, na.rm = TRUE)
  return(sum_rainfall)
}

# Calculate rainfall for previous 1, 3, 6, and 12 months
ausplot_site_dates_long <- ausplot_site_dates_long %>%
  mutate(
    Rainfall_1_month = pmap_dbl(list(Visit_Date, list(FG_daily_rainfall), 1), calculate_previous_rainfall),
    Rainfall_3_months = pmap_dbl(list(Visit_Date, list(FG_daily_rainfall), 3), calculate_previous_rainfall),
    Rainfall_6_months = pmap_dbl(list(Visit_Date, list(FG_daily_rainfall), 6), calculate_previous_rainfall),
    Rainfall_12_months = pmap_dbl(list(Visit_Date, list(FG_daily_rainfall), 12), calculate_previous_rainfall)
  )

# Print the resulting data frame
print(ausplot_site_dates_long)


# Given date
given_date <- as.Date("2024-02-01", format = "%Y-%m-%d")

# Calculate rainfall for previous 1, 3, 6, and 12 months
rainfall_for_given_date <- data.frame(
  Date = given_date,
  Rainfall_1_month = calculate_previous_rainfall(given_date, FG_daily_rainfall, 1),
  Rainfall_3_months = calculate_previous_rainfall(given_date, FG_daily_rainfall, 3),
  Rainfall_6_months = calculate_previous_rainfall(given_date, FG_daily_rainfall, 6),
  Rainfall_12_months = calculate_previous_rainfall(given_date, FG_daily_rainfall, 12)
)


# Function to calculate the number of days with rainfall greater than 0 for previous n months
calculate_previous_rainfall_days <- function(end_date, rainfall_data, n_months) {
  start_date <- end_date %m-% months(n_months)
  
  subset_data <- rainfall_data %>%
    filter(Short.Date >= start_date, Short.Date <= end_date, Rainfall > 0)
  
  num_rainfall_days <- nrow(subset_data)
  return(num_rainfall_days)
}

# Calculate number of rainfall days for each visit
ausplot_site_dates_long <- ausplot_site_dates_long %>%
  mutate(
    RainfallDays_1_month = map_int(Visit_Date, ~ calculate_previous_rainfall_days(.x, FG_daily_rainfall, 1)),
    RainfallDays_3_months = map_int(Visit_Date, ~ calculate_previous_rainfall_days(.x, FG_daily_rainfall, 3)),
    RainfallDays_6_months = map_int(Visit_Date, ~ calculate_previous_rainfall_days(.x, FG_daily_rainfall, 6)),
    RainfallDays_12_months = map_int(Visit_Date, ~ calculate_previous_rainfall_days(.x, FG_daily_rainfall, 12))
  )

ggplot(ausplot_site_dates_long, aes(x = Plot.Name, y = Rainfall_1_month, fill = as.factor(Visit_Date))) +
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Reshape the data for facet_wrap
ausplot_site_dates_long_long <- pivot_longer(
  ausplot_site_dates_long,
  cols = c(Rainfall_1_month, Rainfall_3_months, Rainfall_6_months, Rainfall_12_months),
  names_to = "Rainfall_Type",
  values_to = "Rainfall"
)

# Create the plot
ggplot(ausplot_site_dates_long_long, aes(x = Plot.Name, y = Rainfall, fill = as.factor(Visit_Date))) +
  geom_col(position = "dodge") +
  facet_wrap(~ Rainfall_Type, scales = "free_y", nrow = 2, ncol = 2) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = 'Visit Date')

# Reshape the data for facet_wrap
ausplot_site_dates_long_long <- pivot_longer(
  ausplot_site_dates_long,
  cols = c(Rainfall_1_month, Rainfall_3_months, Rainfall_6_months, Rainfall_12_months,
           RainfallDays_1_month, RainfallDays_3_months, RainfallDays_6_months, RainfallDays_12_months),
  names_to = "Rainfall_Type",
  values_to = "Rainfall"
)


library(stringr)  
# Define custom labels based on Rainfall_Type
custom_labels <- function(x) {
  ifelse(
    grepl("^Rainfall_", x),
    gsub("^(Rainfall_[0-9]+_months)$", "Rainfall (mm) \\1", x),
    gsub("_", " ", x)
  )
}

ausplot_site_dates_long_long$Rainfall_Type <- factor(
  ausplot_site_dates_long_long$Rainfall_Type,
  levels = c("Rainfall_1_month", "Rainfall_3_months", "Rainfall_6_months", "Rainfall_12_months",
             "RainfallDays_1_month", "RainfallDays_3_months", "RainfallDays_6_months", "RainfallDays_12_months"),
  labels = custom_labels,
  ordered = TRUE
)

ausplot_site_dates_long_long %>%
  filter(!(Plot.Name %in% c('NSABHC0028', 'NSABHC0029'))) %>% 
  ggplot(aes(x = Plot.Name, y = Rainfall, fill = as.factor(substr(Visit_Date, 1, 4)))) +
  geom_col(position = "dodge") +
  facet_wrap(~ Rainfall_Type, scales = "free_y", nrow = 2, ncol = 4) +
  scale_fill_manual(values = c("#169c9c", "#b39eb5", "#9ec35d", "#f08080")) +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.background = element_rect(fill = "white"), 
    plot.background = element_rect(fill = "white"),   
    legend.background = element_rect(fill = "white")
  ) +
  labs(fill = 'Visit Year', x = 'Site Name')


# plant census ------------------------------------------------------------

install.packages("remotes")

install.packages("arrow")
library(arrow)

remotes::install_github("traitecoevo/APCalign")
library(APCalign)

resources <- load_taxonomic_resources()

tax_lookup <- create_taxonomic_update_lookup( 
  taxa = unique(c(vegPI$herbarium_determination,
           vegPI22$field_name,
           vegPI24$field_name
  )),
  resources = resources
)

threat_tax_lookup <- create_taxonomic_update_lookup( 
  taxa = threatflora$scientificName,
  resources = resources
)


# downloading thomas' inat observations -----------------------------------


galah_config(email = 'adelegemmell@hotmail.com')

FG_area <- st_read('data/unsw-fowlers.kml')

thomas_obs <- galah_call() |>
  galah_identify('plantae') |>
  galah_geolocate(FG_area) |>
  galah_filter(userId == 'thebeachcomber', year == 2024, month == 3) |>
  galah_select(genus, month, year, family, userId, group = 'basic') |>
  atlas_occurrences() 

write.csv(thomas_obs, 'data_out/thomas_obs.csv')


# combining rasters -------------------------------------------------------


library(raster)

# file directory
mosaics_dir <- "data/2024_mosaics"

# define order of bands
desired_band_order <- c("blue", "green", "red", "red_edge", "nir")

# folder list
folders <- list.dirs(mosaics_dir, full.names = FALSE)
#folders <- folders[folders != ""]
folders <- folders[folders == 'cons_reflectance']


# loop thru each folder 
for (folder in folders) {
  # create path
  folder_path <- file.path(mosaics_dir, folder)
  
  # list of tif files
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
  #output_filename <- file.path("data_out", paste0(folder, "_combined_image.tif"))
  #writeRaster(combined_image, filename = output_filename, format = "GTiff", options="INTERLEAVE=BAND", overwrite = TRUE)
  
  plot(combined_image)



# optimum nir/ndvi values -------------------------------------------------

ndvi_values <- read.csv('data/ndvi_less20_both_table.csv')
nir_values <- read.csv('data/nir_less_05_table.csv')

ggplot(nir_values, aes(x = class, y = nir)) +
  geom_boxplot() +
  facet_wrap(~ location) +
  labs(x = "Cover Type", y = "Reflectance Value") +
  theme_minimal()

ggplot(nir_values %>% filter(class != "unclear"), aes(x = nir, fill = class)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of NIR Reflectance Values by Class",
       x = "NIR Reflectance",
       y = "Density",
       fill = "Class") +
  facet_wrap(~location)

ggplot(ndvi_values %>% filter(class != "unclear"), aes(x = ndviPix, fill = class)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of NDVI Reflectance Values by Class",
       x = "NIR Reflectance",
       y = "Density",
       fill = "Class") +
facet_wrap(~location) 


# biodivmapR take 2 -------------------------------------------------------

# load biodivMapR and useful libraries 
remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR', force = T)
library(biodivMapR)
library(sf)
library(stars)
library(utils)
library(terra)
library(raster)
library(sf)

#clip cons raster to polygon size
cons_stack <- stack('data_out/cons_reflectance_combined_image.tif')
cons_polygon <- st_read('data/cons_shapefile/cons_poly_24.shp') %>% st_zm()

# match crs
cons_polygon <- st_transform(cons_polygon, crs(cons_stack))

# clip raster to polygon 
clipped_cons <- terra::mask(cons_stack, cons_polygon)

# reflectance values
c_reflectance_values <- getValues(clipped_cons)

# convert to 16 bit values
c_reflectance_values_16 <- round(c_reflectance_values * 10000)

#checking no values are > 65535 (eg. are unsigned, eg floating value = pos)
max(c_reflectance_values_16 %>% na.omit())

# cons raster with 16 bit values
clipped_cons_16 <- setValues(clipped_cons, c_reflectance_values_16)

names(clipped_cons_16) <- c("blue", "green", "red", "red_edge", "nir")

plot(clipped_cons_16)

# save
writeRaster(clipped_cons_16, filename = "data_out/cons_reflectance_combined_image_16", format = "GTiff", options="INTERLEAVE=BAND", datatype='INT2U', overwrite=TRUE)

rasttif <- terra::rast('data_out/cons_reflectance_combined_image_16.tif')
names(rasttif) <- c("blue", "green", "red", "red_edge", "nir")
terra::writeRaster(x = rasttif, filename = "ENVI/cons_reflectance_combined_image_16", filetype = 'ENVI', overwrite = T)

tmpdir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/ENVI'
NameRaster <- 'cons_reflectance_combined_image_16'
destfiletiff <- file.path(tmpdir,NameRaster)

Input_Image_File <- destfiletiff


# Set to FALSE if no mask available
Input_Mask_File <- FALSE

Output_Dir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/biodivmapR/RESULTS'

NDVI_Thresh <- 0.05
Blue_Thresh <- 1500
NIR_Thresh <- (0.02*10000)

# continuum removal is a normalisation procedure which reduces multiplicative effects
Continuum_Removal <- TRUE

TypePCA <- 'SPCA'

# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- FALSE

window_size <- 10

nbCPU <- 10
MaxRAM <- 0.5

# number of clusters
nbclusters <- 50

Excluded_WL <- NA

Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = Input_Mask_File,
                                                 Output_Dir = Output_Dir, TypePCA = TypePCA,
                                                 NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                 NIR_Thresh = NIR_Thresh,
                                                 Blue = 450,
                                                 Red = 650,
                                                 NIR = 840)

PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File, 
                          Input_Mask_File = Input_Mask_File,
                          Output_Dir = Output_Dir, 
                          TypePCA = TypePCA, 
                          FilterPCA = FilterPCA,
                          nbCPU = nbCPU, 
                          MaxRAM = MaxRAM, 
                          Continuum_Removal = Continuum_Removal)

# path of the raster resulting from dimensionality reduction
PCA_Files <- PCA_Output$PCA_Files
# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File,
                                Output_Dir = Output_Dir, 
                                PCA_Files = PCA_Output$PCA_Files,
                                TypePCA = PCA_Output$TypePCA, 
                                File_Open = TRUE)

Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File,
                                    Input_Mask_File = Input_Mask_File,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output, 
                                    nbclusters = nbclusters, 
                                    nbCPU = nbCPU, MaxRAM = MaxRAM)



Index_Alpha <- c('Shannon')

window_size <- 1250
map_alpha_div(Input_Image_File = Input_Image_File, 
              Output_Dir = Output_Dir, 
              TypePCA = TypePCA,
              window_size = window_size, 
              nbCPU = nbCPU, 
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha, 
              nbclusters = nbclusters)
??map_alpha_div

con_shan_10 <- rast('biodivmapR/RESULTS/cons_reflectance_combined_image_16/SPCA/ALPHA/Shannon_10')
con_simp_10 <- rast('biodivmapR/RESULTS/cons_reflectance_combined_image_16/SPCA/ALPHA/Simpson_10')

plot(con_shan_10)
plot(con_simp_10)

con_shan_1250 <- rast('biodivmapR/RESULTS/cons_reflectance_combined_image_16/SPCA/ALPHA/Shannon_1250')

plot(con_shan_1250)

# could i add a red filter for bare ground to 
# perform_radiometric_filtering via create_mask_from_threshold function?


# biodivmapR again - unclipped --------------------------------------------

# load biodivMapR and useful libraries 
remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR', force = T)
library(biodivMapR)
library(sf)
library(stars)
library(utils)
library(terra)
library(raster)

tmpdir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/ENVI'
NameRaster <- 'cons_reflectance_combined_image_2'
destfiletiff <- file.path(tmpdir,NameRaster)

Input_Image_File <- destfiletiff

# Set to FALSE if no mask available
Input_Mask_File <- FALSE

Output_Dir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/biodivmapR/RESULTS'

NDVI_Thresh <- 0.05
Blue_Thresh <- 1500/10000
NIR_Thresh <- 0.02

# continuum removal is a normalisation procedure which reduces multiplicative effects
Continuum_Removal <- TRUE

TypePCA <- 'SPCA'

# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- FALSE

window_size <- 10

nbCPU <- 10
MaxRAM <- 0.5

# number of clusters
nbclusters <- 50

Excluded_WL <- NA

Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = Input_Mask_File,
                                                 Output_Dir = Output_Dir, TypePCA = TypePCA,
                                                 NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                 NIR_Thresh = NIR_Thresh,
                                                 Blue = 450,
                                                 Red = 650,
                                                 NIR = 840)

PCA_Output <- perform_PCA(Input_Image_File = Input_Image_File, 
                          Input_Mask_File = FALSE,
                          Output_Dir = Output_Dir, 
                          TypePCA = TypePCA, 
                          FilterPCA = FilterPCA,
                          nbCPU = nbCPU, 
                          MaxRAM = MaxRAM, 
                          Continuum_Removal = Continuum_Removal)

# path of the raster resulting from dimensionality reduction
PCA_Files <- PCA_Output$PCA_Files
# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

Sel_PC <- select_PCA_components(Input_Image_File = Input_Image_File,
                                Output_Dir = Output_Dir, 
                                PCA_Files = PCA_Output$PCA_Files,
                                TypePCA = PCA_Output$TypePCA, 
                                File_Open = TRUE)

Kmeans_info <- map_spectral_species(Input_Image_File = Input_Image_File,
                                    Input_Mask_File = F,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output, 
                                    nbclusters = nbclusters, 
                                    nbCPU = nbCPU, MaxRAM = MaxRAM)



Index_Alpha <- c('Simpson')

window_size <- 10
map_alpha_div(Input_Image_File = Input_Image_File, 
              Output_Dir = Output_Dir, 
              TypePCA = TypePCA,
              window_size = window_size, 
              nbCPU = nbCPU, 
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha, 
              nbclusters = nbclusters)


con_shan_10 <- rast('biodivmapR/RESULTS/cons_reflectance_combined_image_2/SPCA/ALPHA/Shannon_10')
con_simp_10 <- rast('biodivmapR/RESULTS/cons_reflectance_combined_image_2/SPCA/ALPHA/Simpson_10')

plot(con_shan_10)

writeRaster(con_shan_10, 'biodivmapR/RESULTS/cons_reflectance_combined_image_2/SPCA/ALPHA/Shannon_10.tif', NAflag = 0)
writeRaster(con_simp_10, 'biodivmapR/RESULTS/cons_reflectance_combined_image_2/SPCA/ALPHA/Simpson_10.tif', NAflag = 0)

?writeRaster


# mEd ---------------------------------------------------------------------
# clip conservation raster by conservation polygon ... let's see how we go

library(terra)

cons_raster <- stack('data_out/cons_reflectance_combined_image.tif')
cons_polygon <- st_read('data/cons_shapefile/cons_poly_24.shp') %>%
  st_zm()

??crop

clipped_cons <- terra::mask(cons_raster, cons_polygon)
plot(clipped_cons)

st_area(cons_polygon) #10085.4m2

plot(clipped_cons)

#70832525 of 107476308 NA?
c_spectral_data <- as.data.frame(getValues(clipped_cons)) %>%
  na.omit()

calculate_mean_euclidean_distance <- function(spectral_data) {
  dist_matrix <- as.matrix(dist(spectral_data))
  mean_distance <- mean(dist_matrix)
  return(mean_distance)
}
#  can't do this because it will require 5000000 GB hahhaahhahahahaha 
mean_distance_entire_raster <- calculate_mean_euclidean_distance(c_spectral_data)
print(mean_distance_entire_raster)


# co-efficient of variance metrics ------------------------------------------------

#read in subplot dimensions
subplot_dimensions <- read.csv('data/cons_fishnet_simpson_nd0.csv') %>%
  dplyr::select(c(1:3, 6, 15:20))

subplot_ids <- unlist(lapply(1:5, function(i) paste(i, 1:5, sep="_")))

subplot_dimensions <- subplot_dimensions %>% mutate(subplot_id = subplot_ids)

raster_data <- stack('data_out/cons_reflectance_combined_image.tif')

extract_pixel_values <- function(raster_data, min_x, max_x, min_y, max_y) {
  # define the extent
  extent <- extent(min_x, max_x, min_y, max_y)
  
  # crop the raster to the extent
  cropped_raster <- crop(raster_data, extent)
  
  # extract pixel values
  pixel_values <- getValues(cropped_raster)
  
  # remove NA values (if any)
  pixel_values <- pixel_values[!is.na(pixel_values)]
  
  return(pixel_values)
}

# empty df for results
results <- data.frame(subplot_id = character(), CV = numeric(), stringsAsFactors = FALSE)

# CV function
calculate_cv <- function(values) {
  if (length(values) > 0) {
    return(sd(values) / mean(values))
  } else {
    return(NA)
  }
}

# loop through each subplot
for (i in 1:nrow(subplot_dimensions)) {
  subplot_id <- subplot_dimensions$subplot_id[i]
  min_x <- subplot_dimensions$Min_x_coord[i]
  max_x <- subplot_dimensions$Max_x_coord[i]
  min_y <- subplot_dimensions$Min_y_coord[i]
  max_y <- subplot_dimensions$Max_y_coord[i]
  
  # extract pixel values
  pixel_values <- extract_pixel_values(raster_data, min_x, max_x, min_y, max_y)
  
  # calculate CV
  cv_value <- calculate_cv(pixel_values)
  
  # store the result
  results <- rbind(results, data.frame(subplot_id = subplot_id, CV = cv_value))
}

# Print the results
print(results)

## NOW FOR THE MASKED RASTER 
install.packages('stars')
library(stars)
library(raster)

mask <- read_stars('biodivmapR/RESULTS/cons_reflectance_combined_image/SPCA/ShadeMask_Update')

mask_raster <- as(mask, "Raster")
mask_raster[mask_raster == 0] <- NA

raster_data_masked <- mask(raster_data, mask_raster)

plot(raster_data_masked)

extract_pixel_values_m <- function(raster_data_masked, min_x, max_x, min_y, max_y) {
  # define the extent
  extent <- extent(min_x, max_x, min_y, max_y)
  
  # crop the raster to the extent
  cropped_raster_m <- crop(raster_data_masked, extent)
  
  # extract pixel values
  pixel_values_m <- getValues(cropped_raster_m)
  
  # remove NA values (if any)
  pixel_values_m <- pixel_values_m[!is.na(pixel_values_m)]
  
  return(pixel_values_m)
}

# empty df for results
results_masked <- data.frame(subplot_id = character(), CV = numeric(), stringsAsFactors = FALSE)

# CV function
calculate_cv <- function(values) {
  if (length(values) > 0) {
    return(sd(values) / mean(values) * 100)
  } else {
    return(NA)
  }
}

# loop through each subplot
for (i in 1:nrow(subplot_dimensions)) {
  subplot_id <- subplot_dimensions$subplot_id[i]
  min_x <- subplot_dimensions$Min_x_coord[i]
  max_x <- subplot_dimensions$Max_x_coord[i]
  min_y <- subplot_dimensions$Min_y_coord[i]
  max_y <- subplot_dimensions$Max_y_coord[i]
  
  # extract pixel values
  pixel_values_m <- extract_pixel_values(raster_data_masked, min_x, max_x, min_y, max_y)
  
  # calculate CV
  cv_value <- calculate_cv(pixel_values_m)
  
  # store the result
  results_masked <- rbind(results_masked, data.frame(subplot_id = subplot_id, CV = cv_value))
}

# Print the results
print(results_masked)

# field diversity ---------------------------------------------------------
library(vegan)

survey_data <- read.csv('data/ausplots_march_24.csv')
cons_survey_data <- survey_data %>%
  filter(site_unique == 'NSABHC012')

# Extract only direction of the transect (no numbers)
cons_survey_data$transect_direction <- gsub('[[:digit:]]+', '', cons_survey_data$transect)

# Extract only number of the transect (no direction)
cons_survey_data$transect_number <- as.numeric(gsub(".*?([0-9]+).*", "\\1", cons_survey_data$transect))

# Create variable for fixed transect direction (to order them all transects in the same direction)
cons_survey_data$transect_direction2 <- NA 

# Create variable for fixed point number (inverse in some cases as if they had been collected in the same direction)
cons_survey_data$point_number2 <- NA 

# Create XY empty variables for plot XY coordinates
cons_survey_data$X_plot <- NA 
cons_survey_data$Y_plot <- NA

# For loop to homogenize transects and numbers. It converts all E-W to W-E and all S-N to N-S
for (i in 1:nrow(cons_survey_data)){
  if (cons_survey_data[i, "transect_direction"] == "E-W") {
    cons_survey_data[i, "point_number2"] <- 100 - cons_survey_data[i, "point_number"] # If transect E-W, transect fixed is W-E and inverse numbers
    cons_survey_data[i, "transect_direction2"] <- "W-E"
  }
  if (cons_survey_data[i, "transect_direction"] == "W-E") {
    cons_survey_data[i, "point_number2"] <- cons_survey_data[i, "point_number"] # If transect W-E, all stays the same
    cons_survey_data[i, "transect_direction2"] <- "W-E"
  }
  if (cons_survey_data[i, "transect_direction"] == "N-S") {
    cons_survey_data[i, "point_number2"] <- cons_survey_data[i, "point_number"] # If transect N-S, all stays the same
    cons_survey_data[i, "transect_direction2"] <- "N-S"
  }
  if (cons_survey_data[i, "transect_direction"] == "S-N") {
    cons_survey_data[i, "point_number2"] <- 100 - cons_survey_data[i, "point_number"] # If transect S-N, transect fixed is N-S and inverse numbers
    cons_survey_data[i, "transect_direction2"] <- "N-S"
  }
}

# For loop to assign plotXY coordinates to each point intercept
for (i in 1:nrow(cons_survey_data)){
  if (cons_survey_data[i, "transect_direction2"] == "W-E") {
    if (cons_survey_data[i, "transect_number"] == 1){
      cons_survey_data[i, "Y_plot"] <- 10
      cons_survey_data[i, "X_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 2){
      cons_survey_data[i, "Y_plot"] <- 30
      cons_survey_data[i, "X_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 3){
      cons_survey_data[i, "Y_plot"] <- 50
      cons_survey_data[i, "X_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 4){
      cons_survey_data[i, "Y_plot"] <- 70
      cons_survey_data[i, "X_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 5){
      cons_survey_data[i, "Y_plot"] <- 90
      cons_survey_data[i, "X_plot"] <- cons_survey_data[i, "point_number2"]
    }
  }
  if (cons_survey_data[i, "transect_direction2"] == "N-S") {
    if (cons_survey_data[i, "transect_number"] == 1){
      cons_survey_data[i, "X_plot"] <- 10
      cons_survey_data[i, "Y_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 2){
      cons_survey_data[i, "X_plot"] <- 30
      cons_survey_data[i, "Y_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 3){
      cons_survey_data[i, "X_plot"] <- 50
      cons_survey_data[i, "Y_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 4){
      cons_survey_data[i, "X_plot"] <- 70
      cons_survey_data[i, "Y_plot"] <- cons_survey_data[i, "point_number2"]
    }
    if (cons_survey_data[i, "transect_number"] == 5){
      cons_survey_data[i, "X_plot"] <- 90
      cons_survey_data[i, "Y_plot"] <- cons_survey_data[i, "point_number2"]
    }
  }
}

cons_survey_data %>%
  drop_na(standardised_name) %>%
  ggplot(aes(x = X_plot, y = Y_plot, color = growth_form)) +
  geom_point() +
  theme_classic()

# subplot rows and columns - +1 ensures 0 point values fall into correct subplot, 
# pmin ensures 100 point values falls in correct subplot given +1

cons_survey_data$subplot_row <- pmin(ceiling((cons_survey_data$Y_plot + 1) / 20), 5)
cons_survey_data$subplot_col <- pmin(ceiling((cons_survey_data$X_plot + 1) / 20), 5)

# single ID for subplot row and column
cons_survey_data$subplot_id <- paste(cons_survey_data$subplot_row, cons_survey_data$subplot_col, sep = "_")

subplot_diversity <- cons_survey_data %>%
  drop_na(standardised_name) %>%
  group_by(subplot_id) %>%
  summarise(species_richness = n_distinct(standardised_name))

community_matrix <- cons_survey_data %>%
  drop_na(standardised_name) %>%
  count(subplot_id, standardised_name) %>%
  spread(standardised_name, n, fill = 0)

#column V1? what is zis??? 
if ("V1" %in% colnames(community_matrix)) {
  community_matrix <- community_matrix %>% 
    dplyr::select(-V1)
}

# calculate diversity indices
shannon_diversity <- diversity(community_matrix[, -1], index = "shannon")
simpson_diversity <- diversity(community_matrix[, -1], index = "simpson")

#read in biodivmapR values from arcgis table
cons_simpson_biodivmapR <- read.csv('data/cons_fishnet_simpson_nd0.csv') %>%
  dplyr::select(10)

cons_shannon_biodivmapR <- read.csv('data/shannon_nd0_cons.csv') %>%
  dplyr::select(7)

#read in biodivmapR values from arcgis table for NON MASKED raster
cons_nm_simpson_biodivmapR <- read.csv('data/simpson_not_masked_cons.csv') %>%
  dplyr::select(16)

cons_nm_shannon_biodivmapR <- read.csv('data/shannon_not_masked_cons.csv') %>%
  dplyr::select(16)


# need to run the CV code which is further down for this to work! 
# (the cv and biodivmapRSimpson lines)
# tidy up your code bishhh

subplot_diversity <- subplot_diversity %>%
  mutate(
    shannon_diversity = shannon_diversity,
    simpson_diversity = simpson_diversity,
    cv = results$CV,
    cv_masked = results_masked$CV,
    shannon_biodivmapR = cons_shannon_biodivmapR$MEAN,
    simpson_biodivmapR = cons_simpson_biodivmapR$MEAN,
    shannon_non_mask_biodivmapR = cons_nm_shannon_biodivmapR$MEAN,
    simspon_non_mask_biodivmapR = cons_nm_simpson_biodivmapR$MEAN
  )

subplot_coordinates <- cons_survey_data %>%
  dplyr::select(subplot_id, X_plot, Y_plot) %>%
  distinct()

# merge coords 
subplot_diversity_c <- subplot_diversity %>%
  left_join(subplot_coordinates, by = "subplot_id")


# shannon's diversity
ggplot(subplot_diversity_c, aes(x = X_plot, y = Y_plot)) +
  geom_tile(aes(fill = shannon_diversity), color = "white") +
  scale_fill_viridis_c() +
  theme_classic() +
  labs(fill = "Shannon's Diversity")

# simpson's diversity
ggplot(subplot_diversity_c, aes(x = X_plot, y = Y_plot)) +
  geom_tile(aes(fill = simpson_diversity), color = "white") +
  scale_fill_viridis_c() +
  theme_classic() +
  labs(fill = "Simpson's Diversity")


df_long <- subplot_diversity %>%
  dplyr::select(shannon_diversity, shannon_biodivmapR, shannon_non_mask_biodivmapR) %>%
  pivot_longer(cols = c(shannon_biodivmapR, shannon_non_mask_biodivmapR), 
               names_to = "variable", 
               values_to = "value")


#shannons compared
ggplot(df_long, aes(x = shannon_diversity, y = value, color = variable)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(x = "Shannon Diversity", y = "Value") +
  theme_minimal()


# biodivmapR values -------------------------------------------------------

ggplot(data = subplot_diversity, aes(x = shannon_diversity, y = shannon_biodivmapR)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(
    title = "Shannon Diversity vs. Shannon BiodivmapR",
    x = "Shannon Diversity - field",
    y = "Shannon Diversity - biodivmapR"
  ) +
  theme_minimal()

model <- lm(shannon_biodivmapR ~ shannon_diversity, data = subplot_diversity)
summary(model)

cor_test <- cor.test(subplot_diversity$shannon_diversity, subplot_diversity$shannon_biodivmapR)
cor_test$estimate


cv_long <- subplot_diversity %>%
  select(shannon_diversity, cv, cv_masked) %>%
  pivot_longer(cols = c(cv, cv_masked), 
               names_to = "variable", 
               values_to = "value")


ggplot(cv_long, aes(x = shannon_diversity, y = value, color = variable, group = variable)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F) + 
  labs(x = "Shannon Diversity", y = "CV value") +
  theme_minimal()

ggplot(data = subplot_diversity, aes(x = shannon_diversity, y = cv)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = "Shannon Diversity vs. CV",
    x = "Shannon Diversity - field",
    y = "Co-efficient of variance"
  ) +
  theme_minimal()

model <- lm(cv ~ shannon_diversity, data = subplot_diversity)
summary(model)

cor_test <- cor.test(subplot_diversity$shannon_diversity, subplot_diversity$cv)
cor_test$estimate


# convex hull volume ------------------------------------------------------

# from crofts - https://github.com/ALCrofts/CABO_SVH_Forest_Sites/tree/v1.0

calc.CHV <- function(spectral_df, # Dataframe with spectral reflectance values
                     areas_of_interest, # What you want to calculate spectral diversity for, ie. grouping variable, here it is plots. 
                     wavelengths, # Cols where spectral reflectance values are
                     rarefraction, # If TRUE, spectral observations are standardized and randomly resampled. If FALSE, uses all spectral observations as is. 
                     n # Number of random resampling events, if rarefraction = T.
){  
  
  if(rarefraction == TRUE){
    
    # a) Determine min number of spectral points/pixels observed in a plot  
    spectral_points <- spectral_df %>% 
      group_by({{areas_of_interest}}) %>%
      summarise(points = n())
    
    min_points = min(spectral_points$points)
    
    # 1b) Ordinate the data  
    PCA <- spectral_df %>%
      select(c({{wavelengths}})) %>%
      rda(scale = F) 
    
    # 2b) Create a dataframe of the first 3 PC axes
    PCA_scores <- data.frame(PCA$CA$u) %>%
      select(c(1:3)) %>%
      cbind({{spectral_df}}) %>%
      select(c('PC1', 'PC2', 'PC3', {{areas_of_interest}}))
    
    # c) Convert to quoted areas_of_interest object, needed to name list objects created in next step.
    areas_of_interest2 <- deparse(substitute(areas_of_interest))
    
    # d) calculate average CHV for each plot across all resampling events
    CHV <- PCA_scores %>%
      group_split({{areas_of_interest}}) %>% # split by plot
      set_names(map(., ~unique(.[[areas_of_interest2]]))) %>% # assign plot names to list objects 
      map(~replicate(n,
                     .x %>%
                       group_by({{areas_of_interest}}) %>%
                       sample_n(min_points, replace = F), F) %>% # resample PC scores n-times and retain the minimum number of observed spectral observations, where resampling occurs without replacement  
            map(~convhulln(.x[-4], option = 'FA')) %>% # fit convex hulls to each plot and each resampling event
            map_dbl('vol')) %>% # retain the convex hull volume for each plot and each resampling event
      map_df(~mean(.x)) %>% # calculate mean CHV for each plot across all resampling events
      do.call(rbind, .) %>% # collapse to dataframe of plots and CHV
      as.data.frame() %>%
      rownames_to_column('plot_field_id') %>%
      rename(CHV = V1)
    
    rm(areas_of_interest2)
    
    return(CHV)
  }
  
  if(rarefraction == FALSE){
    
    areas_of_interest2 <- deparse(substitute(areas_of_interest))
    
    PCA2 <- spectral_df %>%
      select(c({{wavelengths}})) %>%
      rda(scale = F) 
    
    CHV2 <- data.frame(PCA2$CA$u) %>%
      select(c(1:3)) %>%
      cbind({{spectral_df}}) %>%
      select(c('PC1', 'PC2', 'PC3', {{areas_of_interest}})) %>%
      group_split({{areas_of_interest}}) %>%
      set_names(map(., ~unique(.[[areas_of_interest2]]))) %>%
      map(~convhulln(.x[-4], option = 'FA')) %>%
      map_dbl('vol') %>%
      as.data.frame() %>%
      rownames_to_column('plot_field_id') %>%
      rename(CHV = '.')
    
    rm(areas_of_interest2) 
    
    return(CHV2)
  }
  
}

names(raster_data) <- c("blue", "green", "red", "red_edge", "nir")

c_spectral_data <- as.data.frame(getValues(raster_data)) %>%
  na.omit()

raster_coords <- coordinates(raster_data)
raster_values <- getValues(raster_data)

c_spectral_data <- as.data.frame(cbind(raster_coords, raster_values)) %>%
  na.omit()








## ALL OF THESE TAKE A REALLY REALLY LONG TIME

# approach 1

# determine which subplot a point belongs to
get_subplot_id <- function(x, y, subplot_df) {
  subplot <- subplot_df %>%
    filter(x >= Min_x_coord & x <= Max_x_coord & y >= Min_y_coord & y <= Max_y_coord)
  
  if (nrow(subplot) == 1) {
    return(subplot$subplot_id)
  } else {
    return(NA)  # Return NA if the point does not belong to any subplot
  }
}

# apply to each pixel in c_spectral_data
c_spectral_data$subplot_id <- mapply(get_subplot_id, c_spectral_data$x, c_spectral_data$y, MoreArgs = list(subplot_df = subplot_dimensions))

# approach 2
# Convert raster data to sf object
c_spectral_data_sf <- st_as_sf(c_spectral_data, coords = c("x", "y"), crs = st_crs(raster_data))

# Function to create closed polygons
create_polygon <- function(min_x, min_y, max_x, max_y) {
  st_polygon(list(matrix(c(min_x, min_y, 
                           max_x, min_y, 
                           max_x, max_y, 
                           min_x, max_y, 
                           min_x, min_y), # Closing the polygon
                         ncol = 2, byrow = TRUE)))
}

# Convert subplot dimensions to sf object
subplot_dimensions_sf <- subplot_dimensions %>%
  rowwise() %>%
  mutate(geometry = st_sfc(create_polygon(Min_x_coord, Min_y_coord, Max_x_coord, Max_y_coord))) %>%
  st_as_sf(crs = st_crs(raster_data))

# Perform spatial join
c_spectral_data_sf <- st_join(c_spectral_data_sf, subplot_dimensions_sf, join = st_within)

# Convert back to dataframe if needed
c_spectral_data_2 <- as.data.frame(c_spectral_data_sf)

# to do list
# apply subplot ids to the spectral values data frame 
# create chv formula - DONE (copied from crofts - need to do above to be able to use it)
# consider sv formula?
# mEd method? 
# compare simp, shan for unmasked stuff - DONE
# pielou's evenness?
# functional diversity traits? which ones - think, causal, foliar? 
# check austraits for said traits

