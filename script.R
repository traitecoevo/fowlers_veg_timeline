
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
folders <- folders[folders != ""]

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
  output_filename <- file.path("data_out", paste0(folder, "_combined_image.tif"))
  writeRaster(combined_image, filename = output_filename, format = "GTiff", overwrite = TRUE)
  
  plot(combined_image)
}

emu_exc_raster <- raster('data_out/emu_exc_reflectance_combined_image.tif')
emu_raster <- raster('data_out/emu_reflectance_combined_image.tif')
emu_exc_blue <- raster('data/2024_mosaics/emu_exc_reflectance/blue.tif')
emu_blue <- raster('data/2024_mosaics/emu_reflectance/blue.tif')


# biodivmapR --------------------------------------------------------------

# load biodivMapR and useful libraries 
remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR', force = T)
install.packages("sf")
library(sf)
library(stars)
library(utils)

tmpdir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/biodivmapR'
NameRaster <- 'cons_reflectance_combined_image.tif'
destfiletiff <- file.path(tmpdir,NameRaster)

BandName <- c('Band_1', 'Band_2', 'Band_3', 'Band_4', 'Band_5')
SpectralBands <- c(450, 560, 650, 730, 840)
WLunits <- 'Nanometers'

## CANT DO THIS - no template for phantom 4. need to create this.
create_hdr(ImPath = destfiletiff, Sensor = 'Phantom 4', 
           SpectralBands = SpectralBands,BandName = BandName, WLunits = WLunits)

Input_Image_File <- destfiletiff

Output_Dir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/biodivmapR/RESULTS'

