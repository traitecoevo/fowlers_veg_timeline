# load biodivMapR and useful libraries 
remotes::install_github('cran/dissUtils')
remotes::install_github('jbferet/biodivMapR', force = T)

install.packages("sf")
library(sf)
library(stars)
library(utils)


# url for the S2 subset
url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset'
# create a temporary directory (choose your own data directory)
tmpdir <- tempdir()
# name your binary raster with the same name as the online file
NameRaster <- 'S2A_T33NUD_20180104_Subset'
destfile <- file.path(tmpdir,NameRaster)
download.file(url = url, destfile = destfile, method = 'auto', quiet = FALSE, mode = "wb")

# url for the S2 subset header
urlhdr <-  'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/RASTER/S2A_T33NUD_20180104_Subset.hdr'
# name your raster HDR with the same name as the binary raster, with .hdr extension
destfile_HDR <- get_HDR_name(destfile,showWarnings = FALSE)
download.file(url = urlhdr, destfile = destfile_HDR, method = 'auto', quiet = FALSE, mode = "w")



# read ENVI file with stars
Stars_S2 <- stars::read_stars(destfile, along = 'band',proxy = FALSE)
# write it as a tiff image
# create a specific directory for the tiff image and name your raster
desttiff <- file.path(tmpdir,'TIFF',fsep = '\\')
dir.create(desttiff,showWarnings = FALSE)
destfiletiff <- file.path(desttiff,'S2_Subset.tif',fsep = '\\')
r <- write_stars(Stars_S2, dsn=destfiletiff, driver =  'GTiff', type='Int16')

# read ENVI file with stars
create_hdr(ImPath = destfiletiff, Sensor = 'SENTINEL_2A', 
           SpectralBands = NULL, BandName = NULL, WLunits = NULL)

# library
library(zip)
# name zip file including plots located on the tile
destzip <- file.path(tmpdir,'S2A_T33NUD_Plots.zip',fsep = '\\')
# url for the zip file
url <- 'https://gitlab.com/jbferet/myshareddata/-/raw/master/biodivMapR_S2_Sample/VECTOR/S2A_T33NUD_Plots.zip'
download.file(url = url, destfile = destzip)
destunz <- file.path(tmpdir,'S2A_T33NUD_Plots',fsep = '\\')
unzip(zipfile = destzip,exdir = destunz)

Input_Image_File <- destfiletiff

# Set to FALSE if no mask available
Input_Mask_File <- FALSE

Output_Dir <- 'C:/Users/adele/Documents/fowlers_veg_timeline/data_out/RESULTS'

NDVI_Thresh <- 0.8
Blue_Thresh <- 500
NIR_Thresh <- 1500

Continuum_Removal <- TRUE
TypePCA <- 'SPCA'

# PCA FILTERING:        Set to TRUE if you want second filtering based on PCA outliers to be processed.
# Slower process
# Automatically set to FALSE if TypePCA     = 'MNF'
FilterPCA <- FALSE

window_size <- 10

nbCPU <- 4
MaxRAM <- 0.5
nbclusters <- 50


Input_Mask_File <- perform_radiometric_filtering(Image_Path = Input_Image_File, Mask_Path = Input_Mask_File,
                                                 Output_Dir = Output_Dir, TypePCA = TypePCA,
                                                 NDVI_Thresh = NDVI_Thresh, Blue_Thresh = Blue_Thresh,
                                                 NIR_Thresh = NIR_Thresh)

Excluded_WL <- c(0, 400)
Excluded_WL <- rbind(Excluded_WL, c(895, 1005))
Excluded_WL <- rbind(Excluded_WL, c(1180, 1480))
Excluded_WL <- rbind(Excluded_WL, c(1780, 2040))


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
