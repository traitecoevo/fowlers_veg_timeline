Information prepared by Carly Steen. 9th January 2024.

Note: plot/site used interchangeably = 1ha monitoring plot, surveyed using AusPlots methodology.
______________________________________________________________________________
## Fowlers_Gap_summary_site_data.xlsx ##

Site Location tab: overview of each plot at Fowlers Gap, location summary, SW corner coordinates, site comments

Collection at plot tab: quick overview at what has been collected and when at each plot at Fowlers Gap

Plot corners tab: location coordinates for each plot, includes each corner, centre point and transect start/finish

Structure description tab: vegetation structural description at each plot and plot visit. Note, 2022 revisits do not have full species identification yet.

Point intercept tab: percentage cover (count/1010*100) for each species within the plot based on point intercept observations.

Growth Form Data tab: percentage cover (count/1010*100) of each growth form within the plot based on point intercept observations.

Substrate tab: percentage cover (count/1010*100) of each substrate type within the plot based on point intercept observations.

Plant list: list of plant species vouchered and identified at each plot. Awaiting herbarium determination for the 2022 revisit.

Soil characterisation tab: Soil pit observations and Australian Soil Classification at a plot.

Bulk Density tab: Bulk density results for each plot, samples taken at soil pit.

______________________________________________________________________________
## Fowlers Gap_PI_GrowthForm_data-ALL.xlsx ##

Point intercept data for all plots at Fowlers Gap, including revisits. Growth form data.
This data includes unpublished sites (2022 revisits), growth form and veg barcode have only been supplied.  The veg barcode can be used in the future when herbarium determinations have been returned to include species name if required.  I have included the field names for the 2022 revisits, as discussed, this may help with growth form interpretation.

Attributes:

site_location_name - name of a plot, follows AusPlots alpha-numeric naming convention STATE (2 letters) PLOT TYPE (1 letter) BIROEGION CODE (3 letters) PLOT NUMBER (4 numbers)

site_location_id - unique numeric id for each site location (plot)

site_location_visit_id - unique numeric id for each site location (plot) visit (allows for identifying multiple visits to the one plot).

visit_start_date - date visited / plot observation

transect - name for each of the transects within a plot, spaced in 10m intervals running E-W or S-N.  The transect name reflects the direction the point intercept was taken i.e. the observer may walk it E-W or W-E, the numeric value denotes which transect the observer is on. https://linkeddata.tern.org.au/viewer/ausplots/id/http://linked.data.gov.au/def/ausplots-cv/c5a32483-bf2f-421d-b03d-6d81e1195de2

point_number - 0-100. each point where a point intercept observation is taken along the transect.

growth_form - form or shape of individual plant or Australian broad floristic land cover types, LUT: https://linkeddata.tern.org.au/viewer/ausplots/id/http://linked.data.gov.au/def/ausplots-cv/eae155c7-669c-463a-8d01-01b090472732

height - metres. Height of the plant described as the distance between the soil and the plant canopy, recorded at the uppermost height at which the laser beam intersects with a plant.

veg_barcode - vegetation vouchering, unique barcode assigned to each specimen collected at a plot. Can be used as a primary key to link to herbarium determiniation.

herbarium_determination - The species name as identified by a herbarium botanist. Each specimen collected in the field is assigned a 'field name', the specimen is then sent to herbarium to be identified. 

dead - TRUE/FALSE. All vegetation ‘hits' on the point intercept transect will be classified as either alive or dead.

in_canopy_sky - TRUE/FALSE. If foliage or branches are not intercepted, it is recorded as "in-canopy sky". If foliage or branches are intercepted, the species is recorded as a canopy hit and an estimate of the height is recorded. https://linkeddata.tern.org.au/viewer/ausplots/id/http://linked.data.gov.au/def/ausplots-cv/bc8156c2-c2a7-4b2b-8ece-3f1959734d6e

field_name - name given to uniquely identify veg voucher in the field, if species is unknown to observer it may be just a descriptive summary of the veg voucher e.g. pale purple daisy.  Field names are not generally published as they are internal identifiers which are then replaced with the herbarium botanists determination.  Only the 2022 unpublished site data field names have been included here.

______________________________________________________________________________

## TERNSites_visits_peryear.xlsx ##

Lists all TERN Ecosystem Surveillance plots as of December 2023. Plots where there has been limited data collected have been removed from this list.  There are plots on the list that are not publically available through AusplotsR or EcoPlots.  
Plots not publically available are those that:
(1) field data entry has not been complete.
(2) field data entry has been complete, but awaiting complete species identification by herbarium (i.e. herbarium determination).
These plots have been identified in the spreadsheet, highlighted with colour.  Colour key has been included in the spreadsheet.

Attributes:
site_location_id - unique id for each plot.
Site Name - name of a plot, follows AusPlots alpha-numeric naming convention STATE (2 letters) PLOT TYPE (1 letter) BIROEGION CODE (3 letters) PLOT NUMBER (4 numbers)
Visit 1 / Visit 2 / Visit 3 / Visit 4 / Visit 5 - Identifies what year a plot has been visited and each year it has subsequently been visited (where revisit has occurred)

______________________________________________________________________________
## TERN Contacts: ##

Ben Sparrow
ASSOCIATE PROFESSOR AND PROGRAM LEAD, TERN ECOSYSTEM SURVEILLANCE 
P +61 (0) 8 8313 1201
ben.sparrow@adelaide.edu.au

Carly Steen
DATA AND CURATION, TERN ECOSYSTEM SURVEILLANCE
M 0402 028 665 
carly.steen@adelaide.edu.au

www.tern.org.au
