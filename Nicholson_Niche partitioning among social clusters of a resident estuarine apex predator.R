
#############################################################################################

# Niche partitioning among social clusters of an estuarine resident apex predator
# Krista Nicholson, Neil Loneragan and Lars Bejder
# Published in Behavioral Ecology and Sociobiology 2021

#############################################################################################
#####################    Stable isotope analyses - Mixing models    #########################
#############################################################################################


# Set working directory
setwd("C:/")

#Download package: MixSIAR

library(MixSIAR)

# Load consumer, sources and discrimination factor data

mix.filename <- "consumers21.csv"
source.filename <- "sources_raw21_FG.csv"
discr.filename <- "TDF21_FG.csv"

# Models
# Model 1: NULL
# Model 2: group (=social cluster), random effect
# Model 3: sex (male, female), fixed effect
# Model 4: age (2 classes: juvenile, adult), fixed effect
# Model 5: sex + age, fixed effects
# Model 6: sex_age (created a new factor = sex * age), fixed effect
# Model 7: sex_age + group, fixed and random effect, respectively
# Model 8: group + ID, random effects, ID nested within group
# Model 9: group + age, random and fixed effect, respectively/ Age nested within group
# Model 10: group + sex, random and fixed effect, respectively, Sex nested within group


# Model 1
mix1 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors=NULL,
                      fac_random=NULL,
                      fac_nested=NULL,
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix1)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix1)

plot_data(filename="isospace_plot1", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix1,source,discr)

if(mix1$n.iso==2) calc_area(source=source,mix=mix1,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_1.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix1, source)

jags.1 <- run_model(run="very long",mix1,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.1, mix1, source)


# Model 2
mix2 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors="group",
                      fac_random=TRUE,
                      fac_nested=FALSE,
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix2)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix2)

plot_data(filename="isospace_plot2", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix2,source,discr)

if(mix2$n.iso==2) calc_area(source=source,mix=mix2,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_2.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix2, source)

jags.2 <- run_model(run="very long",mix2,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.2, mix2, source)



# Model 3

mix3 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors="sex",
                      fac_random=FALSE,
                      fac_nested=FALSE,
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix3)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix3)

plot_data(filename="isospace_plot3", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix3,source,discr)

if(mix3$n.iso==2) calc_area(source=source,mix=mix3,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_3.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix3, source)

jags.3 <- run_model(run="very long",mix3,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.3, mix3, source)


# Model 4

mix4 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors="age",
                      fac_random=FALSE,
                      fac_nested=FALSE,
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix4)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix4)

plot_data(filename="isospace_plot4", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix4,source,discr)

if(mix4$n.iso==2) calc_area(source=source,mix=mix4,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_4.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix4, source)

jags.4 <- run_model(run="very long",mix4,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.4, mix4, source)


# Model 5
mix5 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors=c("age","sex"),
                      fac_random=c(FALSE,FALSE),
                      fac_nested=c(FALSE,FALSE),
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix5)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix5)

plot_data(filename="isospace_plot5", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix5,source,discr)

if(mix5$n.iso==2) calc_area(source=source,mix=mix5,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_5.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix5, source)

jags.5 <- run_model(run="very long",mix5,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.5, mix5, source)


# Model 6
mix6 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors="sex_age",
                      fac_random=FALSE,
                      fac_nested=FALSE,
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix6)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix6)


plot_data(filename="isospace_plot6", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix6,source,discr)

if(mix6$n.iso==2) calc_area(source=source,mix=mix6,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_6.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix6, source)

jags.6 <- run_model(run="very long",mix6,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.6, mix6, source)


# Model 7

mix7 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors=c("group","sex_age"),
                      fac_random=c(TRUE,FALSE),
                      fac_nested=c(FALSE,TRUE),
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix7)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix7)


plot_data(filename="isospace_plot7", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix7,source,discr)

if(mix7$n.iso==2) calc_area(source=source,mix=mix7,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_7.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix7, source)

jags.7 <- run_model(run="very long",mix7,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.7, mix7, source)


# Model 8
mix8 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors=c("group","ID"),
                      fac_random=c(TRUE,TRUE),
                      fac_nested=c(FALSE,TRUE),
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix8)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix8)

plot_data(filename="isospace_plot8", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix8,source,discr)

if(mix8$n.iso==2) calc_area(source=source,mix=mix8,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_8.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix8, source)

jags.8 <- run_model(run="very long",mix8,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.8, mix8, source)


# Model 9

mix9 <- load_mix_data(filename=mix.filename,
                      iso_names=c("d13C","d15N"),
                      factors=c("group","age"),
                      fac_random=c(TRUE,FALSE),
                      fac_nested=c(FALSE,TRUE),
                      cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix9)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix9)


plot_data(filename="isospace_plot9", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix9,source,discr)

if(mix9$n.iso==2) calc_area(source=source,mix=mix9,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_9.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix9, source)

jags.9 <- run_model(run="very long",mix9,source,discr,model_filename,
                    alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.9, mix9, source)


# Model 10

mix10 <- load_mix_data(filename=mix.filename,
                       iso_names=c("d13C","d15N"),
                       factors=c("group","sex"),
                       fac_random=c(TRUE,FALSE),
                       fac_nested=c(FALSE,TRUE),
                       cont_effects=NULL)

source <- load_source_data(filename = "sources_raw21_FG.csv",
                           source_factors = NULL,
                           conc_dep=FALSE,
                           data_type = "raw", mix10)

discr <- load_discr_data(filename = "TDF21_FG.csv",mix10)


plot_data(filename="isospace_plot10", plot_save_pdf=TRUE,
          plot_save_png=FALSE, mix10,source,discr)

if(mix10$n.iso==2) calc_area(source=source,mix=mix10,discr=discr)

plot_prior(alpha.prior=1,source)

model_filename <- "MixSIAR_model_10.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix10, source)

jags.10 <- run_model(run="very long",mix10,source,discr,model_filename,
                     alpha.prior = 1,resid_err,process_err)

output_JAGS(jags.10, mix10, source)


#############################################################################################
##########    Activity space analyses and utilization distribtion overlaps    ###############
#############################################################################################

# Download packages: 'OpenStreetMap', 'adehabitatHR', 'ggmap' and 'XML', 'maps'

library(sp)
library(OpenStreetMap)
library(rgdal)
library(ggplot2)
library(ggmap)
library(XML)
library(adehabitatHR)
library(raster)
library(dismo)
library(XML)
library(rgeos)

# Set working directory
setwd("C://")

# Input locations of individuals
PH <- read.csv("PH_homeranges_070520_GIS_UTM.csv",h=T)
PH <- subset(PH,select = c(Long,Lat,SC))
class (PH)
head(PH)


# Extract social cluster names from PH
IDs <- sort(as.vector(unique(PH[,"SC"])))
length(IDs)
dfID <- data.frame(IDs)

# Input study area
PH_studyarea <- readOGR(dsn = "C:/", layer = "PH_Study_Area_UTM" ) 
PH_studyarea_UTM <- spTransform(PH_studyarea, CRS("+init=epsg:32750"))
plot(PH_studyarea)


# Create a boolean raster file of the study area (water=1,land=0)
rgrid <- raster(extent(PH_studyarea))
res(rgrid) <- c(100, 100)                # set grid resolution
rgrid[] <- 1                             # assign a value of 1 to each grid cell
rgrid_msk <- mask(rgrid,PH_studyarea)    # clip the grid layer with the shape file
rgrid_msk[is.na(rgrid_msk)] <- 0         # assign a 0 to all cells that are not water
plot(rgrid_msk)
proj4string(rgrid_msk) <- CRS(proj4string(PH_studyarea_UTM)) ## set layer CRS to UTM zone 50 south


# Convert mask to spatial points data frame
grid_ae <- as(rgrid_msk, 'SpatialPointsDataFrame')
grid_ae <- grid_ae[!is.na(grid_ae@data$layer), ]
gridded(grid_ae) <- TRUE
summary(grid_ae)
hab <- grid_ae                          # assign to a new object hab

# Convert Lat and Long to a spatialpixeldataframe
xy_PH = PH[c("Long", "Lat")]
coordinates(xy_PH)=c("Long","Lat")
PH_sp<-SpatialPointsDataFrame(xy_PH, PH)
proj4string(PH_sp) <- CRS("+init=epsg:32750") # set coordinate system as UTM50

# Run kernel density estimates using the habitat as grid. 
# Epanechnikov kernel and smoothing factor 'h' set as href
ud_epa <- kernelUD(PH_sp[,3],
                   h="href",
                   grid = grid_ae,
                   kern="epa",
                   extent = 2)


# Assign to new object
ud_epa_new <- ud_epa


# Check home ranges for both 95% and 50% utilization distributions
homerange <- getverticeshr(ud_epa_new[["XX"]],percent = 95)
plot(homerange,col = 1:3)

homerange <- getverticeshr(ud_epa_new[["XX"]],percent = 50)
plot(homerange,col = 1:3)


# Change to spatial pixels data frame
udspdf <- estUDm2spixdf(ud_epa_new)
fullgrid(udspdf) <- TRUE
fullgrid(hab)<-TRUE


# Multiply each UD with the 1/0 (hab) and rescale so that the sum of the new UD sums up to 0.00001 
resu <- lapply(1:ncol(udspdf), function(i) {udspdf[[i]] * hab[[1]]/sum(udspdf[[i]] * hab[[1]])/10000}) 
resu <- as.data.frame(resu)
names(resu) <- names(udspdf@data)
udspdf@data <- resu
fullgrid(udspdf) <- FALSE

# Transfer back into a object of class estUDm
re <- lapply(1:ncol(udspdf), function(i) { 
  so <- new("estUD", udspdf[,i]) 
  so@h <- list(h=0, meth="specified")
  so@vol <- FALSE 
  return(so) 
}) 

names(re) <- names(udspdf)
class(re) <- "estUDm" 


# Save object
save(re, file="Kernel_densities_epa_first five.RData")
load("Kernel_densities_epa_first five.RData")

# Calculate home range overlaps for 95% and 50% KDEs
overlaps_UDOI_epa <- kerneloverlaphr(re, method="PHR", percent=50, conditional=TRUE) 

# Write objects as csv files (repeat for 50% KDEs)
write.csv(overlaps_UDOI_epa, file="overlaps_epa_PHR_SocialClusters_95.csv")



