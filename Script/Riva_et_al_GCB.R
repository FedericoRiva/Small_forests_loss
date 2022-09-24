#Riva, Martin, Millard and Fahrig. 2022, Loss of the world's smallest forests. Global Change Biology

##################################################
### 0 - SESSION INFO
################################################## 

# R version 4.1.0 (2021-05-18)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_Canada.1252  LC_CTYPE=English_Canada.1252    LC_MONETARY=English_Canada.1252
# [4] LC_NUMERIC=C                    LC_TIME=English_Canada.1252    
# 
# attached base packages:
#   [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] fasterize_1.0.3        sf_1.0-6               landscapemetrics_1.5.4 dplyr_1.0.8           
# [5] rgdal_1.5-28           raster_3.5-15          sp_1.4-6              
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.8         magrittr_2.0.2     units_0.8-0        tidyselect_1.1.2   lattice_0.20-44   
# [6] R6_2.5.1           rlang_1.0.1        fansi_1.0.2        tools_4.1.0        grid_4.1.0        
# [11] KernSmooth_2.23-20 utf8_1.2.2         e1071_1.7-9        terra_1.5-21       cli_3.2.0         
# [16] DBI_1.1.2          class_7.3-19       ellipsis_0.3.2     assertthat_0.2.1   tibble_3.1.6      
# [21] lifecycle_1.0.1    crayon_1.5.0       purrr_0.3.4        vctrs_0.3.8        codetools_0.2-18  
# [26] glue_1.6.2         proxy_0.4-26       compiler_4.1.0     pillar_1.7.0       generics_0.1.2    
# [31] classInt_0.4-3     pkgconfig_2.0.3  


##################################################
### 1 - PREPARATION
################################################## 

# open R packages
library(raster)
library(rgdal)
library(dplyr)
library(sp)
library(landscapemetrics)
library(sf)
library(fasterize)
library(parallel)
library(glmmTMB)
library(sjPlot)
library(effects)
library(ggeffects)
library(ggplot2)
library(ggthemes)
library(ggpubr)


#
# list of dataset necessary available in supplementary information with the paper
#

##################################################
### 2 - OPEN FILES PREPARED FOR ANALYSIS (code for data preparation at the end of the script)
################################################## 

# reclassified raster for analysis
lc1992 <- raster("lc1992_waterna.tif")
lc2020 <- raster("lc2020_waterna.tif")
biome_rast <- raster("biome_rast30.tif")

# random sample of forest cells 
for1992 <- read.csv("data_forest1992.csv")
for2020 <- read.csv("data_forest2020.csv")

# area of the forest in which each randomly sampled point fell
res1992 <- read.csv("results1992.csv")
res2020 <- read.csv("results2020.csv")

# final dataset for models
data <- read.csv("results_final.csv")

# aggregated raster at 30km for map in fig. 1
raster_plot1 <- raster("agg1992_30.tif")
raster_plot2 <- raster("agg2020_30.tif")

# open raster ESA-CCI-LC reclassified as human land cover vs not
pop <- raster("lc2020_humandominated.tif")

##################################################
### 3 - CALCULATE STATISTICS FOR FOREST CHANGE
##################################################

# calculate number of forest cells between 1992 and 2020
cellStats(lc1992, sum) # 492440932
cellStats(lc2020, sum) # 487172235

# each cell is 300*300m, i.e., 9 ha
492440932 * 9 # 4 431 968 388 ha
487172235 * 9 # 4 384 550 115 ha
# 47418273 ha lost, 474.000 sq km

# merge random points with respective forest patch area
merged1992 <- cbind(res1992, for1992)
merged2020 <- cbind(res2020, for2020)

# log-transform patch size
merged1992$logarea <- log10(merged1992$value)
merged2020$logarea <- log10(merged2020$value)

# amount of habitat in "small" forests 
sum((merged1992$logarea) < 3.54)/length(merged1992$logarea)
sum((merged2020$logarea) < 3.54)/length(merged2020$logarea)

# calculate number of forest cells between 1992 and 2020
cellStats(lc1992, sum) # 492440932
cellStats(lc2020, sum) # 487172235

492440932 * 9 # 4 431 968 388 ha
487172235 * 9 # 4 384 550 115 ha
487172235/492440932
# we lost 1.1% of forests globally,  (492440932 - 487172235)*9
# 47418273 ha lost, 474.000 sq km

# check the proportion of patches smaller that different sizes
sum((merged1992$logarea) > 5.99)/length(merged1992$logarea)
sum((merged2020$logarea) > 5.99)/length(merged2020$logarea)

sum((merged1992$logarea) < 2)/length(merged1992$logarea) # 3% forest in patches smaller than 100 ha
sum((merged1992$logarea) < 3)/length(merged1992$logarea) # 8% forest in patches smaller than 1000 ha
sum((merged1992$logarea) < 4)/length(merged1992$logarea) # 12% forest in patches smaller than 10.000 ha
sum((merged1992$logarea) < 5)/length(merged1992$logarea) # 16% forest in patches smaller than 100.000 ha

sum((merged2020$logarea) < 2)/length(merged2020$logarea)
sum((merged2020$logarea) < 3)/length(merged2020$logarea)
sum((merged2020$logarea) < 4)/length(merged2020$logarea)
sum((merged2020$logarea) < 5)/length(merged2020$logarea)

# 3.5% of the world forests is in patches smaller than 1 km2 (100 ha)
# 7.5% smaller than 10 km2 (1000 ha)
# 11.5# of forests globally are smaller than 100 km2 (10^4 ha) 

sum((merged1992$logarea) > 5)/length(merged1992$logarea)
sum((merged2020$logarea) > 5)/length(merged2020$logarea)

##################################################
### 4 - MODELING
################################################## 

###patch size models
table(data$biome)/1000000 # proportion of points in different biomes; see metadata of original biome shapefile 
names(data)[names(data) == "value"] <- "PatchArea"

#data <- data[1:10000,] # this selects the first 10.000 rows, to speed up calculation

#patch area in ha: PatchArea
#prob_loss: binary, 1 means forest was lost, 0 means no forest loss
#log_area: log of PatchArea

## biome as factor
data$biome <- as.factor(data$biome)


model <- (glmmTMB(prob_loss ~ logarea + (logarea | biome), 
                  data = data, 
                  family = binomial, 
                  control = glmmTMBControl(parallel = 10)))

#sensitivity patches smaller than 1000 ha removed or more than 800000 ha; qualitatively equal results
model2 <- (glmmTMB(prob_loss ~ logarea + (logarea | biome), 
                   data = data[!(data$PatchArea>800000),], # change to 1000 ha to sensitivity to classification error
                   family = binomial, 
                   control = glmmTMBControl(parallel = 10)))

# report model coefficients
tab_model(model)
tab_model(model2)

# preliminary visualization
plot(allEffects(model))
plot(allEffects(model2))

# paper figure ggpredict
modelgg <- ggpredict(model,
                     c("logarea"),
                     type = "re"
)
plot(modelgg)

modelgg2 <- ggpredict(model2,
                      c("logarea"),
                      type = "re"
)
plot(modelgg2)

##################################################
### 5 - PLOT FIGURE 1
################################################## 

#plot: patch area by biome
plot_model <- ggpredict(model, 
                        c("logarea", "biome" ),
                        #c("log_area [all]", "biome" ),
                        type = "re") 

#reorder biomes
levels(plot_model$group)
levels(plot_model$group) <- c("Tropical & Subtropical Moist Broadleaf Forests",
                              "Tropical & Subtropical Dry Broadleaf Forests",
                              "Tropical & Subtropical Coniferous Forests",
                              "Temperate Broadleaf & Mixed Forests",
                              "Temperate Conifer Forests",
                              "Boreal Forests/Taiga",
                              "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                              "Temperate Grasslands, Savannas & Shrublands",
                              "Flooded Grasslands & Savannas",
                              "Montane Grasslands & Shrublands",
                              "Tundra",
                              "Mediterranean Forests, Woodlands & Scrub",
                              "Deserts & Xeric Shrublands",
                              "Mangroves")

plot_model$group <- factor(plot_model$group, levels = c("Tropical & Subtropical Moist Broadleaf Forests",
                                                        "Tropical & Subtropical Dry Broadleaf Forests",
                                                        "Tropical & Subtropical Coniferous Forests",
                                                        "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                        "Temperate Broadleaf & Mixed Forests",
                                                        "Temperate Conifer Forests",
                                                        "Temperate Grasslands, Savannas & Shrublands",
                                                        "Boreal Forests/Taiga",
                                                        "Montane Grasslands & Shrublands",
                                                        "Tundra",
                                                        "Mediterranean Forests, Woodlands & Scrub",
                                                        "Deserts & Xeric Shrublands",
                                                        "Flooded Grasslands & Savannas",
                                                        "Mangroves"))

#for labeller function to wrap text in facet labels
plot_model$group <- factor(plot_model$group, labels = c("Tropical & Subtropical Moist Broadleaf Forests",
                                                        "Tropical & Subtropical Dry Broadleaf Forests",
                                                        "Tropical & Subtropical Coniferous Forests",
                                                        "Tropical & Subtropical Grasslands, Savannas & Shrublands",
                                                        "Temperate Broadleaf & Mixed Forests",
                                                        "Temperate Conifer Forests",
                                                        "Temperate Grasslands, Savannas & Shrublands",
                                                        "Boreal Forests/Taiga",
                                                        "Montane Grasslands & Shrublands",
                                                        "Tundra",
                                                        "Mediterranean Forests, Woodlands & Scrub",
                                                        "Deserts & Xeric Shrublands",
                                                        "Flooded Grasslands & Savannas",
                                                        "Mangroves"))

#create data frame to add n= text to each facet #when I include the "group" argument, it changes the order of the data; without this argument (or when "cyl" it writes 14 labels on each facet on top of one another but in the right order
dat_text <- data.frame(
  label = c("1 n = ", "2 n = ", "3 n = ", "4 n = ", "5 n = ", "6 n = ", "7 n = ",
            "8 n = ", "9 n = ", "10 n = ", "11 n = ", "12 n = ", "13 n = ", "14 n = "),
  cyl   = c("Tropical & Subtropical Moist Broadleaf Forests",
            "Tropical & Subtropical Dry Broadleaf Forests",
            "Tropical & Subtropical Coniferous Forests",
            "Tropical & Subtropical Grasslands, Savannas & Shrublands",
            "Temperate Broadleaf & Mixed Forests",
            "Temperate Conifer Forests",
            "Temperate Grasslands, Savannas & Shrublands",
            "Boreal Forests/Taiga",
            "Montane Grasslands & Shrublands",
            "Tundra",
            "Mediterranean Forests, Woodlands & Scrub",
            "Deserts & Xeric Shrublands",
            "Flooded Grasslands & Savannas",
            "Mangroves"),
  x = 4,
  y = 0.8)

P2 <- plot(plot_model, 
           ci = TRUE,
           ci.style = "ribbon",
           limit.range = FALSE,
           collapse.group = FALSE,
           alpha = 0.08,
           line.size = 1.5,
           colors = c("black", # 30
                      "black", # 2
                      "black", # 1
                      "black", # 12
                      "black", # 6
                      "black", # 26
                      "black", # 14
                      "black", # 2
                      "black", # < 1
                      "black", # 1
                      "black", # 3
                      "black", # 1
                      "black", # 1
                      "black"))+# < 1
  labs(
    x = "Patch size (ha)", 
    y = "Probability of habitat loss", 
    title = ""  ) +
  facet_wrap(~group, 
             ncol = 7,
             labeller = label_wrap_gen(width = 25, multi_line = TRUE)) + #change width for num of characters for facet title wrapping (14-25)
  theme(strip.text.x = element_text(size = 14, family = "serif")) +
  theme_few() +
  scale_x_continuous(limits = c(2,6),
                     breaks = c(2, 3, 4, 5, 6), # 2, 4, 6
                     labels = c(expression(10^{2}),expression(10^{3}),expression(10^{4}),expression(10^{5}), expression(10^{6}))) +
  #100, 10,000, 1,000,000
  ylim(0,0.9) + #change so that full error bars are plotted
  theme(legend.position="none") +
  theme(axis.title=element_text(size = 16,  family="serif"),
        axis.text=element_text(size = 14, family = "serif"),
        strip.text = element_text(size = 12, family = "serif"),
        panel.spacing = unit(1.3, "lines"))#



P2

### plot map
raster_plot1 <- as.factor(raster_plot1)

#both dataframes have to have the same number of rows, but 1992 has more NAs than 2020 (only 19 NA rows)
#mask 2020 with 1992 to get rid of 19 extra cells
mask <- raster_plot1
raster_plot2 <- raster_plot2 * mask
raster_plot2 <- as.factor(raster_plot2)

r1_df <- as.data.frame(rasterToPoints(raster_plot1))
r2_df <- as.data.frame(rasterToPoints(raster_plot2))
r3_df <- r1_df

#create column with 1 1, 0 1 and 1 0, to plot as categories
r3_df$agg1992_30 <- as.factor(paste(r1_df$agg1992_30, r2_df$layer)) 

P1 <- ggplot() + 
  geom_tile(data = r3_df, aes(x = x, y = y, fill = agg1992_30)) +
  theme_few() +
  coord_fixed()+
  scale_fill_manual(labels = c("No forest", "Forest gained since 1992", "Forest lost since 1992", "Forest in 1992 and 2020"),
                    values=c("lavenderblush3", "yellow", "black", "forestgreen"))+ #red and black for science?
  theme(legend.position="bottom", # azure3, lavenderblush3, honeydew3, thistle3
        legend.title = element_blank(),
        axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.border=element_blank(),
        text=element_text(size= 24,  family="serif"))


P1
#put plots together (map and models)

figure <- ggarrange(P1, P2,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2,
                    heights = c(0.9,1),
                    widths = c(2,1),
                    font.label = list(size = 24, face = "bold", color ="black", family = "serif")#,
                    #align = "hv"
)

figure


#ggsave("figure.jpg", path = "Figures", width = 14000, height = 10000, units = "px",  device='jpeg', dpi=900, bg = "white")
#write.csv(plot_model, "predictions_forest_loss.csv") # table with predictions of habitat loss across patch sizes in different biomes


##################################################
### 6 - SMALL PATCHES IN ANTHROPOGENIC LANDSCAPES
################################################## 

# merge forest point coordinates and the patch area of the forest where each point fell
merged1992 <- cbind(res1992, for1992)
merged2020 <- cbind(res2020, for2020)

# calculate log10 of forest patch area
merged1992$logarea <- log10(merged1992$value)
merged2020$logarea <- log10(merged2020$value)


# aggregate to a 100 km cell to create "landscapes" and their degree of human activiy
pop <- raster::aggregate(pop, 333, fun = mean) 

# extract "pop" raster representing the proportion of anthropogenic land cover in a landscape
extract_pop <- raster::extract(pop, merged2020[,10:11])
extract_pop[is.na(extract_pop)] <- 0 # 0.002 of observations are NA, replace with 0

# convert into a dataframe
global_pop <- as.data.frame(rasterToPoints(pop))
# check the values of top 75%, 90% and 95% most anthropogenic landscapes
quantile(global_pop$lc2020_humandominated, probs = c(0.75, 0.9, 0.95))
# 0.504 is the 90th percentile, used in the paper

# convert into "huamn dominated" landscapes (> 0.504) or not
extract_pop_bin <- extract_pop
extract_pop_bin[extract_pop_bin <= 0.504] <- 0
extract_pop_bin[extract_pop_bin > 0.504] <- 1
extract_pop_bin <- as.data.frame(extract_pop_bin)
extract_pop_bin$logarea <- merged2020$logarea

densely_pop <- extract_pop_bin[extract_pop_bin$extract_pop_bin > 0,]
sum((densely_pop$logarea) < 3.5)/length(densely_pop$logarea)
sum((merged2020$logarea) < 3.5)/length(merged2020$logarea)

extract_pop_bin <- extract_pop
extract_pop_bin[extract_pop_bin <= 0.722] <- 0
extract_pop_bin[extract_pop_bin > 0.722] <- 1
extract_pop_bin <- as.data.frame(extract_pop_bin)
extract_pop_bin$logarea <- merged2020$logarea

densely_pop <- extract_pop_bin[extract_pop_bin$extract_pop_bin > 0,]
sum((densely_pop$logarea) < 3.5)/length(densely_pop$logarea)
sum((merged2020$logarea) < 3.5)/length(merged2020$logarea)

# ##################################################
# ### SUPPLEMENTARY CODE
# ################################################## 
# 
# ##################################################
# ### S1 RECLASSIFY ORIGINAL LAND COVER (ESA-CCI-LC)
# ################################################## 
# 
# #open land cover rasters
# lc1992 <- raster("Data/LC_1992.tif") # original ESA LULC from 1992, projected to equal Earth 
# lc2020 <- raster("Data/LC_2020.tif") # same as above, for 2020  
# 
# #reclassify 1992 into 1 for forest, 0 for non forest
# m <- c(0, 48, 0, # values between 0 and 48 become zeros
#        49, 101, 1, # values between 49 and 101 become ones
#        109, 158, 0,
#        159, 171, 1,
#        179, 208, 0,
#        209, 211, NA, #water as NA
#        212, 230, 0)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rc1 <- reclassify(lc1992, rclmat)
# 
# writeRaster(rc1, filename="lc1992_waterna", format="GTiff") #save raster
# 
# #reclassify lc2020 (same as for lc1992)
# m <- c(0, 48, 0,
#        49, 101, 1,
#        109, 158, 0,
#        159, 171, 1,
#        179, 208, 0,
#        209, 211, NA, #water as NA
#        212, 230, 0)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rc2 <- reclassify(lc2020, rclmat)
# 
# 
# writeRaster(rc2, filename="lc2020_waterna", format="GTiff") #save raster
# 
# 
# ##################################################
# ### S2 PREPARE BIOME RASTER
# ################################################## 
# 
# #aggregate projected raster to use as aggregation template
# agg30 <- aggregate(lc1992, 100) #aggregate to 30km
# writeRaster(agg30, filename="agg30", format="GTiff") #value between 0 and 1
# 
# #open biomes and reproject
# biomes <- readOGR("Data/wwf_terr_ecos.shp") # world biome shapefile; see supplementary for link to this original file
# mycrs <- "+proj=eqearth +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"
# biome_proj <- spTransform(biomes, crs(mycrs))
# 
# #rasterize biomes using fasterize
# biome_sf <- st_as_sf(biome_proj, crs = mycrs)
# 
# biome_rast <- fasterize(
#   biome_sf,
#   agg30,
#   field = "BIOME",
#   fun = "last"
# )
# 
# writeRaster(biome_rast, filename="biome_rast30", format="GTiff")
#
#
##################################################
### S3 GENERATE ANTHROPOGENIC LANDSCAPES
################################################## 
## reclassify land cover into anthropogenic vs natural land cover
#
# reclassify 1992 into 1 for human-dominated (urban and cropland), 0 for natural
# m <- c(0, 8, 0,
#        9, 28, 1, #cropland
#        29, 188, 0,
#        189, 198, 1, #urban
#        199, 208, 0,
#        209, 211, NA, #water as NA
#        212, 230, 0)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rc_human_1992 <- reclassify(lc1992, rclmat)
# rc_human_1992
# 
# writeRaster(rc_human_1992, filename="lc1992_humandominated", format="GTiff") #save raster
# 
# #reclassify lc2020 (same as for lc1992)
# m <- c(0, 8, 0,
#        9, 28, 1, #cropland
#        29, 188, 0,
#        189, 198, 1, #urban
#        199, 208, 0,
#        209, 211, NA, #water as NA
#        212, 230, 0)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# rc_human_2020 <- reclassify(lc2020, rclmat)
# rc_human_2020
# writeRaster(rc2, filename="lc2020_reclass", format="GTiff")
# writeRaster(rc_human_2020, filename="lc2020_humandominated", format="GTiff") #save raster
# 
# ##################################################
# ### S4 SAMPLE FOREST PLOTS
# ################################################## 
# 
# #create random points across globe for land cover 1992 or 2020 raster
# x <- sampleRandom(lc2020, 4000000, na.rm=TRUE, xy=TRUE) #e.g., 30% of terrestrial systems are forests, 29% of planet is land
# x <- as.data.frame(x)
# #only keep forest points
# x <- x[!(x$lc2020_waterna==0),] 
# 
# write.csv(x, "data_forest2020.csv")
# 
# #function to get patch size (50km in all directions from point, end up with a maximum patch area of 100km x 100km)
# PATCH_SIZE <- function(number, list_points, raster_cover){
#   prova <- list_points[number,1:2]
#   prova_raster <- crop(raster_cover,
#                        extent((prova$x - 50000), (prova$x + 50000), (prova$y - 50000), (prova$y + 50000)))
#   object <- extract_lsm(prova_raster, y = as.matrix(prova), what = "lsm_p_area", directions = 4)
#   return(object)
# }
# 
# 
# # parallelize
# numCores <- detectCores(logical = TRUE) - 1
# numCores
# ## [1] we have 7 cores
# cl <- makeCluster(numCores)
# 
# start_time <- Sys.time()
# 
# parallel::clusterEvalQ(cl, {
#   library(raster)
#   library(rgdal)
#   library(dplyr)
#   library(sp)
#   library(landscapemetrics)
#   
# })
# 
# parallel::clusterExport(cl, varlist = c("x", "lc2020"))
# 
# maps <- parallel::parLapply(cl, # the list of parallel R running
#                             seq_len(nrow(x)),
#                             # replicate the function on the name of the species, from table_species; the parentheses limit this to 20 species
#                             PATCH_SIZE, # the function to replicate
#                             list_points = x,
#                             raster_cover = lc2020)
# 
# stopCluster(cl)
# 
# end_time <- Sys.time()
# end_time - start_time
# 
# # save as a csv
# results <- do.call(rbind, maps)
# write.csv(results,"results2020.csv")
# 
# 
# ##################################################
# ### S5 DATA PREPARATION FOR MODELING
# ################################################## 
# 
# # random sample of forest cells 
# for1992 <- read.csv("data_forest1992.csv")
# for2020 <- read.csv("data_forest2020.csv")
# 
# # area of the forest in which each randomly sampled point fell
# res1992 <- read.csv("results1992.csv")
# res2020 <- read.csv("results2020.csv")
# 
# # merge random points with respective forest patch area
# merged1992 <- cbind(res1992, for1992)
# merged2020 <- cbind(res2020, for2020)
# 
# #remove patches smaller than 100 ha to minimize classification error
# merged1992 <- merged1992[merged1992$logarea > 2, ]                          
# merged2020 <- merged2020[merged2020$logarea > 2, ]   
# 
# #extract if point falls within forest in 2020
# extract_2020 <- raster::extract(lc2020, merged1992[,10:11])
# merged1992$extract_2020 <- extract_2020 #create new column with forest values for 2020
# 
# #create new column to calculate probability that forest patch was lost in 2019
# merged1992$prob_loss <- abs(merged1992$extract_2020 - 1)
# 
# #extract biome for sampled points
# extract_biome <- raster::extract(biome_rast, merged1992[,10:11])
# merged1992$biome <- extract_biome #add biome column
# 
# #drop rows with biome = 98 and 99, not terrestrial
# results_final <- merged1992[!(merged1992$biome==98),]
# results_final <- results_final[!(results_final$biome==99),]
# results_final <- na.omit(results_final) #make sure there are no NA rows
# 
# # keep a million observations
# results_final <- results_final[1:1000000,]
# 
# #save final csv for analysis
# #write.csv(results_final,"results_final.csv")
# 
# names(results_final)[names(results_final) == "value"] <- "PatchArea"