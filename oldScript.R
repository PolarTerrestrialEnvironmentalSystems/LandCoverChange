#################################
### Template: LandCoverChange ###
### 22.10.2023                ###
### author: Simeon Lisovski   ###
###         Peter Ewald       ###
###         Thomas Böhmer     ###
#################################

source("/bioing/data/LandCoverChange/Template/LandCoverChange_mod.R")

## Packages
library(tensorflow)
library(tidyr)
library(dplyr)
library(tibble)
library(readxl)
library(fs)
library(rgee)
library(geojson)
library(sf); sf_use_s2(FALSE)
library(stars)
library(vegan)
library(zoo)
library(patchwork)
library(scales)
library(transport)
library(corit)
library(ggplot2)

#############
### Setup ###
#############

check_files <- TRUE

#### Google Earth Engine ####
{
  ## User: slisovsk
  {
    # rgee::ee_install_set_pyenv(
    #  py_path = "/home/slisovsk/.local/share/r-miniconda/envs/ee/bin/python"
    # )
    # ee_Initialize(drive = TRUE, user = "simeon.lisovski@gmail.com")
  }
  
  ## User: tboehmer
  {
    ee_Initialize(drive = TRUE, user = "th.boehmer77@gmail.com") ## this function needs to run without errors
  }
  
}

#### selection of grid cells ####
{
  
  grid_folder1 = "//bioing/data/LandCoverChange/"
  
  # region_grid_global <- read.csv(paste0(grid_folder1,"gridcells_with_pollen_centers_50km_sort.csv"), header=TRUE)
  # 
  # # region_grid_nam <- subset(region_grid_global, Long < -50 & Lat >= 13)
  # region_grid_asia <- subset(region_grid_global, Long > 43 & Lat >= 18)
  
  
  
  PSA_sites_global <- read.csv(paste0(grid_folder1,"PSA_locations_northern_hemisphere.csv"), header=TRUE)
  
  PSA_sites_global <- PSA_sites_global %>% filter(Continent == "Europe") %>% filter(Latitude >= 52 & Latitude <= 55) ## remove Lat-filter to run all sites
  
}


runindex <- sort(PSA_sites_global$Dataset_ID)


for(rr in 1:length(runindex)){
  
  cat("\n")
  print(paste0("region grid: ",runindex[rr]," (",rr,"/",length(runindex),")"))
  
  tryCatch(
    expr={
      
      # current_gridcell <- region_grid_global[region_grid_global$grid_ID_new == runindex[rr], ]
      
      # current_long <- ifelse(current_gridcell$Long > 0, paste0(round(current_gridcell$Long, 2),"E"), paste0(round(abs(current_gridcell$Long), 2),"W")) 
      # current_lat  <- ifelse(current_gridcell$Lat > 0, paste0(round(current_gridcell$Lat, 2),"N"), paste0(round(abs(current_gridcell$Lat), 2),"S")) 
      
      current_psa  <- PSA_sites_global[PSA_sites_global$Dataset_ID == runindex[rr], ]
      
      current_long <- ifelse(current_psa$Long > 0, paste0(round(current_psa$Long, 2),"E"), paste0(round(abs(current_psa$Long), 2),"W")) 
      current_lat  <- ifelse(current_psa$Lat > 0, paste0(round(current_psa$Lat, 2),"N"), paste0(round(abs(current_psa$Lat), 2),"S")) 
      
      source("/bioing/data/LandCoverChange/Template/LandCoverChange_mod.R")
      
      #### Directory Setup ####
      {
        data_dir  <- "/bioing/data/LandCoverChange/LandCoverChangeProject_data_new/"
        dataFiles <- tibble(path = list.files(data_dir, recursive = T)) %>%
          mutate(folder = sapply(strsplit(gsub("data/", "", path), "/"), function(x)  x[[1]]),
                 file = sapply(strsplit(gsub("data/", "", path), "/"), function(x)  x[[2]])) %>%
          dplyr::select(folder, file)
        # knitr::kable(dataFiles)
        
        # dir_out <- paste0("/bioing/data/LandCoverChange/TestRuns_gridded_regions_Pollen_Asia_newTransport/gridcell_",current_gridcell$grid_ID_new,"__",current_lat,"_",current_long,"/")
        dir_out <- paste0("/bioing/data/LandCoverChange/TestRuns_PSA_maps_Reveals_Europe_HumanInfluence_v2/PSA_",current_psa$Dataset_ID,"__",current_lat,"_",current_long,"/")
        if(!file.exists(dir_out)) {dir.create(dir_out)
          createFolders(dir_out)
          fs::dir_tree(path = dir_out, recurse = TRUE)
        }
      }
      
      #### Spatial extent and resolution ####
      {
        # define map and projection, extent and resolution
        # centre      <- c(13.43, 55.93)
        # centre      <- c(round(current_gridcell$Long,2), round(current_gridcell$Lat,2))
        # extent      <- 50               # [km]
        
        centre      <- c(round(current_psa$Long,2), round(current_psa$Lat,2))
        extent      <- ifelse(current_psa$Pollen_Source_Radius > 250000, 250000, current_psa$Pollen_Source_Radius) / 1000  # [km]
        
        buffer      <- 1000             # buffer to include PSAs for calibration [km]
        
        ## Native resolution = 309.2208 meter
        resolution   <- 300            # [m] possible options: 300, 3000, 5000, 10000 
        
        # Calibration Map (define separate calibration groups)
        calibrationMap = suppressWarnings(
          rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
            dplyr::select("continent") %>% rename(groups = continent)
        )
      }
      
      
      #### Create vegetation object for PSAs ####
      {
        # load data, filter pollen source areas in region and/or calibration subset
        if(!file.exists(glue::glue("{dir_out}/summaryResults/psaVeg.rda")) & check_files) {
          psaVeg <- make_psaVeg(centre, extent, buffer, resolution, proj = "laea",
                                calibrationMap = calibrationMap,
                                veg_folder = file.path(data_dir, "data", "vegetation_cover"), 
                                interpolate = TRUE, 
                                fc = 1/500, dt = 500, k = 5, cluster = NULL)
          save(psaVeg, file = glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
          
          plot(psaVeg)
          ggsave(glue::glue("{dir_out}/summaryResults/Map.png"), width = 10, height = 10, units = "cm")
          
        } else load(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
      }
      
      #### Plot vegetation time series for specific sites ####
      ## options: plot only n_taxa taxa with largest mean abundance
      # plotVegSite(psaVeg, ID = unique(psaVeg@vegetation$Dataset_ID)[10], n_taxa = 10, nrow = 2, ncol = 5)
      ## plot
      
      #################
      ### Landcover ###
      #################
      
      #### Modern LandCover share ####
      {
        modern_age_threshold <- 1000 ### modern time threshold for calibration
        
        modernVegetation <- psaVeg@vegetation %>% filter(Interp, Age <= modern_age_threshold) %>%
          dplyr::select(-c(Interp, Age)) %>%
          group_by(Dataset_ID) %>%
          summarise(across(everything(), function(x) mean(x, na.rm = TRUE))) %>% filter(!apply(., 1, function(x) all(is.nan(x[-1]))))
        
        class_defs <- readxl::read_excel(glue::glue("{data_dir}settings/landcover_class_definitions_mod.xlsx")) 
      }
      
      ### CCI LandCover map of study area ####
      {
        roi = sf_as_ee(psaVeg@regionMap %>% st_transform(4326) %>% st_shift_longitude())
        
        dataset = ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(roi)
        
        if(!file_exists(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff")) & check_files) {
          dataset_rast <- ee_as_stars(dataset, via = 'drive', quiet = TRUE)
          write_stars(dataset_rast, glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
        }
        
        # extract number of pixel per class
        class_areas  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
          collection = roi,
          reducer    = ee$Reducer$count()$group(groupField = 1),
          scale      = dataset$projection()$nominalScale()$getInfo()
        )$getInfo()
        
        regLCov = lapply(class_areas$features[[1]][[4]][[1]], function(x) {
          tibble(Class_Code = x$group, count = x$count)
        }) %>% Reduce("rbind",.) %>% left_join(class_defs %>% dplyr::select(Class_Code, Variable_Name, Class_Label), by = join_by(Class_Code)) %>%
          mutate(Percent = round((count/sum(count))*100,2)) %>% relocate(Percent, .after = Variable_Name) %>%
          relocate(count, .after = Variable_Name) %>% rename(Count = count) %>% arrange(desc(Count))
        
        # knitr::kable(regLCov)
      }
      
      ### Modern LandCover shares in PSAs ####
      {
        ## remove classes for calibration that have very little pixels in study area
        class_defs[class_defs$Include_Class_In_Calibration & is.na(class_defs$Merge_To_Class) & 
                     class_defs$Class_Code %in% regLCov$Class_Code[regLCov$Percent<0.25], 'Include_Class_In_Calibration'] <- FALSE
        
        lcovShares <- modernLC(psaVeg, psaIDs = unique(modernVegetation$Dataset_ID), resolution = resolution, class_defs)
        
        
        ## plot
        {
          # example plot of 3 ids
          # pieTab <- lcovShares %>% filter(Dataset_ID%in%unique(Dataset_ID)[sample(1:nrow(.), 3)]) %>%
          #   pivot_longer(cols = -Dataset_ID, values_to = "Percentage") %>% mutate(Class_Code = as.numeric(name)) %>%
          #   left_join(class_defs %>% dplyr::select(Class_Code, Class_Label, Color_Code), by = join_by(Class_Code))
          # 
          # ggplot(pieTab, aes(x = '', y = Percentage, fill = as.factor(Class_Code))) +
          #   geom_bar(stat="identity", width=1) +
          #   scale_fill_manual(values = pieTab$Color_Code, breaks = pieTab$Class_Code, name = "Landcover class", labels = pieTab$Class_Label) +
          #   coord_polar("y", start=0) +
          #   facet_wrap(~Dataset_ID, ncol = 3) +
          #   theme_light()
        }
      }
      
      ### Translate Vegetation to Landcover ####
      {
        if(!file.exists(glue::glue("{dir_out}/summaryResults/veg_lc_model.rda")) & check_files) {
          
          veg_lc_model <- regressionModel(psaVeg, lcovShares, vegetation_modern_years = 1000, prefilter_taxa = T) %>% 
            suppressWarnings()
          
          save(veg_lc_model, file = glue::glue("{dir_out}/summaryResults/veg_lc_model.rda"))
        } else load(glue::glue("{dir_out}/summaryResults/veg_lc_model.rda"))
        
        ## summary
        {
          # lapply(1:length(veg_lc_model), function(x) tibble(Region = names(veg_lc_model)[x]) %>% bind_cols(
          #   matrix(veg_lc_model[[x]]$sd_per_class, nrow = 1) %>%
          #     as_tibble() %>% setNames(names(veg_lc_model[[x]]$sd_per_class)))) %>%
          #   Reduce("rbind",.) %>%
          #   knitr::kable(caption = "Standard deviation (sd) per class and region") %>% suppressWarnings()
        }
        
        if(!file.exists(glue::glue("{dir_out}/summaryResults/psaLCover.rda")) & check_files) {
          psaLCover <- make_psaLCover(psaVeg, lcovShares, veg_lc_model)
          save(psaLCover, file = glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
        } else load(glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
        
        # knitr::kable(psaLCover@landcover_ts[1:10,])
        
      }
      
      
      ### stratigraphic plots ###
      
      {
        
        if(file.exists(glue::glue("{dir_out}/summaryResults/psaLCover.rda"))){
          
          lc_ts_df <- as.data.frame(psaLCover@landcover_ts)
          colnames(lc_ts_df)[-c(1:2)] <- paste0("LC_",names(lc_ts_df[,-c(1:2)]))
          
          # lc_ts_df1 <- lc_ts_df %>% 
          #   filter(Age > 1000)
          # 
          # psa_lc_index <- unique(lc_ts_df1$Dataset_ID)
          psa_lc_index <- current_psa$Dataset_ID
          
          if(length(psa_lc_index) > 0){
            
            for(j in 1:length(psa_lc_index)){
              
              # print(paste0(j,"/",length(psa_lc_index)))
              
              lc_ts_df_sub <- lc_ts_df %>% filter(Dataset_ID == psa_lc_index[j])
              
              lc_ts_df_plot  <- lc_ts_df_sub %>% 
                pivot_longer(cols = -c(Dataset_ID, Age), values_to = "Percentage")
              
              lc_ts_df_plot$name_sort = factor(lc_ts_df_plot$name, levels=stringr::str_sort(unique(lc_ts_df_plot$name), numeric = TRUE)) 
              
              
              p1 <- ggplot(lc_ts_df_plot %>% filter(!is.na(Percentage))) +
                geom_path(aes(x = Percentage, y = Age, group = as.factor(name)), color = "#E69F00") +  
                facet_wrap(~name_sort, nrow = 4, ncol = 5) +
                scale_y_reverse(breaks=seq(0,max(lc_ts_df_plot$Age),2000)) +
                labs(title = glue::glue("Pollen source area ID: {psa_lc_index[j]}")) +
                theme_light() +
                theme(legend.position = "top")
              
              
              ggsave(p1, file = glue::glue("{dir_out}/summaryResults/stratplot_LC_PSA_",psa_lc_index[j],".png"), width = 20, height = 20, units = "cm", dpi = 300)
              
            } # end for
            
          } # end if
          
        } # end if
        
        # -------------------------------------------------------------------------------------------------
        
        if(file.exists(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))){
          
          if(length(psa_lc_index) > 0){
            
            for(j in 1:length(psa_lc_index)){
              
              # print(paste0(j,"/",length(psa_lc_index)))
              
              n_taxa = 10
              meta_cols = c("Dataset_ID", "Age", "Interp")
              
              df_plot = psaVeg@vegetation %>% filter(Dataset_ID == psa_lc_index[j])
              
              # select n_taxa taxa
              subtitle <- ifelse(is.null(n_taxa), "All taxa available", glue::glue("Selection of {n_taxa}")) 
              n_taxa   <- ifelse(is.null(n_taxa), ncol(df_plot)-length(meta_cols), n_taxa)
              df_plot  <- df_plot %>% dplyr::select(c(meta_cols, 
                                                      names(sort(apply((df_plot %>% dplyr::select(-meta_cols)), 2, mean, na.rm = T), decreasing = T))[1:n_taxa])) %>%
                pivot_longer(cols = -c(meta_cols), values_to = "Percentage")
              
              
              p2 <- ggplot(df_plot %>% filter(!is.na(Percentage))) +
                geom_path(aes(x = Percentage, y = Age, group = as.factor(Interp), color = as.factor(Interp))) +
                scale_colour_manual(values = c("#E69F00", "grey60"), name = "Interpolated") +
                facet_wrap(~name, nrow = 2, ncol = 5) +
                scale_y_reverse(breaks=seq(0,max(lc_ts_df_plot$Age),2000)) +
                labs(title = glue::glue("Pollen source area ID: {psa_lc_index[j]}"),
                     subtitle = subtitle) +
                theme_light() +
                theme(legend.position = "top")
              
              
              ggsave(p2, file = glue::glue("{dir_out}/summaryResults/stratplot_major_taxa_PSA_",psa_lc_index[j],".png"), width = 20, height = 20, units = "cm", dpi = 300)
              
            } # end for
            
          } # end if
          
        } # end if
        
        # -------------------------------------------------------------------------------------------------
        
      }
      
      
      
      
      ############
      ### RDAs ###
      ############
      
      # psaLcovTS <- as.data.frame(psaLCover@landcover_ts)
      # psaLcovTS <- psaLcovTS %>% 
      #   filter(Age > 1000)
      
      
      # IDs             <- unique(psaLcovTS$Dataset_ID)   #unique(psaLCover@psas$Dataset_ID)
      IDs             <- current_psa$Dataset_ID 
      
      pxl_extract_PSA     <- 500       # max. points for RDA to search within the PSA   
      pxl_extract_PSA_min <- 50        # min. points for RDA to search in buffer rings if they are not found in the PSA 
      
      # pxl_extract_max <- 100  # max. points for RDA                              500
      max_buffer      <- 500  # buffer to include PSAs for calibration [km]      1000
      method          <- "getInfo" ## Options "ee_as_sf", "getInfo"
      
      ### Landcover classes: coordinate selection (TIME CONSUMING) ####
      ### Make sure PC/Mac does not go into sleep modus.
      {
        for(ID in IDs) {
          
          # cat(glue::glue("\rPSA {which(unique(psaLCover@psas$Dataset_ID)==ID)} of {length(unique(psaLCover@psas$Dataset_ID))}"))
          # cat(glue::glue("\rPSA {which(unique(psaLcovTS$Dataset_ID)==ID)} of {length(unique(psaLcovTS$Dataset_ID))}"))
          
          
          # buffer  <- lapply(bufferS, function(x) psa %>% st_buffer(x))
          # 
          # ggplot() +
          #   geom_stars(data = initRast[1], mapping = aes(fill = as.factor(Landcover))) +
          #   geom_sf(psaVeg@psas %>% filter(Dataset_ID %in% IDs) %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="cyan") +
          # 
          #   geom_sf(buffer[[1]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black", linetype="dotted") +
          #   geom_sf(buffer[[2]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[3]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[4]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[5]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[6]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[7]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[8]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[9]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          #   geom_sf(buffer[[10]] %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="black") +
          # 
          #   scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
          #                     breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
          #                     labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
          #   xlab("") + ylab("") +
          #   theme_light() +
          #   theme(legend.position = "none")
          
          
          if(!file_exists(glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda")) & check_files) {
            
            psa      <- psaVeg@psas %>% filter(Dataset_ID==ID)
            bufferS  <- seq(0, max_buffer*1000, length = 10)
            
            nPxl    <- pixelsBuffer(psa, psaLCover@resolution, class_defs, bufferS, cutoff = pxl_extract_PSA) %>%
              
              full_join(tibble(lcov   = as.factor(rep(as.numeric(names(psaLCover@landcover_ts)[-c(1,2)]), length(bufferS))), 
                               radius = rep(bufferS, each = length(as.numeric(names(psaLCover@landcover_ts)[-c(1,2)]))), count = 0), 
                        by = join_by(lcov, radius, count)) %>%
              arrange(lcov, radius, count) %>%
              group_split(lcov) %>% lapply(., function(x) if(any(x$count>=pxl_extract_PSA)) x[min(which(x$count>pxl_extract_PSA)),] else x[which.max(x$count),]) %>% 
              Reduce("rbind",.) %>% mutate(count = ifelse(count>pxl_extract_PSA, pxl_extract_PSA, count)) %>% filter(!duplicated(lcov))
            
            
            #################################
            
            # original approach:
            {
              # psaCrds <- lapply(nPxl$lcov, function(class) {
              #   if((nPxl %>% filter(lcov == class) %>% pull(count)) > 0) {
              # 
              #     bbox_ee   <- psa %>% st_buffer(nPxl$radius[nPxl$lcov==class]) %>% sf_as_ee()
              #     dataset   <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(bbox_ee)
              #     # Map$addLayer(dataset) +
              #     #   Map$addLayer(bbox_ee)
              # 
              #     # A list of pixel values to replace.
              #     fromList = class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
              #       filter(Predicted_In_Past, Merge_To_Class == class) %>% pull(Class_Code) %>% as.list()
              # 
              #     toList = rep(1, length(fromList)) %>% as.list()
              # 
              #     classImage = dataset$remap(
              #       from = fromList,
              #       to = toList,
              #       defaultValue = -1,
              #       bandName = 'b1')
              # 
              #     dts <- dataset$updateMask(classImage$gt(0))
              # 
              #     if(method == "ee_as_sf") {
              # 
              #       if((nPxl %>% filter(lcov == class) %>% pull(count)) < 1000) {
              #         vectors = dts$sample(
              #           region = bbox_ee,
              #           geometries = TRUE)
              #       } else {
              #         vectors = dts$sample(
              #           region = bbox_ee,
              #           geometries = TRUE,
              #           scale = (dts$projection()$nominalScale()$getInfo()),
              #           numPixels = (nPxl %>% filter(lcov == class) %>% pull(count)))
              #       }
              # 
              #       suppressMessages({
              #         t <- ee_as_sf(vectors$limit(nPxl %>% filter(lcov == class) %>% pull(count)),
              #                       via = 'drive', container = "rgee_backup", quiet = T, add_metadata = FALSE)
              #       })
              # 
              #       t %>% mutate(Dataset_ID = ID,
              #                    lcov = as.numeric(as.character(class)),
              #                    dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
              #         dplyr::select(Dataset_ID, lcov, dist)
              # 
              #     } else {
              # 
              #       Map$addLayer(dts) +
              #         Map$addLayer(bbox_ee)
              # 
              #       table = dts$sample(
              #         region = bbox_ee,
              #         scale = (dts$projection()$nominalScale()$getInfo()),
              #         geometries = TRUE
              #       )$limit(nPxl %>% filter(lcov == class) %>% pull(count))$getInfo()
              # 
              #       lapply(table[[4]], function(f) {
              #         tibble(Dataset_ID = ID,
              #                lcov       = as.numeric(as.character(class)),
              #                dist       = NA,
              #                lon        = f$geometry$coordinates[1],
              #                lat        = f$geometry$coordinates[2])
              #       }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
              #         mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
              #         suppressWarnings()
              # 
              #     }
              # 
              #   } else NULL
              # 
              # }) %>% Reduce("rbind",.)
              }
            
            # approach with rings:
            {
            # psaCrds <- lapply(nPxl$lcov, function(class) {
            #   if((nPxl %>% filter(lcov == class) %>% pull(count)) > 0) {
            #     
            #     PxlToSample <- nPxl %>% filter(lcov == class) %>% pull(count)
            #     sampledPxl  <- 0
            #     
            #     Pxl_LCClass <- NULL
            #     
            #     repetition <- 0
            #     
            #     while(sampledPxl < nPxl[nPxl$lcov==class,"count"]){
            #       
            #       if(repetition == 0){
            #         
            #         bbox_ee   <- psa %>% st_buffer(bufferS[repetition + 1]) %>% sf_as_ee()
            #         
            #         cat(paste("\rLC-class:",nPxl %>% filter(lcov == class) %>% pull(lcov)," | repetition:",repetition+1))
            #         # cat("\n")
            #         
            #       }else{
            #         
            #         bbox_ee_inner   <- psa %>% st_buffer(bufferS[repetition + 1]) %>%  st_as_sfc() # %>% sf_as_ee()
            #         bbox_ee_outer   <- psa %>% st_buffer(bufferS[repetition + 2]) %>%  st_as_sfc() # %>% sf_as_ee()
            #         
            #         # remove bbox_ee_inner from bbox_ee_outer
            #         bbox_ee <- st_difference(bbox_ee_outer, bbox_ee_inner) %>% sf_as_ee()
            #         
            #         cat(paste("\rLC-class:",nPxl %>% filter(lcov == class) %>% pull(lcov)," | repetition:",repetition+2))
            #         # cat("\n")
            #         
            #       }
            #       
            #       # bbox_ee   <- psa %>% st_buffer(bufferS[repetition + 1]) %>% sf_as_ee()
            #       dataset   <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(bbox_ee)
            #       
            #       # A list of pixel values to replace.
            #       fromList = class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
            #         filter(#Predicted_In_Past, 
            #           Merge_To_Class == class) %>% pull(Class_Code) %>% as.list()
            #       
            #       toList = rep(1, length(fromList)) %>% as.list()
            #       
            #       classImage = dataset$remap(
            #         from = fromList,
            #         to = toList,
            #         defaultValue = -1,
            #         bandName = 'b1')
            #       
            #       dts <- dataset$updateMask(classImage$gt(0))
            #       
            #       
            #       Map$addLayer(dts) +
            #         Map$addLayer(bbox_ee)
            #       
            #       table = dts$sample(
            #         region = bbox_ee,
            #         scale = (dts$projection()$nominalScale()$getInfo()),
            #         geometries = TRUE
            #       )$limit(PxlToSample)$getInfo()
            #       
            #       
            #       if(length(table[[4]]) > 0){
            #         Pxlfound <-  lapply(table[[4]], function(f) {
            #         tibble(Dataset_ID = ID, 
            #                lcov       = as.numeric(as.character(class)), 
            #                dist       = NA,
            #                lon        = f$geometry$coordinates[1],
            #                lat        = f$geometry$coordinates[2])
            #       }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
            #         mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
            #         suppressWarnings()
            #       
            #     
            #       Pxl_LCClass <- rbind(Pxl_LCClass, Pxlfound)
            #     
            #     
            #       PxlToSample <- (nPxl %>% filter(lcov == class) %>% pull(count)) - dim(Pxl_LCClass)[1]
            #       sampledPxl  <- dim(Pxl_LCClass)[1]
            #       
            #       } # end if (table > 0)
            #       
            #       repetition <- repetition + 1
            #       
            #   
            #       
            #     } # end while
            #     
            #     cat("\n")
            #    
            #   } else NULL
            #    
            #   return(Pxl_LCClass)
            #   
            # }) %>% Reduce("rbind",.)
           }
               
            #################################
            
            
            psaCrds <- lapply(nPxl$lcov, function(class) {
              if((nPxl %>% filter(lcov == class) %>% pull(count)) > 0) {
                
                cat(paste("\rLC-class:",nPxl %>% filter(lcov == class) %>% pull(lcov)))
                
                Pxl_LCClass <- NULL
                
                # area of the PSA
                bbox_ee   <- psa %>% sf_as_ee()
                dataset   <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(bbox_ee)
                
                # Map$addLayer(dataset) +
                #   Map$addLayer(bbox_ee)

                # A list of pixel values to replace.
                fromList = class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
                  filter(Predicted_In_Past, Merge_To_Class == class) %>% pull(Class_Code) %>% as.list()

                toList = rep(1, length(fromList)) %>% as.list()

                classImage = dataset$remap(
                  from = fromList,
                  to = toList,
                  defaultValue = -1,
                  bandName = 'b1')

                dts <- dataset$updateMask(classImage$gt(0))

                if(method == "ee_as_sf") {

                  if((nPxl %>% filter(lcov == class) %>% pull(count)) < 1000) {
                    vectors = dts$sample(
                      region = bbox_ee,
                      geometries = TRUE)
                  } else {
                    vectors = dts$sample(
                      region = bbox_ee,
                      geometries = TRUE,
                      scale = (dts$projection()$nominalScale()$getInfo()),
                      numPixels = (nPxl %>% filter(lcov == class) %>% pull(count)))
                  }

                  suppressMessages({
                    t <- ee_as_sf(vectors$limit(nPxl %>% filter(lcov == class) %>% pull(count)),
                                  via = 'drive', container = "rgee_backup", quiet = T, add_metadata = FALSE)
                  })

                  t %>% mutate(Dataset_ID = ID,
                               lcov = as.numeric(as.character(class)),
                               dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
                    dplyr::select(Dataset_ID, lcov, dist)

                } else {

                  # Map$addLayer(dts) +
                  #   Map$addLayer(bbox_ee)

                  table = dts$sample(
                    region = bbox_ee,
                    scale = (dts$projection()$nominalScale()$getInfo()),
                    geometries = TRUE
                  )$limit(pxl_extract_PSA)$getInfo()

                  
                  if(length(table[[4]]) > 0){
                  Pxlfound <-  lapply(table[[4]], function(f) {
                    tibble(Dataset_ID = ID, 
                           lcov       = as.numeric(as.character(class)), 
                           dist       = NA,
                           lon        = f$geometry$coordinates[1],
                           lat        = f$geometry$coordinates[2])
                  }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
                    mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
                    suppressWarnings()

                 Pxl_LCClass <- rbind(Pxl_LCClass, Pxlfound)
                  } # end if (table > 0)
                  
                } # end else

                
                Pxl_LCClass_found <- ifelse(length(dim(Pxl_LCClass)[1]) == 0, 0, dim(Pxl_LCClass)[1])
                
                
                if(Pxl_LCClass_found < pxl_extract_PSA_min){
                  
                      PxlToSample <- pxl_extract_PSA_min - Pxl_LCClass_found
                      sampledPxl  <- 0

                      repetition <- 1
                   
                      while(sampledPxl < pxl_extract_PSA_min){
                        
                        if(bufferS[repetition] == tail(bufferS, 1)){
                          
                         
                          break
                          
                          }
                        
                        bbox_ee_inner   <- psa %>% st_buffer(bufferS[repetition]) %>%  st_as_sfc() # %>% sf_as_ee()
                        bbox_ee_outer   <- psa %>% st_buffer(bufferS[repetition + 1]) %>%  st_as_sfc() # %>% sf_as_ee()
                        
                        # remove bbox_ee_inner from bbox_ee_outer
                        bbox_ee <- st_difference(bbox_ee_outer, bbox_ee_inner) %>% sf_as_ee()
                        
                        cat(paste("\rLC-class:",nPxl %>% filter(lcov == class) %>% pull(lcov)," | ring:",repetition+1))
                        # cat("\n")
                  
                  
                        dataset   <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(bbox_ee)
                        
                        # A list of pixel values to replace.
                        fromList = class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
                          filter(#Predicted_In_Past,
                            Merge_To_Class == class) %>% pull(Class_Code) %>% as.list()
                        
                        toList = rep(1, length(fromList)) %>% as.list()
                        
                        classImage = dataset$remap(
                          from = fromList,
                          to = toList,
                          defaultValue = -1,
                          bandName = 'b1')
                        
                        dts <- dataset$updateMask(classImage$gt(0))
                        
                        
                        # Map$addLayer(dts) +
                        #   Map$addLayer(bbox_ee)
                        
                        table = dts$sample(
                          region = bbox_ee,
                          scale = (dts$projection()$nominalScale()$getInfo()),
                          geometries = TRUE
                        )$limit(PxlToSample)$getInfo()  
                        
                        
                        if(length(table[[4]]) > 0){
                          Pxlfound_Buf <-  lapply(table[[4]], function(f) {
                            tibble(Dataset_ID = ID,
                                   lcov       = as.numeric(as.character(class)),
                                   dist       = NA,
                                   lon        = f$geometry$coordinates[1],
                                   lat        = f$geometry$coordinates[2])
                          }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
                            mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
                            suppressWarnings()
                          
                    
                          Pxl_LCClass <- rbind(Pxl_LCClass, Pxlfound_Buf)
                          
                          PxlToSample <- ifelse(pxl_extract_PSA_min - dim(Pxl_LCClass)[1] < 0, 0, pxl_extract_PSA_min - dim(Pxl_LCClass)[1])   ###  (nPxl %>% filter(lcov == class) %>% pull(count)) - dim(Pxl_LCClass)[1]
                          sampledPxl  <- dim(Pxl_LCClass)[1]
                          
                        } # end if (table > 0)
                        
                        repetition <- repetition + 1
                          
                      } # end while
                        
                } # end if (Pixel < pxl_extract_PSA)
                
                cat("\n")
                
              } else NULL
              
              return(Pxl_LCClass)
              
              cat("\n")
              
            }) %>% Reduce("rbind",.)  
                      
                 
            psaCrds_summary <- psaCrds %>% st_drop_geometry() %>% group_by(lcov) %>% summarise(pixel_found = n())
            
            # save pixelflow table:
            openxlsx::write.xlsx(psaCrds_summary, file = glue::glue("{dir_out}/summaryResults/psaCrds_summary_{ID}.xlsx"), asTable = TRUE)
            
            
            ### save results
            save(psaCrds, file = glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda"))
            
          }
          
          cat("\b\b\b\b\b\b\b\b\b")
          
        }
      }
      
      
      ### Add Environmental variables
      {
        clim_variables <- read_xlsx(glue::glue("{data_dir}/settings/bioclim_vars_definitions.xlsx"))
        clim_dir       <- list.files(glue::glue("{data_dir}/data/global_maps"), pattern = "*.tif", full.names = T)
        
        clim_stars <- lapply(clim_variables$Var_Name, function(x) {
          ind <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
          read_stars(clim_dir[ind]) %>% setNames(x)
        })
        
        for(ID in IDs) {
          if(!file.exists(glue::glue("{dir_out}/psaEnv/psaEnv_{ID}.rda")) & check_files) {
            
            load(glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda"))
            
            psaEnv <- lapply(clim_variables$Var_Name, function(x) {
              ind <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
              read_stars(clim_dir[ind]) %>% setNames(x) %>% st_extract(., psaCrds) %>% pull(x) %>% unlist()
            }) %>% Reduce("cbind",.) %>% as_tibble() %>% setNames(clim_variables$Var_Name) %>%
              bind_cols(read_stars(clim_dir[which(sapply(clim_dir, function(y) grepl("ETOPO", y, fixed = T)))]) %>%
                          st_set_crs(4326) %>% setNames("elev") %>% st_extract(., psaCrds) %>% pull("elev") %>% unlist() %>% as_tibble() %>% setNames("elev"))
            
            save(psaEnv, file = glue::glue("{dir_out}/psaEnv/psaEnv_{ID}.rda"))
          }
        }
      }
      
      
      
      
      ### Initiate map
      {
        dataset_rast <- read_stars(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
        
        mergeTab    <- class_defs %>% mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
          dplyr::select(Class_Code, lcov, Class_Plotlabel, Predicted_In_Past, Color_Code)
        
        initRast    <- dataset_rast %>% st_as_stars() %>% setNames("Landcover") %>% 
          mutate(Landcover = mergeTab$lcov[match(Landcover, mergeTab$Class_Code)]) %>%
          mutate(predict = mergeTab$Predicted_In_Past[match(Landcover, mergeTab$lcov)])
        
        ggplot() +
          geom_stars(data = initRast[1], mapping = aes(fill = as.factor(Landcover))) +
          geom_sf(psaVeg@psas %>% filter(Dataset_ID %in% IDs) %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="cyan") + 
          scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                            breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                            labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
          xlab("") + ylab("") +
          theme_light()
        
        ggsave(glue::glue("{dir_out}/summaryResults/map_init.png"), width = 20, height = 12, units = "cm")
        
      }
      
      ### Climate variable Maps
      {
        if(!file.exists(glue::glue("{dir_out}/summaryResults/climRegion.tiff")) & check_files) {
          clim_variables <- read_xlsx(glue::glue("{data_dir}/settings/bioclim_vars_definitions.xlsx"))
          clim_dir       <- list.files(glue::glue("{data_dir}/data/global_maps"), pattern = "*.tif", full.names = T)
          
          climRegion <- c(lapply(clim_variables$Var_Name, function(x) {
            ind     <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
            st_clim <- read_stars(clim_dir[ind]) %>% setNames(x)
            # st_warp(st_clim, initRast, use_gdal = T, method = "mode") %>% st_as_stars() %>% setNames(x)
            st_warp(st_clim, initRast, use_gdal = F, method = "near") %>% st_as_stars() %>% setNames(x)
          }) %>% do.call("c",.), st_warp(read_stars(clim_dir[grepl("ETOPO", clim_dir)]), initRast, use_gdal = F, method = "near") %>% setNames("elev")) %>% merge()
          write_stars(climRegion, glue::glue("{dir_out}/summaryResults/climRegion.tiff"))
        } else climRegion <- read_stars(glue::glue("{dir_out}/summaryResults/climRegion.tiff"))
      }
      
      
      
      
      
      
      ### RDA model PSAs
      threshold <- 0.15 ### 15 percentile of sample size across samples in one PSA
      {
        for(ID in IDs) {
          
          if(!file.exists(glue::glue("{dir_out}/psaOvlp/psaOvlp_{ID}.rda")) & check_files) {
            
            load(glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda")) ## psaCrds
            load(glue::glue("{dir_out}/psaEnv/psaEnv_{ID}.rda"))   ## psaEnv
            
            datst <- psaCrds %>% st_drop_geometry() %>% dplyr::select(-Dataset_ID) %>%
              bind_cols(psaEnv) %>% mutate(lcov = as.factor(lcov))
            
            # datst1 <- psaCrds %>% st_drop_geometry() %>% dplyr::select(-Dataset_ID) %>% 
            #   bind_cols(psaEnv) %>% mutate(lcov = as.factor(lcov))
            # 
            # datst <- datst1 %>% 
            #   mutate(across(starts_with(c("bio", "gdd", "elev")),~.*(1/(dist*dist))))
            
            
            sampleSize  <- datst %>% group_by(lcov) %>% summarise(sample = n())
            thrsh       <- quantile(sampleSize$sample, probs = threshold)
            
            dataset_sub <- sampleSize %>% mutate(sample = ifelse(sample>thrsh, thrsh, sample)) %>% 
              group_split(lcov) %>% lapply(., function(x) {
                datst %>% filter(lcov == x$lcov) %>% arrange(dist) %>% dplyr::select(-dist) %>% slice(1:x$sample)
              }) %>% Reduce("rbind", .)
            
            dataset_sub$source <- "calibration"
            
            #####
            
            pxlTab_init <- st_as_sf(climRegion) %>% mutate(lcov = (initRast %>% st_as_sf() %>% st_drop_geometry())[,1])
            
            # landcover classes that are kept fix:
            del    <- c(0,140,205,210,220)   ### 15,190,
            #del_lc <- paste0("LC_",del)
            
            envTab_human <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>% 
              # filter(lcov==15 | lcov==190) %>%
              filter(!(lcov %in% del)) %>%
              select(-pxlInd) %>% 
              as_tibble() %>% suppressWarnings()
            
            envTab_human <- envTab_human[,c(which(names(envTab_human) == "lcov"),1:(ncol(envTab_human)-1))]
            envTab_human$lcov <- as.factor(as.character(envTab_human$lcov))
            
            envTab_human$source <- "initmap"
            
            
            
            dataset_envHuman <- rbind(dataset_sub, envTab_human)
            #dataset_envHuman <- unique(dataset_envHuman)
            
            
            # z-standardization
            dataset_sub_z <- bind_cols(lcov = dataset_envHuman$lcov, 
                                       source = dataset_envHuman$source,
                                       dataset_envHuman %>% reframe(
                                         across(
                                           .cols  = c(-lcov, -source),
                                           .fns   = function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T),
                                           .names = "{col}"
                                         )
                                       )) %>% filter(apply(., 1, function(x) all(!is.na(x))))
            
            
            save(dataset_sub_z, file = glue::glue("{dir_out}/summaryResults/z_transformed_df_{ID}.rda"))
            
            
            dataset_sub_z_calibration <- dataset_sub_z %>% filter(source == "calibration") %>% select(-source)
            
            rdaOut <- rdaMod(dataset_sub_z_calibration, ID = ID)
            
            # # z-standardization
            # dataset_sub_z <- bind_cols(lcov = dataset_sub$lcov, 
            #                            dataset_sub %>% reframe(
            #                              across(
            #                                .cols  = c(-lcov),
            #                                .fns   = function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T),
            #                                .names = "{col}"
            #                              )
            #                            )) %>% filter(apply(., 1, function(x) all(!is.na(x))))
            
            # rdaOut <- rdaMod(dataset_sub_z, ID = ID)
            
            save(rdaOut, file = glue::glue("{dir_out}/summaryResults/rdaOut_{ID}.rda"))
            
            ### Overlap and centers
            psaOvlp <- make_psaOvlp(rdaOut)
            
            # plot(psaOvlp@centers)
            # plot(psaOvlp@densities)
            
            save(psaOvlp, file = glue::glue("{dir_out}/psaOvlp/psaOvlp_{ID}.rda"))
          }
        }
      }
      
      ### plot RDA output
      plotRDAsummary(ID, dir_out, class_defs)
      ggsave(glue::glue("{dir_out}/summaryResults/rda_summary.png"), width = 30, height = 22, units = "cm")
      
    
      
      
      ###################################
      ### LandCoverChange: Pixel flow ###
      ###################################
      
      
      
      ### Calculate LandCover Distribution Change (transport)
      ### Transport of proportional distribution
      IDs <- as.numeric(gsub(".rda", "", gsub("psaOvlp_", "", list.files(glue::glue("{dir_out}/psaOvlp")))))
      {
        if(!file.exists(glue::glue("{dir_out}/summaryResults/flowPSA.rda")) & check_files) {
          flowPSA <- lapply(IDs, function(ID) {
            load(glue::glue("{dir_out}/psaOvlp/psaOvlp_{ID}.rda"))
            probs <- psaOvlp@overlapp[,,1]
            colnames(probs) <- rownames(probs) <- psaOvlp@lcov
            
            psaLCover_psa <- psaLCover@landcover_ts %>% filter(Dataset_ID==ID)
            if(!all(colnames(probs) %in% glue::glue("LC_{names(psaLCover_psa)[-c(1:2)]}"))) {
              add <- (matrix(0, ncol = ncol(probs), nrow = nrow(psaLCover_psa)) %>%
                        as_tibble() %>% setNames(gsub("LC_","",colnames(probs))))[!colnames(probs) %in% 
                                                                                    glue::glue("LC_{names(psaLCover_psa)[-c(1:2)]}")]
              psaLCover_psa <- psaLCover_psa %>% bind_cols(add) %>% 
                dplyr::select(Dataset_ID, Age, gsub("LC_","",colnames(probs)))
            }
            
            psaLCover_psa <- psaLCover_psa %>% filter(!Age == 500)
            
            lapply(1:(nrow(psaLCover_psa)-1), function(r) {
              
              flow_in <- psaLCover_psa[r:(r+1),]
              
              trans <- transport(flow_in[1,] %>% dplyr::select(gsub("LC_", "", psaOvlp@lcov)) %>% as.numeric()/100, 
                                 flow_in[2,] %>% dplyr::select(gsub("LC_", "", psaOvlp@lcov)) %>% as.numeric()/100, 
                                 as.matrix(probs), fullreturn = TRUE)
              
              trans$default %>% as_tibble() %>%
                mutate(from = psaOvlp@lcov[from], to = psaOvlp@lcov[to], mass = mass*100) %>%
                mutate(ID = ID, .before = from) %>% mutate(year = flow_in$Age[2], .before = from)
              
            }) %>% Reduce("rbind", .)
            
          })
          names(flowPSA) <- IDs
          save(flowPSA, file = glue::glue("{dir_out}/summaryResults/flowPSA.rda"))
        } else load(glue::glue("{dir_out}/summaryResults/flowPSA.rda"))
      }
      
      ### Pixel change/flow [little time consuming]
      ## Method: one-by-one, highest overlap
      ## parallel computation
      cl <- parallel::makeCluster(min(c(length(IDs), parallel::detectCores()-25)))
      {
        parallel::clusterExport(cl, c("dir_out", "initRast"))
        
        invisible(parallel::clusterEvalQ(cl, {
          library(tidyr)
          library(glue)
          library(dplyr)
          library(tibble)
          library(sf); sf_use_s2(FALSE)
          library(stars)
          library(vegan)
          library(stringr)
          library(ggplot2)
          climRegion <- read_stars(glue::glue("{dir_out}/summaryResults/climRegion.tiff")) %>% st_as_stars()
          pxlTab     <- st_as_sf(climRegion) %>% mutate(lcov = (initRast %>% st_as_sf() %>% st_drop_geometry())[,1])
          load(glue::glue("{dir_out}/summaryResults/flowPSA.rda"))
          
          ##########################################################################################################
          ####### /!\ MAKE SURE THAT data_dir AND class_defs ARE THE SAME AS AT THE START OF THE SCRIPT /!\ ########
          data_dir    <- "/bioing/data/LandCoverChange/LandCoverChangeProject_data_new/"
          class_defs  <- readxl::read_excel(glue::glue("{data_dir}settings/landcover_class_definitions_mod.xlsx")) 
          mergeTab    <- class_defs %>% mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
            dplyr::select(Class_Code, lcov, Class_Plotlabel, Predicted_In_Past, Color_Code)
          check_files <- TRUE
          
          # Function to find inverse pairs
          find_inverse_pairs <- function(data, row_number) {
            target_row <- data[row_number, ]
            inverse_rows <- which(data$from == target_row$to & 
                                    data$to == target_row$from)
            return(data[inverse_rows, ])
          }
          
          ##########################################################################################################
          
        }))
        
        parallelFlow <- parallel::parLapply(cl, IDs, function(ID) {
          
          if(!file.exists(glue::glue("{dir_out}/psaChange/psaChange_{ID}.rda")) & check_files) {
            
            load(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
            load(glue::glue("{dir_out}/psaOvlp/psaOvlp_{ID}.rda"))
            
            IDs <- as.numeric(gsub(".rda", "", gsub("psaOvlp_", "", list.files(glue::glue("{dir_out}/psaOvlp")))))
            
            flow <- flowPSA[as.character(ID)][[1]]
            # flow <- flow %>% filter(!year == 500)
            lcov <- tibble(init = pxlTab$lcov)
            
            # initialization:
            currentRast      <- initRast
            lc_recon         <- NULL
            pxl_transport_df <- NULL
            
        
            # -------------------------------------------------------------------------------------------------
            
            # load(glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda")) ## psaCrds
            # load(glue::glue("{dir_out}/psaEnv/psaEnv_{ID}.rda"))   ## psaEnv
            # 
            # datst <- psaCrds %>% st_drop_geometry() %>% dplyr::select(-Dataset_ID) %>%
            #   bind_cols(psaEnv) %>% mutate(lcov = as.factor(lcov))
            # 
            # sampleSize  <- datst %>% group_by(lcov) %>% summarise(sample = n())
            # thrsh       <- quantile(sampleSize$sample, probs = threshold)
            # 
            # envTab_pixel <- sampleSize %>% mutate(sample = ifelse(sample>thrsh, thrsh, sample)) %>% 
            #   group_split(lcov) %>% lapply(., function(x) {
            #     datst %>% filter(lcov == x$lcov) %>% arrange(dist) %>% dplyr::select(-dist) %>% slice(1:x$sample)
            #   }) %>% Reduce("rbind", .)
            # 
            # 
            # 
            # # replace pixel with human influence in the init map with LC-classes with highest density:
            # 
            pxlTab_init <- st_as_sf(climRegion) %>% mutate(lcov = (initRast %>% st_as_sf() %>% st_drop_geometry())[,1])

            # landcover classes that are kept fix:
            del    <- c(0,140,205,210,220)   ### 15,190,
            #del_lc <- paste0("LC_",del)

            envTab_human <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>%
              # filter(lcov==15 | lcov==190) %>%
              filter(!(lcov %in% del)) %>%
              #select(-pxlInd) %>%
              as_tibble() %>% suppressWarnings()

            # envTab_human <- envTab_human[,c(which(names(envTab_human) == "lcov"),1:(ncol(envTab_human)-1))]
            # envTab_human$lcov <- as.factor(as.character(envTab_human$lcov))


            # 
            # dataset_envHuman <- rbind(envTab_pixel, envTab_human)
            # dataset_envHuman <- unique(dataset_envHuman)
            # 
            # 
            # 
            # # z-standardization
            # dataset_envHuman_z <- bind_cols(lcov = glue::glue("LC_{dataset_envHuman$lcov}"),    #lcov = dataset_envHuman$lcov,   
            #                            dataset_envHuman %>% reframe(
            #                              across(
            #                                .cols  = c(-lcov),
            #                                .fns   = function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T),
            #                                .names = "{col}"
            #                              )
            #                            )) %>% filter(apply(., 1, function(x) all(!is.na(x))))
            # 
            # 
            
            
            
            load(glue::glue("{dir_out}/summaryResults/z_transformed_df_{ID}.rda")) ## z-transformed dataset
            
            dataset_sub_z_initmap <- dataset_sub_z %>% filter(source == "initmap") %>% select(-source) %>% 
              mutate(pxlInd = envTab_human$pxlInd)
            
            # select only cropland and urban area classes:
            dataset_sub_z_initmap_human <- dataset_sub_z_initmap %>% filter(lcov==15 | lcov==190) 
            
            
            rda_scores_human <- as.data.frame(predict(psaOvlp@rda_mod[[1]], newdata = dataset_sub_z_initmap_human[,-ncol(dataset_sub_z_initmap_human)], type = "wa", scaling = 1)) %>% dplyr::select("RDA1", "RDA2") %>%
              mutate(pxlInd = dataset_sub_z_initmap_human$pxlInd, lcov = paste0("LC_",dataset_sub_z_initmap_human$lcov)) %>% 
              filter(!is.na(RDA1)) %>% st_as_sf(coords = c("RDA1", "RDA2"))  
            
            
            
            ###### plot
            
            {
              
              # colorTab <- class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>% 
              #   filter(!is.na(Merge_To_Class), !duplicated(Merge_To_Class), !is.na(Color_Code)) %>% filter(Merge_To_Class%in%c(as.numeric(gsub("LC_", "", psaOvlp@lcov)),15,190))
              # 
              # colorRDA <- colorTab %>% left_join(psaOvlp@scores %>% group_by(lcov) %>% summarise(sample = n()) %>% 
              #                                      mutate(Merge_To_Class = as.numeric(gsub("LC_", "", lcov))) %>% dplyr::select(Merge_To_Class, sample), by = join_by(Merge_To_Class))
              # 
              # lc_centers1 <- psaOvlp@centers %>% st_as_sf(coords = c("RDA1", "RDA2")) %>% st_drop_geometry()
              # 
              # rda_scores_human_plot <- as.data.frame(predict(psaOvlp@rda_mod[[1]], newdata = dataset_sub_z_initmap_human[,-ncol(dataset_sub_z_initmap_human)], type = "wa", scaling = 1)) %>% dplyr::select("RDA1", "RDA2") %>%
              #   mutate(pxlInd = dataset_sub_z_initmap_human$pxlInd, lcov = paste0("LC_",dataset_sub_z_initmap_human$lcov)) %>% 
              #   filter(!is.na(RDA1)) %>% 
              #   mutate(lcov = gsub("LC_","",lcov))
              # 
              # 
              # 
              # pl_human <- psaOvlp@scores %>% mutate(lcov_code = gsub("LC_", "", lcov)) %>%
              #   ggplot(aes(x = RDA1, y = RDA2, color = lcov_code)) +
              #   
              #   geom_point(data=rda_scores_human_plot, aes(x = RDA1, y = RDA2, color= lcov), size = 1, alpha = 0.8) +
              #   
              #   geom_point(size = 3, alpha = 0.8) +
              #   scale_color_manual(values = colorRDA$Color_Code, breaks = colorRDA$Merge_To_Class, name = "",
              #                      labels = glue::glue("{colorRDA$Variable_Name} - {colorRDA$Class_Plotlabel} (n = {colorRDA$sample})")) +
              #   
              #   geom_point(data=psaOvlp@centers, aes(x = RDA1, y = RDA2), size = 2, shape = 15, alpha = 0.8, color="black") +
              #   
              #   theme_light() +
              #   theme(legend.position = "bottom") +
              #   guides(color=guide_legend(nrow=3, byrow=TRUE))
              # 
              # pl_human
              
              }
            
            ######
            
            
            
            # determine distance to the center of each landcover class:
            
            lc_centers <- psaOvlp@centers %>% st_as_sf(coords = c("RDA1", "RDA2"))
            
            densities_human_center <- as.data.frame(st_distance(rda_scores_human, lc_centers) %>% st_drop_geometry()) %>%
              mutate(pxlInd = rda_scores_human$pxlInd,
                     old_lcov = rda_scores_human$lcov) %>% 
              setNames(c(psaOvlp@lcov, "pxlInd", "old_lcov")) 
            
            densities_human_center_DT <- data.table::data.table(densities_human_center %>% 
                                                                  select(-c(pxlInd, old_lcov))) 
            
            # determine the column with minimum value
            densities_human_center_DT[, new_lcov := colnames(.SD)[max.col(-.SD, ties.method = "first")]]
            
            converttable_center <- data.frame(pxlInd = densities_human_center$pxlInd,
                                              old_lcov = densities_human_center$old_lcov,
                                              new_lcov = densities_human_center_DT$new_lcov) 
            
            
            pxlTab_init_df <- pxlTab_init %>% rownames_to_column(var = 'pxlInd')
            
            
            
            
            
            
            
            
            # replace pixel with human influence in the init map with LC-classes with highest density:
            
            # pxlTab_init <- st_as_sf(climRegion) %>% mutate(lcov = (initRast %>% st_as_sf() %>% st_drop_geometry())[,1])
            # 
            # envTab_human <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>%
            #   filter(lcov==15 | lcov==190) %>%
            #   as_tibble() %>% suppressWarnings()

            # dataset_sub_z_human <- bind_cols(pxlInd = dataset_envHuman$pxlInd, lcov = glue::glue("LC_{dataset_envHuman$lcov}"),
            #                                  dataset_envHuman %>% reframe(
            #                                    across(
            #                                      .cols  = c(-lcov, -pxlInd),
            #                                      .fns   = function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T),
            #                                      .names = "{col}"
            #                                    )
            #                                  ))
            # 
            # dataset_sub_z_human_noNA <- lapply(1:ncol(dataset_sub_z_human), function(x) {
            #   out <- dataset_sub_z_human[,x]
            #   if(is.numeric(out[,1]) & any(is.na(out[,1]))) {
            #     out[is.na(out[,1]),1] <- median(out[is.na(out[,1]),1], na.rm = T)
            #   }
            #   out
            # }) %>% Reduce("cbind",.) %>% as_tibble()
            
            
            # rda_scores_human <- as.data.frame(predict(psaOvlp@rda_mod[[1]], newdata = dataset_sub_z_human_noNA[,-1], type = "wa", scaling = 1)) %>% dplyr::select("RDA1", "RDA2") %>%
            #   mutate(pxlInd = envTab_human$pxlInd, lcov = paste0("LC_",envTab_human$lcov)) %>% filter(!is.na(RDA1)) %>% st_as_sf(coords = c("RDA1", "RDA2"))  
            # 
            # densities_human <- as.data.frame(st_extract(psaOvlp@densities, rda_scores_human) %>% st_drop_geometry() %>% pull(densities)) %>%
            #   mutate(pxlInd = envTab_human$pxlInd,
            #          old_lcov = envTab_human$lcov) %>% 
            #   setNames(c(psaOvlp@lcov, "pxlInd", "old_lcov")) 
            #   
            # densities_human_noNA <- densities_human[complete.cases(densities_human), ]
            #   
            # densities_human_noNA_allZero <- densities_human_noNA %>% rownames_to_column() %>%  
            #   filter_at(vars(-c(rowname, pxlInd, old_lcov)), all_vars(. == 0))
            # 
            # # for those pixel where the densities aren't entirely zero:
            # 
            # densities_human_noNA_keep <- densities_human_noNA %>% filter(!(pxlInd %in% densities_human_noNA_allZero$pxlInd))
            # 
            # 
            # densities_human_noNA_DT <- data.table::data.table(densities_human_noNA_keep %>% 
            #                                                     select(-c(pxlInd, old_lcov)))
            # 
            # # determine the column with maximum value
            # densities_human_noNA_DT[, new_lcov := colnames(.SD)[max.col(.SD, ties.method = "first")]]
            # 
            # converttable_keep <- data.frame(pxlInd = densities_human_noNA_keep$pxlInd,
            #                                 old_lcov = densities_human_noNA_keep$old_lcov,
            #                                 new_lcov = densities_human_noNA_DT$new_lcov) 
            # 
            # 
            # # for those pixel where the densities are entirely zero:
            # 
            # rda_scores_human_allZero <- rda_scores_human %>% filter(pxlInd %in% densities_human_noNA_allZero$pxlInd)
            # 
            # lc_centers <- psaOvlp@centers %>% st_as_sf(coords = c("RDA1", "RDA2"))
            # 
            # densities_human_allZero <- as.data.frame(st_distance(rda_scores_human_allZero, lc_centers) %>% st_drop_geometry()) %>%
            #   mutate(pxlInd = rda_scores_human_allZero$pxlInd,
            #          old_lcov = rda_scores_human_allZero$lcov) %>% 
            #   setNames(c(psaOvlp@lcov, "pxlInd", "old_lcov")) 
            # 
            # densities_human_allZero_DT <- data.table::data.table(densities_human_allZero %>% 
            #                                                     select(-c(pxlInd, old_lcov)))
            # 
            # # determine the column with minimum value
            # densities_human_allZero_DT[, new_lcov := colnames(.SD)[max.col(-.SD, ties.method = "first")]]
            # 
            # converttable_allZero <- data.frame(pxlInd = densities_human_allZero$pxlInd,
            #                                    old_lcov = densities_human_allZero$old_lcov,
            #                                    new_lcov = densities_human_allZero_DT$new_lcov) 
            # 
            # 
            # # combine tables
            # converttable_human <- rbind(converttable_keep, converttable_allZero)
            # 
            # 
            # pxlTab_init_df <- pxlTab_init %>% rownames_to_column(var = 'pxlInd')
            # 
            # # replace old human impact classes with new ones:
            # 
            # pxlTab_init_df_keep <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>% 
            #   filter(!(pxlInd %in% converttable_human$pxlInd)) 
            # 
            # pxlTab_init_df_replace <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>% 
            #   filter(pxlInd %in% converttable_human$pxlInd) %>% 
            #   left_join(converttable_human, by="pxlInd") %>% 
            #   mutate(lcov_corr = as.numeric(gsub("LC_", "", new_lcov))) %>% 
            #   select(pxlInd:elev, lcov_corr) %>% 
            #   rename(lcov = lcov_corr)
            # 
            # 
            # pxlTab_new <- rbind(pxlTab_init_df_keep, pxlTab_init_df_replace)
            # pxlTab_new <- pxlTab_new %>% as_tibble() %>% arrange(as.numeric(pxlInd))
            # 
            
            #############################################################################
            # 
            # # determine distance to the center of each landcover class:
            # 
            # lc_centers <- psaOvlp@centers %>% st_as_sf(coords = c("RDA1", "RDA2"))
            # 
            # densities_human_center <- as.data.frame(st_distance(rda_scores_human, lc_centers) %>% st_drop_geometry()) %>%
            #   mutate(pxlInd = rda_scores_human$pxlInd,
            #          old_lcov = rda_scores_human$lcov) %>% 
            #   setNames(c(psaOvlp@lcov, "pxlInd", "old_lcov")) 
            # 
            # densities_human_center_DT <- data.table::data.table(densities_human_center %>% 
            #                                                        select(-c(pxlInd, old_lcov)))
            # 
            # # determine the column with minimum value
            # densities_human_center_DT[, new_lcov := colnames(.SD)[max.col(-.SD, ties.method = "first")]]
            # 
            # converttable_center <- data.frame(pxlInd = densities_human_center$pxlInd,
            #                                    old_lcov = densities_human_center$old_lcov,
            #                                    new_lcov = densities_human_center_DT$new_lcov) 
            # 
            # 
            # pxlTab_init_df <- pxlTab_init %>% rownames_to_column(var = 'pxlInd')
            # 
            
            
            # replace old human impact classes with new ones:
            
            pxlTab_init_df_keep <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>% 
              filter(!(pxlInd %in% converttable_center$pxlInd)) 
            
            pxlTab_init_df_replace <- pxlTab_init %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>% 
              filter(pxlInd %in% converttable_center$pxlInd) %>% 
              left_join(converttable_center, by="pxlInd") %>% 
              mutate(lcov_corr = as.numeric(gsub("LC_", "", new_lcov))) %>% 
              select(pxlInd:elev, lcov_corr) %>% 
              rename(lcov = lcov_corr)
            
            
            pxlTab_new <- rbind(pxlTab_init_df_keep, pxlTab_init_df_replace)
            pxlTab_new <- pxlTab_new %>% as_tibble() %>% arrange(as.numeric(pxlInd))
            
            
            lcov_human <- tibble(init_human = pxlTab_new$lcov)
            lcov       <- bind_cols(lcov, lcov_human) %>% suppressMessages()
            
            # convert new initmap_human into raster
            initRast_human <- currentRast[1]
            initRast_human$Landcover <- matrix(c(lcov_human)[[1]], ncol = dim(initRast_human$Landcover)[1], nrow = dim(initRast_human$Landcover)[2])
            # initRast_human
            initRast_human <- initRast_human %>% setNames("Landcover")
            
            currentRast <- initRast_human
            
            # -------------------------------------------------------------------------------------------------
            
            initRast_human_plot <- initRast_human %>% setNames("LandcoverChange")
            
            # ggplot() +
            #   #geom_stars(data = initRast_human_plot, mapping = aes(fill = as.factor(Landcover))) +
            #   geom_stars(data = initRast_human_plot, mapping = aes(fill = as.factor(LandcoverChange))) +
            #   #facet_wrap(~attributes) +
            #   scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
            #                     breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
            #                     labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
            #   theme_void()
            #
            # ggsave(glue::glue("{dir_out}/summaryResults/psaChange_maps_{ID}_1ka.png"), width = 30, height = 20, units = "cm")
            
            
            # save init map with replaced human influence:
            ggplot() +
              geom_stars(data = initRast_human_plot, mapping = aes(fill = as.factor(LandcoverChange))) +
              geom_sf(psaVeg@psas %>% filter(Dataset_ID %in% IDs) %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="cyan") + 
              scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                                breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                                labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
              xlab("") + ylab("") +
              theme_light()
            
            ggsave(glue::glue("{dir_out}/summaryResults/map_init_human.png"), width = 20, height = 12, units = "cm")
            
            # -------------------------------------------------------------------------------------------------
            
            # landcover classes that are kept fix:
            del    <- c(0, 15,140,190,205,210,220)
            del_lc <- paste0("LC_",del)
            
            # number of movable pixel in init map:
            pxlTab_human <- st_as_sf(climRegion) %>% mutate(lcov = (currentRast %>% st_as_sf() %>% st_drop_geometry())[,1])
            
            # pxlTab_move <- pxlTab[!pxlTab$lcov %in% del, ]
            pxlTab_move <- pxlTab_human[!pxlTab_human$lcov %in% del, ]
            pxlTab_move_tab <- table(pxlTab_move$lcov)
            
            
           
            
            for(y in unique(flow$year)) {
              
              cat(glue::glue("\rpixel flow: year {y} ({which(y==unique(flow$year))} of {length(unique(flow$year))})"))
              
              pxlTab <- st_as_sf(climRegion) %>% mutate(lcov = (currentRast %>% st_as_sf() %>% st_drop_geometry())[,1])
              
              flow_y_init <- flow %>% filter(year == y) %>%
                full_join(tibble(from = glue::glue("LC_{unlist(lcov[,ncol(lcov)])}")) %>% filter(!from %in% del_lc) %>% 
                            group_by(from) %>% summarise(pxl = n()), by = 'from') %>%
                mutate(movable_Pxl = sum(pxlTab_move_tab)) 
             
               
              # flow_y_init %>% select(ID,year,from,pxl) %>% distinct()
              
              flow_y_init1       <- flow_y_init
              flow_y_init1_clean <- NULL
              
              
              for(ff in 1:dim(flow_y_init1)[1]){
                
                row_sub <- flow_y_init1[1, ]
                
                row_sub$mass <- ifelse(row_sub$from == row_sub$to, 0, row_sub$mass) 
                
                inverse_row_sub <- find_inverse_pairs(flow_y_init1, 1)
                
                if(dim(inverse_row_sub)[1] > 0){
                  
                  pairing_df <- rbind(row_sub, inverse_row_sub)
                  
                  pairing_df <- pairing_df %>% arrange(mass)
                  
                  pairing_df[2,"mass"] <- pairing_df[2,"mass"] - pairing_df[1,"mass"]
                  pairing_df[1,"mass"] <- 0
                  
                  flow_y_init1 <- flow_y_init1 %>% filter(!(from %in% pairing_df$from & to %in% pairing_df$to))
                  
                }else{
                  
                  pairing_df <- row_sub
                  
                  flow_y_init1 <- flow_y_init1 %>% filter(!(from == pairing_df$from & to == pairing_df$to))
                  
                }
                
                flow_y_init1_clean <- rbind(flow_y_init1_clean, pairing_df)
                
                if(dim(flow_y_init1)[1] == 0){break}
                
              }
              
              
              flow_y <- flow_y_init1_clean %>%
                mutate(nrPxl = round(movable_Pxl*(mass/100),0)) %>% filter(nrPxl>0) %>% arrange(pxl) %>% 
                mutate(pxl = ifelse(is.na(pxl),0,pxl)) %>% 
                mutate(useAll = ifelse(any(nrPxl > round(pxl+pxl*0.1,0)), TRUE, FALSE))
              
              
              
              ovlpTab <- lapply(1:nrow(flow_y), function(r) {
                
                envTab <- pxlTab %>% st_drop_geometry() %>% rownames_to_column(var = 'pxlInd') %>%
                  filter(lcov==as.numeric(gsub("LC_", "", flow_y$from[r]))) %>%
                  as_tibble() %>% suppressWarnings()
                  
                ### project into RDA space
                ovlp <- psaOvlp@densities[,,,unlist(flow_y[r,c("from", "to")])] %>% st_apply(., 1:2, function(x) x[1]*x[2])
                {
                  dist_rast <- psaOvlp@densities[,,,unlist(flow_y[r,"to"])] %>% setNames("distance")
                  dist_vect <- dist_rast %>% st_as_sf()
                  dist_vect$distance <- st_distance(dist_vect, dist_vect[which.max((dist_vect %>% st_drop_geometry())[,1])[1],])
                  dist_rast[[1]][] <- dist_vect$distance[,1]
                }
                
                
           #     if(dim(envTab)[1] > 10){
                  
                  # z-standardization
                  dataset_sub_z <- bind_cols(pxlInd = envTab$pxlInd, lcov = glue::glue("LC_{envTab$lcov}"), 
                                             envTab %>% reframe(
                                               across(
                                                 .cols  = c(-lcov, -pxlInd),
                                                 .fns   = function(x) (x - mean(x, na.rm = T))/sd(x, na.rm = T),
                                                 .names = "{col}"
                                               )
                                             ))
                  
                  dataset_sub_z_noNA <- lapply(1:ncol(dataset_sub_z), function(x) {
                    out <- dataset_sub_z[,x]
                    if(is.numeric(out[,1]) & any(is.na(out[,1]))) {
                      out[is.na(out[,1]),1] <- median(out[is.na(out[,1]),1], na.rm = T)
                    }
                    out
                  }) %>% Reduce("cbind",.) %>% as_tibble()
                  
                  if(nrow(dataset_sub_z_noNA)>0) {
                    as.data.frame(predict(psaOvlp@rda_mod[[1]], newdata = dataset_sub_z_noNA[,-1], type = "wa", scaling = 1)) %>% dplyr::select("RDA1", "RDA2") %>%
                      mutate(pxlInd = envTab$pxlInd, lcov = paste0("LC_",envTab$lcov)) %>% filter(!is.na(RDA1)) %>% st_as_sf(coords = c("RDA1", "RDA2")) %>%
                      mutate(ovlp = st_extract(ovlp,.) %>% st_drop_geometry() %>% pull(densities),
                             dist = st_extract(dist_rast,.) %>% st_drop_geometry() %>% pull(distance)) %>% st_drop_geometry() %>%
                      mutate(ovlp = ifelse(is.na(ovlp), 0, ovlp), ovlp = scales::rescale(max(ovlp)-ovlp, c(0,1))) %>%
                      mutate(from = flow_y$from[r], to = flow_y$to[r]) %>% dplyr::select(pxlInd, from, to, dist, ovlp, lcov)
                  } else NULL
       #         } # end if (dim(envTab))
                  
              }) %>% Reduce("rbind",.) %>% 
                group_by(from, to) %>% arrange(desc(ovlp)) %>% ungroup()
              
              
              # outTab <- ovlpTab %>% group_split(from, to) %>%
              #   lapply(., function(x) {
              #     x1 <- x %>% filter(lcov == from) %>% arrange(desc(ovlp), dist) #%>% slice(1:(flow_y %>% filter(from==x$from[1], to==x$to[1]) %>% pull(nrPxl)))
              #     x2 <- x %>% filter(!(lcov == from | lcov %in% del_lc)) %>% arrange(desc(ovlp), dist) #%>% slice(1:(flow_y %>% filter(from==x$from[1], to==x$to[1]) %>% pull(nrPxl)))
              #     rbind(x1, x2) %>% slice(1:(flow_y %>% filter(from==x$from[1], to==x$to[1]) %>% pull(nrPxl)))
              #   }) %>% Reduce("rbind",.) %>%
              #   filter(!duplicated(pxlInd)) %>% ungroup()
              
              
              outTab <- ovlpTab %>% group_split(from, to) %>%
                lapply(., function(x) {
                  x %>% arrange(desc(ovlp), dist) %>% slice(1:(flow_y %>% filter(from==x$from[1], to==x$to[1]) %>% pull(nrPxl)))
                }) %>% Reduce("rbind",.) %>%
                filter(!duplicated(pxlInd)) %>% ungroup()
              
              
              outTab_sums <- outTab %>% group_by(from, to) %>% summarise(nrPxl_found = n(), .groups = 'drop')
              
              flowRemain <- flow_y %>% full_join(outTab_sums, join_by(from, to)) %>%
                mutate(nrPxl_found = ifelse(is.na(nrPxl_found), 0, nrPxl_found)) %>% 
                mutate(discr = nrPxl - nrPxl_found) %>% filter(discr > 0) %>% dplyr::select(from, to, discr) %>% 
                arrange(discr) %>% suppressMessages()
              
              flowRemain_thres <- round(sum(flowRemain$discr) * 0.05, 0)
              
              # -------------------------------------------------------------------------------------------------
              
              if(dim(flowRemain)[1] > 0){  
                
                repetition <- 0
                
                repeat {
                  
                  repetition <- repetition + 1
                  # cat(paste("\rrepetition:",repetition))
                  
                  outNew <- NULL 
                  
                  for(p in 1:dim(flowRemain)[1]){
                    
                    ovlpTab_p <- ovlpTab %>%
                      filter(!(pxlInd %in% outTab$pxlInd)) %>%
                      filter(!(lcov %in% del_lc)) %>% 
                      filter(from == flowRemain$from[p], to == flowRemain$to[p]) %>% 
                      arrange(desc(ovlp), dist) %>%
                      slice(1:min((flowRemain %>% filter(from==flowRemain$from[p], to==flowRemain$to[p]) %>% pull(discr)), flowRemain_thres))
                    
                    outNew <- rbind(outNew, ovlpTab_p)
                    
                    ovlpTab <- ovlpTab %>% filter(!(pxlInd %in% ovlpTab_p$pxlInd))
                    
                  }
                  
                  if(nrow(outNew)<1 || is.null(outNew)) break else outTab <- outTab %>% bind_rows(outNew)
                  
                  
                  outTab_sums <- outTab %>% group_by(from, to) %>% summarise(nrPxl_found = n(), .groups = 'drop')
                  
                  flowRemain <- flow_y %>% full_join(outTab_sums, join_by(from, to)) %>%
                    mutate(nrPxl_found = ifelse(is.na(nrPxl_found), 0, nrPxl_found)) %>% 
                    mutate(discr = nrPxl - nrPxl_found) %>% filter(discr > 0) %>% dplyr::select(from, to, discr) %>% 
                    arrange(discr) %>% suppressMessages()
                  
                  if(nrow(flowRemain)<1) break
                  
                } # end repeat    
                
              } # end if
              
              outTab_moved <- outTab %>% group_by(from, to) %>% summarise(N = n(), .groups = 'drop') 
              
              flow_y_moved <- full_join(flow_y, outTab_moved, by=c("from", "to")) %>% mutate(diff = nrPxl - N)
              
              # -------------------------------------------------------------------------------------------------
              
              pxl_transport_df <- rbind(pxl_transport_df, flow_y_moved)
              
              new_lcov <- lcov[,ncol(lcov)] 
              new_lcov[as.numeric(outTab$pxlInd),] <-  as.numeric(gsub("LC_", "", outTab$to))
              colnames(new_lcov) <- paste0("year_",y)
              
              lcov <- bind_cols(lcov, new_lcov) %>% suppressMessages()
              
              rast_tmp <- currentRast[1]
              rast_tmp$Landcover <- matrix(c(new_lcov)[[1]], ncol = dim(rast_tmp$Landcover)[1], nrow = dim(rast_tmp$Landcover)[2])
              rast_tmp
              rast_tmp <- rast_tmp %>% setNames("Landcover")
              
              currentRast <- rast_tmp
              
            }
            
            # save pixelflow table:
            openxlsx::write.xlsx(pxl_transport_df, file = glue::glue("{dir_out}/summaryResults/pixelflow_table_{ID}.xlsx"), asTable = TRUE)
            
            
            
            
            # modified psaChange without init for stratplot (only init_human):
            lcov_tmp <- lcov[,-1]
            
            psaChange_strat <- lapply(1:ncol(lcov_tmp), function(x) {
              rast <- initRast[1]
              rast$Landcover <- matrix(c(lcov_tmp[,x])[[1]], ncol = dim(rast$Landcover)[1], nrow = dim(rast$Landcover)[2])
              rast
            }) %>% Reduce("c",.) %>% setNames(glue::glue("year_{c('Init', unique(flow$year))}")) %>% merge()
            #setNames(names(lcov_tmp)) %>% merge()
            
            
            
            
            
            # save psaChange:
            psaChange <- lapply(1:ncol(lcov), function(x) {
              rast <- initRast[1]
              rast$Landcover <- matrix(c(lcov[,x])[[1]], ncol = dim(rast$Landcover)[1], nrow = dim(rast$Landcover)[2])
              rast
            }) %>% Reduce("c",.) %>% setNames(glue::glue("year_{c('Init', 'Init_human', unique(flow$year))}")) %>% merge()
            #setNames(names(lcov)) %>% merge()
            
            save(psaChange, file = glue::glue("{dir_out}/psaChange/psaChange_{ID}.rda"))
            
            
            # rast_for_plot <- rast_tmp %>% setNames("LandcoverChange")
            # 
            # ggplot() +
            #   #geom_stars(data = rast_for_plot, mapping = aes(fill = as.factor(Landcover))) +
            #   geom_stars(data = rast_for_plot, mapping = aes(fill = as.factor(LandcoverChange))) +
            #   #facet_wrap(~attributes) +
            #   scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
            #                     breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
            #                     labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
            #   theme_void()
            # 
            # ggsave(glue::glue("{dir_out}/summaryResults/psaChange_maps_{ID}_1ka.png"), width = 30, height = 20, units = "cm")
            
            # -------------------------------------------------------------------------------------------------
            
            # map based on PSA:
            
            psaChange_for_plot <- psaChange %>% setNames("LandcoverChange") 
            
            ggplot() +
              geom_stars(data = psaChange_for_plot, mapping = aes(fill = as.factor(LandcoverChange))) +
              facet_wrap(~attributes) +
              scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                                breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                                labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
              theme_void()
            
            ggsave(glue::glue("{dir_out}/summaryResults/psaChange_maps_{ID}.png"), width = 30, height = 20, units = "cm")
            
            # -------------------------------------------------------------------------------------------------
            
            # map with flooded regions based on PSA:
            
            psaChange_for_plot_flood <- psaChange_for_plot
            
            flooded_lc <- c(165, 180, 210)
            
            `%notin%` <- Negate(`%in%`)
            
            for(k in 1:dim(psaChange_for_plot_flood)[3]){
              psaChange_for_plot_flood[[1]][,,k][which(psaChange_for_plot_flood[[1]][,,k] %notin% flooded_lc)] <- 9999
            }
            
            ggplot() +
              geom_stars(data = psaChange_for_plot_flood, mapping = aes(fill = as.factor(LandcoverChange))) +
              facet_wrap(~attributes) +
              scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                                breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                                labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
              theme_void()
            
            ggsave(glue::glue("{dir_out}/summaryResults/psaChange_maps_flooded_{ID}.png"), width = 30, height = 20, units = "cm")
            
            # -------------------------------------------------------------------------------------------------
            
            # stratigraphic plot:
            
            # timeslices_psa <- names(psaChange)
            timeslices_psa <- c(0, unique(flow$year))
            
            lc_ts_rast_plot <- NULL
            
            
            for(k in 1:dim(psaChange_strat[[1]])[3]){
              
              # print(paste0(k,"/",dim(psaChange[[1]])[3]))
            
              pxl_stat <- table(psaChange_strat[[1]][,,k])
              
              pxl_stat_lc <- pxl_stat[!(names(pxl_stat) %in% del)] %>% 
                data.frame() %>% 
                mutate(Percentage = (Freq / sum(Freq)) * 100) %>% 
                mutate(Var1 = paste0("LC_",Var1)) %>% 
                select(-Freq)
              
              colnames(pxl_stat_lc)[1] <- "name"
              
              # age_stat_lc1 <- unlist(lapply(strsplit(split="_", timeslices_psa),function(x)return(x[k]))) 
              # age_stat_lc  <- gsub(".year","",age_stat_lc1)
              # age_stat_lc  <- ifelse((age_stat_lc == "Init" | age_stat_lc == "init"), 0, as.numeric(age_stat_lc))
              age_stat_lc <- timeslices_psa[k]
              
              
              pxl_stat_lc_ts <- data.frame(Dataset_ID = ID, Age = age_stat_lc, pxl_stat_lc)
              
              lc_ts_rast_plot <- rbind(lc_ts_rast_plot, pxl_stat_lc_ts)
              
            } # end for [k]
            
            lc_ts_rast_plot <- lc_ts_rast_plot %>% mutate(source = "Raster")
            
            # -------------------------------------------------------------------------------------------------
            
            load(glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
            
            lc_ts_df <- as.data.frame(psaLCover@landcover_ts)
            colnames(lc_ts_df)[-c(1:2)] <- paste0("LC_",names(lc_ts_df[,-c(1:2)]))
            
            psa_lc_index <- unique(lc_ts_df$Dataset_ID)
            
            psa_lc <- psa_lc_index[which(psa_lc_index == ID)] 
            
            lc_ts_df_sub <- lc_ts_df %>% filter(Dataset_ID == psa_lc)
            
            lc_ts_df_plot  <- lc_ts_df_sub %>% 
              pivot_longer(cols = -c(Dataset_ID, Age), values_to = "Percentage") %>% 
              mutate(source = "PSA")
            
            # -------------------------------------------------------------------------------------------------
            
            lc_stratplot_df <- rbind(lc_ts_df_plot, lc_ts_rast_plot)
            
            lc_stratplot_df$name_sort = factor(lc_stratplot_df$name, levels=str_sort(unique(lc_stratplot_df$name), numeric = TRUE)) 
            
            p3 <- ggplot(lc_stratplot_df %>% filter(!is.na(Percentage))) +
              geom_path(aes(x = Percentage, y = Age, group = as.factor(source), color = as.factor(source))) +
              scale_colour_manual(values = c("#E69F00", "midnightblue"), name = "") +
              facet_wrap(~name_sort, nrow = 4, ncol = 5) +
              scale_y_reverse(breaks=seq(0,max(lc_ts_df_plot$Age),2000)) +
              labs(title = glue::glue("Pollen source area ID: {psa_lc}")) +
              theme_light() +
              theme(legend.position = "top")
            
            ggsave(p3, file = glue::glue("{dir_out}/summaryResults/stratplot_LC_PSA_",ID,"_raster.png"), width = 20, height = 20, units = "cm", dpi = 300)
            
          } # end if (file exists)
          
        })
        
        parallel::stopCluster(cl)
      }
      
      # ### Merge Maps
      # {
      #   psaChangeTab  <- tibble(Path = list.files(glue::glue("{dir_out}/psaChange")),
      #                           Dataset_ID = as.numeric(gsub(".rda", "", gsub("psaChange_", "", Path))))
      #   
      #   lcovChangeMap <- mergeMaps(psaLCover, initRast, psaChangeTab, path = glue::glue("{dir_out}/psaChange"))
      #   
      #   save(lcovChangeMap, file = glue::glue("{dir_out}/summaryResults/lcovChangeMap.rda"))
      #   
      #   ggplot() +
      #     geom_stars(data = lcovChangeMap, mapping = aes(fill = as.factor(LandcoverChange))) +
      #     facet_wrap(~attributes) +
      #     scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
      #                       breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
      #                       labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
      #     theme_void()
      #   
      #   ggsave(glue::glue("{dir_out}/summaryResults/mergedMaps.png"), width = 30, height = 20, units = "cm")
      #   
      #   
      #   
      #   lcovChangeMap_flood <- lcovChangeMap
      #   
      #   flooded_lc <- c(165, 180, 210)
      #   
      #   `%notin%` <- Negate(`%in%`)
      #   
      #   for(k in 1:dim(lcovChangeMap_flood)[3]){
      #     
      #     lcovChangeMap_flood[[1]][,,k][which(lcovChangeMap_flood[[1]][,,k] %notin% flooded_lc)] <- 9999
      #     
      #   }
      #   
      #   
      #   ggplot() +
      #     geom_stars(data = lcovChangeMap_flood, mapping = aes(fill = as.factor(LandcoverChange))) +
      #     facet_wrap(~attributes) +
      #     scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
      #                       breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
      #                       labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
      #     theme_void()
      #   
      #   ggsave(glue::glue("{dir_out}/summaryResults/mergedMaps_flooded.png"), width = 30, height = 20, units = "cm")
      #   
      #   
      # }
      
      
      
    },error=function(e){}
  )
  
  
} # end for (rr)







# 
# psa_centers      <- psaVeg@psas %>% filter(Dataset_ID==ID)
# 
# 
# psa_centers1 <- sp::SpatialPointsDataFrame(coords=XY, data=psa_centers$geometry)
# psa_centers1 <- as.data.frame(psa_centers$geometry)
# 
# geom_point(data = grid_data, aes(x = grid_data$long, y = grid_data$lat), size = 3, 
#            shape = 23, fill = "orange")
# 
# 


