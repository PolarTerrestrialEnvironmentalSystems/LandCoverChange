#################################
### Template: LandCoverChange ###
### 22.10.2023                ###
### author: Simeon Lisovski   ###
###         Peter Ewald       ###
###         Thomas Böhmer     ###
#################################

library(tensorflow)
library(rgee)
library(geojson)
library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(readxl)
library(sf); sf_use_s2(FALSE)
library(stars)
library(ggplot2)
library(patchwork)
library(vegan)
library(scales)
library(transport)


source("functions/LandCoverChange.R")

#############
### Setup ###
#############

check_files <- TRUE
data_dir    <- "~/Documents/LandCoverChangeProject_data_new"
project_dir <- "~/Documents/LandCoverChange_test"

#### Google Earth Engine ####
{
  library(rgee)
  gee_user <- "simeon.lisovski@gmail.com"
  # gee_user <- "th.boehmer77@gmail.com"
  ee_Initialize(drive = TRUE, user = gee_user)
 
  # rgee::ee_install_set_pyenv(
  #  py_env = "ree",
  #  py_path = "~/anaconda3/envs/rgee/bin/python3.12"
  # )
}


### PSAs
metadata_path <- "/Volumes/projects/bioing/data/LandCoverChange/" 
psa_metadata  <- read_csv(glue::glue("{metadata_path}/PSA_locations_northern_hemisphere.csv"), show_col_types = F)


for(rr in sample(1:nrow(psa_metadata), 10)) {#1:nrow(psa_metadata)) {
  
  cat("\n")
  print(paste0("PSA: ", psa_metadata$Dataset_ID[rr]," (",rr,"/",nrow(psa_metadata),")"))
  
  ID <- psa_metadata$Dataset_ID[rr]
  
  ####################################
  ### Vegetation and Initial Setup ###
  ####################################
  
  #### 1.1 Directory Setup ####
  {
    current_psa  <- psa_metadata[rr,]
    psa_long     <- with(list(Long = current_psa$Longitude), glue::glue("{round(abs(Long),2)}{ifelse(Long>0, 'E', 'W')}"))
    psa_lat      <- with(list(Lat  = current_psa$Latitude),  glue::glue("{round(abs(Lat),2)}{ifelse(Lat>0, 'S', 'N')}"))

    dataFiles <- tibble(path = list.files(data_dir, recursive = T)) %>%
      mutate(folder = sapply(strsplit(gsub("data/", "", path), "/"), function(x)  x[[1]]),
             file = sapply(strsplit(gsub("data/", "", path), "/"), function(x)  x[[2]])) %>%
      dplyr::select(folder, file)
    # knitr::kable(dataFiles)
    
    dir_out <- glue::glue("{project_dir}/PSA_{current_psa$Dataset_ID}__{psa_long}_{psa_lat}/")
    if(!file.exists(dir_out)) {
      dir.create(dir_out)
      createFolders(dir_out)
      fs::dir_tree(path = dir_out, recurse = TRUE)
    }
  }
  
  #### 1.2 Spatial extent and resolution
  {
    ## PSA polygon and output resolution
    psa <- st_point(current_psa %>% dplyr::select(Longitude, Latitude) %>% unlist()) %>% st_sfc(crs = 4326)
    
    ## radius PSA [m]
    radius <- current_psa %>% pull(Pollen_Source_Radius)
    
    ## buffer to include PSAs for calibration [km]
    calibration_buffer  <- 1000            
    
    ## Native resolution = 309.2208 meter, possible options: 300, 3000, 5000, 10000 
    resolution   <- 300    
    
    ## Pixel per class for RDA
    pxl_extract_PSA        <- 800    # max. points for RDA to search within the PSA   
    pxl_extract_PSA_buffer <- 50     # min. points for RDA to search in buffer rings if they are not found in the PSA 
    
    ## Maximal spatial buffer around PSA
    max_buffer             <- 1000        # maximum buffer for search
  }
  
  #### 1.3 Create vegetation object for PSA
  {
    # load data, filter pollen source areas in region and/or calibration subset
    if(!file.exists(glue::glue("{dir_out}/summaryResults/psaVeg.rda")) & check_files) {
      psaVeg <- make_psaVeg(psa, radius, calibration_buffer, resolution, proj = "laea",
                            veg_folder = file.path(data_dir, "data", "vegetation_cover"), 
                            fc = 1/500, dt = 500, k = 5, nrCores = 9)
      save(psaVeg, file = glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
      
      plot(psaVeg)
      ggsave(glue::glue("{dir_out}/summaryResults/Map.png"), width = 10, height = 10, units = "cm")
      
    } else load(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
    
    ## Plot often not usefull since data only contains surface information for recent time slice
    # plotVegSite(psaVeg, ID = unique(psaVeg@vegetation %>% filter(!Modern) %>% pull(Dataset_ID))[6], n_taxa = 10, nrow = 2, ncol = 5)
  }
  
  
  ###################################
  ### Landcover                   ###
  ###################################
  
  #### 2.1 Modern LandCover share
  {
    if(!file.exists(glue::glue("{dir_out}/psaEnv/modernVegetation.rda"))) {
      ### modern time threshold for calibration
      modern_age_threshold <- 1000 
      
      modernVegetation <- psaVeg@vegetation %>% filter((Interp & Age <= modern_age_threshold) | Modern) %>%
        dplyr::select(-c(Interp, Age, Modern)) %>%
        group_by(Dataset_ID) %>%
        summarise(across(everything(), function(x) mean(x, na.rm = TRUE))) %>% filter(!apply(., 1, function(x) all(is.nan(x[-1]))))
      save(modernVegetation, file = glue::glue("{dir_out}/psaEnv/modernVegetation.rda"))
    } else load(glue::glue("{dir_out}/psaEnv/modernVegetation.rda"))
    
    class_defs <- readxl::read_excel(glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx")) 
  } ## end 2.1
  
  #### 2.2 CCI LandCover map of PSA
  {
    roi = sf_as_ee(psaVeg@regionMap %>% st_transform(4326) %>% st_shift_longitude())
    dataset <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff")) & check_files) {
      dataset_rast <- getLCCfromGEE(ee_dataset = dataset, sf_roi = psaVeg@regionMap, y_max = 500, buffer = 0.1)
      write_stars(dataset_rast, glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
    } else dataset_rast <- read_stars(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
       
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
  } ## end 2.2
  
  #### 2.3 Modern LandCover shares in PSA
  { 
    
    ### Classes present in buffer to use in calibration/prediction
    roi_maxbuffer = sf_as_ee(psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_buffer(.$'Pollen_Source_Radius [m]' + (max_buffer*1000)) %>% 
                             st_transform(4326) %>% st_shift_longitude())    
    
    classes_buffer  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
      collection = roi_maxbuffer,
      reducer    = ee$Reducer$count()$frequencyHistogram(),
      scale      = dataset$projection()$nominalScale()$getInfo()
    )$getInfo()
    
    filterLCover = tibble(Class_Code = as.numeric(names(classes_buffer$features[[1]][[4]][[6]])),
                               count = unlist(classes_buffer$features[[1]][[4]][[6]])) %>%
      left_join(class_defs %>% mutate(Merge_To_Class = ifelse(!is.na(Merge_To_Class), Merge_To_Class, Class_Code)) %>%
                  filter(Include_Class_In_Calibration) %>% dplyr::select(Class_Code, Merge_To_Class), by = 'Class_Code') %>%
      filter(!is.na(Merge_To_Class)) %>% group_by(Merge_To_Class) %>% summarise(count = sum(count)) %>% filter(count>=pxl_extract_PSA_buffer) %>%
      rename(Lcov = Merge_To_Class)
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/lcovShares.rda"))) {
      lcovShares <- modernLC(psaVeg, psaIDs = unique(modernVegetation$Dataset_ID), class_defs = class_defs, class_keep = filterLCover, resolution = resolution)
      save(lcovShares, file = glue::glue("{dir_out}/summaryResults/lcovShares.rda"))
    } else load(glue::glue("{dir_out}/summaryResults/lcovShares.rda"))
    
    ## plot
    # {
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
    # }
    
  } ## end 2.3
  
  #### 2.4 Translate Vegetation to Landcover
  {
    if(!file.exists(glue::glue("{dir_out}/summaryResults/veg_lc_model.rda")) & check_files) {
      
      veg_lc_model <- regressionModel(psaVeg, lcovShares, vegetation_modern_years = 1000, prefilter_taxa = T) %>% 
        suppressWarnings()
      
      save(veg_lc_model, file = glue::glue("{dir_out}/summaryResults/veg_lc_model.rda"))
      
    } else load(glue::glue("{dir_out}/summaryResults/veg_lc_model.rda"))
    
    ## summary
    # {
      # lapply(1:length(veg_lc_model), function(x) tibble(Region = names(veg_lc_model)[x]) %>% bind_cols(
      #   matrix(veg_lc_model[[x]]$sd_per_class, nrow = 1) %>%
      #     as_tibble() %>% setNames(names(veg_lc_model[[x]]$sd_per_class)))) %>%
      #   Reduce("rbind",.) %>%
      #   knitr::kable(caption = "Standard deviation (sd) per class and region") %>% suppressWarnings()
    # }
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/psaLCover.rda")) & check_files) {
      psaLCover <- make_psaLCover(ID = ID, psaVeg, lcovShares, veg_lc_model)
      save(psaLCover, file = glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
    } else load(glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
    
    # knitr::kable(psaLCover@landcover_ts[1:10,])
    
  } ## end 2.4
  
  
  ###################################
  ### RDAs                        ###
  ###################################
  
  #### 3.1 Pixel selection
  {
    if(!file.exists(glue::glue("{dir_out}/psaCrds/psaCrds.rda")) & check_files) {
      
      psa      <- psaVeg@psas %>% filter(Dataset_ID==ID) %>% st_buffer(.$'Pollen_Source_Radius [m]')
      bufferS  <- seq(0, max_buffer*1000, length = 10)
      
      nPxl <- pixelsBuffer(psa, psaLCover@resolution, class_defs, bufferS, cutoff_PSA = pxl_extract_PSA, cutoff_buffer = pxl_extract_PSA_buffer) %>% 
        group_split(lcov) %>% lapply(., function(x) {
        if(x %>% filter(radius==0) %>% pull(count) >  pxl_extract_PSA) {
          x %>% filter(radius==0) %>% mutate(cumsum = count)
        } else {
          tmp <- x %>% mutate(cumsum = cumsum(count),
                              over   = cumsum>pxl_extract_PSA_buffer)
          if(sum(tmp$over)>1) {
            tmp[1:min(which(tmp$over)),] %>% filter(count>0) %>% dplyr::select(-over)
            } else {
              tmp %>% filter(count>0)%>% dplyr::select(-over) 
            }
        }
      }) %>% Reduce("rbind",.) 
                      
      psaCrds <- lapply(unique(nPxl$lcov), function(class) {
        
        classTab <- nPxl %>% filter(lcov==class)
        
        lapply(1:nrow(classTab), function(rad) {
          
          if(classTab$radius[rad]==0) {
            bbox_ee  <- psa %>% st_geometry() %>% sf_as_ee()
          } else {
            bbox_ee <- psa %>% st_buffer(classTab$radius[rad]) %>% 
              st_difference(psa %>% st_buffer(bufferS[which(bufferS==classTab$radius[rad])-1])) %>%
              st_geometry() %>% sf_as_ee() %>% suppressWarnings()
          }
          
          dataset  <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))$clip(bbox_ee)
            
          fromList = class_defs %>% mutate(Merge_To_Class = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
              filter(Predicted_In_Past, Merge_To_Class == class) %>% pull(Class_Code) %>% as.list()
            
          toList = rep(1, length(fromList)) %>% as.list()
            
          classImage = dataset$remap(
            from = fromList,
            to = toList,
            defaultValue = -1,
            bandName = 'b1')
            
          dts <- dataset$updateMask(classImage$gt(0))
 
          table = dts$sample(
            region = bbox_ee,
            scale = (dts$projection()$nominalScale()$getInfo()),
            geometries = TRUE
          )$limit(min(c(classTab$count[rad], pxl_extract_PSA)))$getInfo()
            
          lapply(table[[4]], function(f) {
                tibble(Dataset_ID = current_psa$Dataset_ID,
                         lcov       = as.numeric(as.character(class)),
                         dist       = NA,
                         lon        = f$geometry$coordinates[1],
                         lat        = f$geometry$coordinates[2])
                }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
            mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
            suppressWarnings()
          
        }) %>% Reduce("rbind", .)
        
      }) %>% Reduce("rbind",.)  
      
      psaCrds_summary <- psaCrds %>% st_drop_geometry() %>% group_by(lcov) %>% summarise(pixel_found = n())
      
      # save pixel summary
      openxlsx::write.xlsx(psaCrds_summary, file = glue::glue("{dir_out}/summaryResults/psaCrds_summary.xlsx"), asTable = TRUE)
  
      ### save results
      save(psaCrds, file = glue::glue("{dir_out}/psaCrds/psaCrds.rda"))
        
    } 
  
  }## end 3.1
  
  #### 3.2 Add climate variables & z-transformation parameters
  {
    threshold <- 0.15 ### 15 percentile of sample size across samples in one PSA
      
    load(glue::glue("{dir_out}/psaCrds/psaCrds.rda"))
      
      if(!file.exists(glue::glue("{dir_out}/psaEnv/psaEnv.rda")) & check_files) {
        
        clim_variables <- read_xlsx(glue::glue("{data_dir}/settings/bioclim_vars_definitions.xlsx"))
        clim_dir       <- list.files(glue::glue("{data_dir}/data/global_maps"), pattern = "*.tif", full.names = T)
      
        clim_stars <- lapply(clim_variables$Var_Name, function(x) {
          ind <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
          read_stars(clim_dir[ind]) %>% setNames(x)
        })
        
        psaEnv_init <- lapply(clim_stars, function(x) st_extract(x, psaCrds) %>% pull(1) %>% unlist()) %>% Reduce("cbind", .) %>% 
          as_tibble() %>% setNames(clim_variables$Var_Name) %>% mutate(
            elev = st_extract(read_stars(clim_dir[grepl("ETOPO", clim_dir)]) %>% st_set_crs(4326) %>% setNames("elev"), psaCrds) %>% pull(1) 
          )
        
        class_z <- class_defs %>% filter(Include_Class_In_Calibration | Transfer) %>% pull(Class_Code) 
        
        ## transformation parameters
        elev   <-  read_stars(clim_dir[grepl("ETOPO", clim_dir)]) %>% st_set_crs(4326) %>%
          st_crop(psaVeg@psas %>% filter(Dataset_ID==current_psa$Dataset_ID) %>% st_buffer(.$'Pollen_Source_Radius [m]') %>% st_transform(4326)) %>%
          st_as_stars() %>% suppressMessages()
        
        allCim <- append(lapply(clim_variables$Var_Name, function(x) {
          ind  <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
          rast <- read_stars(clim_dir[ind]) %>% setNames(x) %>% 
            st_crop(psaVeg@psas %>% filter(Dataset_ID==current_psa$Dataset_ID) %>% st_buffer(.$'Pollen_Source_Radius [m]') %>% st_transform(4326)) %>%
            st_as_stars() %>% suppressMessages()
          tibble(vals = c(rast[[1]]), class = st_extract(dataset_rast, st_as_sf(st_coordinates(rast), coords = c("x", "y"), crs = 4326)) %>% pull(1)) %>%
            filter(class %in% class_z) %>% pull(vals)
          }), list(
            tibble(vals = c(elev[[1]]), class = st_extract(dataset_rast, st_as_sf(st_coordinates(elev), coords = c("x", "y"), crs = 4326)) %>% pull(1)) %>%
                     filter(class %in% class_z) %>% pull(vals)
          ))
        
        z_tranTab <- lapply(1:ncol(psaEnv_init), function(x) {
          tibble(Var_Name = c(clim_variables$Var_Name, "elev")[x], 
                 mean = mean(c(psaEnv_init %>% pull(x), allCim[[x]]), na.rm = T),
                 sd   = sd(c(psaEnv_init %>% pull(x), allCim[[x]]), na.rm = T))
        }) %>% Reduce("rbind", .)
        
        ## apply z transformation to psaEnv
        psaEnv <- lapply(1:ncol(psaEnv_init), function(x) {
          dat = psaEnv_init %>% pull(x)
          dat[is.na(dat)] <- median(dat, na.rm = T)
          tibble(c(scale(dat, z_tranTab[x,2], z_tranTab[x,3])))
        }) %>% Reduce("cbind",.) %>% setNames(names(psaEnv_init))
        
        
        save(psaEnv, file = glue::glue("{dir_out}/psaEnv/psaEnv.rda"))
        save(z_tranTab, file = glue::glue("{dir_out}/psaEnv/z_tranTab.rda"))
      }
    } ## end 3.2
    
  #### 3.3 RDA models   
  {
  
    ## psaCrds: coordinates of selected pixels per class
    ## psaEnv:  environmental variable (z-transformed)
    
    if(!file.exists(glue::glue("{dir_out}/psaOvlp/psaOvlp.rda")) & check_files) {
        
        load(glue::glue("{dir_out}/psaCrds/psaCrds.rda")) ## psaCrds
        load(glue::glue("{dir_out}/psaEnv/psaEnv.rda"))   ## psaEnv
        
        datst <- psaCrds %>% st_drop_geometry() %>% dplyr::select(-Dataset_ID) %>%
          bind_cols(psaEnv) %>% mutate(lcov = as.factor(lcov))
        
        sampleSize  <- datst %>% group_by(lcov) %>% summarise(sample = n())
        thrsh       <- quantile(sampleSize$sample, probs = threshold)
        
        dataset_sub <- sampleSize %>% mutate(sample = ifelse(sample>thrsh, thrsh, sample)) %>% 
          group_split(lcov) %>% lapply(., function(x) {
            datst %>% filter(lcov == x$lcov) %>% arrange(dist) %>% dplyr::select(-dist) %>% slice(1:x$sample)
          }) %>% Reduce("rbind", .) %>% mutate(source = "calibration")
      
        rdaOut <- rdaMod(dataset_sub %>% dplyr::select(-source), ID = ID)
        
        save(rdaOut, file = glue::glue("{dir_out}/summaryResults/rdaOut.rda"))
        
        ### Overlap and centers
        psaOvlp <- make_psaOvlp(rdaOut)
        
        # plot(psaOvlp@centers)
        # plot(psaOvlp@densities)
        
        save(psaOvlp, file = glue::glue("{dir_out}/psaOvlp/psaOvlp.rda"))
      }
  } ## end 3.3
    
  ### 3.4 plot RDA output
  {
    if(!file.exists(glue::glue("{dir_out}/summaryResults/rda_summary.png")) & check_files) {
      load(glue::glue("{dir_out}/psaOvlp/psaOvlp.rda"))
      plotRDAsummary(ID, dir_out, class_defs)
      ggsave(glue::glue("{dir_out}/summaryResults/rda_summary.png"), width = 30, height = 22, units = "cm")
    }
  }

  
  ###################################
  ### Pixel flow                  ###
  ###################################
  
  #### 4.1 Init Maps
  {

    ## load(z-trans variables)
    load(glue::glue("{dir_out}/psaEnv/z_tranTab.rda"))
      
    ## Climate maps
    clim_variables <- read_xlsx(glue::glue("{data_dir}/settings/bioclim_vars_definitions.xlsx"))
    clim_dir       <- list.files(glue::glue("{data_dir}/data/global_maps"), pattern = "*.tif", full.names = T)
      
    clim_stars <- lapply(clim_variables$Var_Name, function(x) {
      ind <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
      read_stars(clim_dir[ind]) %>% setNames(x) %>% st_crop(psaVeg@regionMap %>% st_buffer(500) %>% st_transform(4326)) %>% st_as_stars() %>% suppressMessages()
    })
      
    elev_stars <- read_stars(clim_dir[grepl("ETOPO", clim_dir)]) %>% st_set_crs(4326)  %>%
      st_warp(., clim_stars[[1]], use_gdal = F, method = "near") %>%
      setNames("elev") %>% suppressMessages() %>% suppressWarnings()
      
    env_stars <- lapply(1:nrow(z_tranTab), function(x) {
      append(clim_stars, list(elev_stars))[[x]] %>% setNames("var") %>%
        mutate(var = scale(var, z_tranTab[x,2], z_tranTab[x,3]))
    }) %>% do.call("c", .) %>% setNames(z_tranTab$Var_Name) %>% merge()
    
   ## LandCover Init Map 
   if(!file.exists(glue::glue("{dir_out}/summaryResults/map_init.png"))  & check_files) {
      dataset_rast <- read_stars(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff"))
        
      mergeTab    <- class_defs %>% mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
          dplyr::select(Class_Code, lcov, Class_Plotlabel, Predicted_In_Past, Color_Code, Transfer)
        
      initRast    <- dataset_rast %>% st_as_stars() %>% setNames("Landcover") %>% 
          mutate(Landcover = mergeTab$lcov[match(Landcover, mergeTab$Class_Code)]) %>%
          mutate(predict = mergeTab$Predicted_In_Past[match(Landcover, mergeTab$lcov)])
    
      
      ### Class transfer
      {
        trans <- class_defs %>% mutate(Merge_To_Class = ifelse(!is.na(Merge_To_Class), Merge_To_Class, Class_Code)) %>% filter(Transfer) %>% pull(Merge_To_Class) %>% unique()
        if(any(unique(c(initRast[[1]]))%in%trans)){
          
          pl_human <- ggplot() +
            geom_stars(data = initRast[1], mapping = aes(fill = as.factor(Landcover))) +
            geom_sf(psaVeg@psas %>% filter(Dataset_ID %in% ID) %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="cyan") + 
            scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                              breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                              labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
            xlab("") + ylab("") +
            theme_light()
          
          ggsave(plot = pl_human, file = glue::glue("{dir_out}/summaryResults/map_init_withHuman.png"), width = 20, height = 12, units = "cm")
          
          transPxl <- initRast[1] %>% st_as_sf(na.rm = FALSE) %>% st_centroid() %>% rownames_to_column(var = "index") %>% suppressWarnings()
          transEnv <- st_extract(env_stars, transPxl %>% filter(Landcover %in% trans)) %>% st_as_sf() %>% st_drop_geometry()  
          
          transRDA <- st_extract(
            psaOvlp@densities,
              predict(
                psaOvlp@rda_mod[[1]],
                newdata = transEnv %>% as.data.frame(),
                type = "wa",
                scaling = 1
              )[,1:2] %>% apply(., 2, zoo::na.approx) %>% as.data.frame() %>% st_as_sf(coords = c("RDA1", "RDA2"))
          ) %>% st_as_sf()
          
          toClass <- as.numeric(gsub("LC_", "", names(transRDA %>% st_drop_geometry())))
          
          transCls <- toClass[abind::abind(as.matrix(transRDA %>% st_drop_geometry()),
            lapply(1:nrow(psaOvlp@centers), function(x) (st_point(as.numeric(psaOvlp@centers[x,])) %>% st_sfc() %>% st_distance(., transRDA %>% st_geometry()))[1,]) %>%
            Reduce("cbind", .) %>% as.matrix(), along = 3) %>% apply(., 1:2, function(x) ifelse(is.na(x[1]), 0, x[1]) + (1 - x[2])) %>% apply(., 1, function(x) which.max(x))]
          
          transPxl$Landcover[as.numeric(transPxl %>% filter(Landcover %in% trans) %>% pull(index))] <- transCls
          initRast[1] <- initRast[1] %>% mutate(Landcover = transPxl %>% pull(Landcover))
        }
      }
          
      ggplot() +
        geom_stars(data = initRast[1], mapping = aes(fill = as.factor(Landcover))) +
        geom_sf(psaVeg@psas %>% filter(Dataset_ID %in% ID) %>% st_transform(4326), mapping = aes(geometry = geometry), fill = NA, color="cyan") + 
        scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                          breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                          labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
        xlab("") + ylab("") +
        theme_light()
      
      save(initRast, file = glue::glue("{dir_out}/summaryResults/initRast.rda"))
      ggsave(glue::glue("{dir_out}/summaryResults/map_init.png"), width = 20, height = 12, units = "cm")
   }
  
  }
  
  #### 4.2 Pixel Flow
  # {
  #   # if(!file.exists(glue::glue("{dir_out}/summaryResults/flowPSA.rda")) & check_files) {
  #     
  #  load(glue::glue("{dir_out}/summaryResults/initRast.rda"))
  #  load(glue::glue("{dir_out}/summaryResults/rdaOut.rda"))
  #  load(glue::glue("{dir_out}/psaOvlp/psaOvlp.rda"))
  # 
  #  psaFlow <- psaLCover@landcover_ts %>% filter(Dataset_ID==ID, Age>500) %>%
  #    group_split(Age) %>% #lapply(., function(x) x %>% dplyr::select(Dataset_ID, Age, gsub("LC_", "", psaOvlp@lcov))) %>%
  #    Reduce("rbind", .)
  #   
  # ### stars as sf & environment
  # init_sf  <- initRast[1] %>% st_as_sf(na.rm = F) %>% st_centroid() %>%
  #       rowid_to_column(var = 'pxlID') %>% relocate("pxlID", .before = "Landcover") %>%
  #       filter(Landcover %in% as.numeric(names(psaFlow)[-c(1:2)])) %>% suppressWarnings()
  # 
  # classesToMove <- init_sf %>% pull(Landcover) %>% unique() %>% sort()
  # 
  # 
  #   init_env <- st_extract(env_stars, init_sf) %>% st_as_sf() %>% st_drop_geometry() %>% 
  #         apply(., 2, function(x) { x[is.na(x)] <- median(x, na.rm = T); x }) %>% suppressMessages() 
  #   
  #   ### Landcover - Age Matrix
  #   classFlow <- matrix(nrow = nrow(init_env), ncol = length(psaFlow$Age %>% unique())+1)
  #   classFlow[,1] <- init_sf$Landcover
  #   
  #   ### classes flowPSA
  #   cls     <- st_dimensions(psaOvlp@densities)$attributes$values
  #   cls_num <- as.numeric(gsub("LC_", "", cls))
  #   
  #   ### class density and distance to centers
  #   dens  <- lapply(1:dim(psaOvlp@densities)[3], function(x) {
  #     split(psaOvlp@densities)[x,] %>% setNames("densities") %>% mutate(densities = abs(densities - 1))
  #   }) %>% do.call("c", .) %>% setNames(cls) %>% merge()
  #   
  #   dists <- lapply(1:dim(psaOvlp@densities)[3], function(x) {
  #     tmp <- split(dens)[x,]
  #     min <- (tmp %>% st_as_sf() %>% pull(1) %>% which.min())[1]
  #     tmp %>% mutate(dists = as.numeric(st_distance(tmp %>% st_as_sf(), tmp %>% st_as_sf() %>% slice(min)))*0.5) %>%
  #       dplyr::select(dists)
  #   }) %>% do.call("c", .) %>% setNames(cls) %>% merge()
  #   
  #   
  #   ### Pixel flow
  #   for(y in 2:ncol(classFlow)) {
  #     
  #     ## get RDA1 & RDA2 of all pixels
  #     rda_xy <- as.data.frame(
  #       predict(
  #         psaOvlp@rda_mod[[1]],
  #         newdata = lapply(classFlow[, y - 1], function(x)
  #           as.numeric(x == cls_num)) %>% do.call("rbind", .) %>%
  #           as.data.frame() %>% setNames(cls) %>% bind_cols(init_env),
  #         type = "wa",
  #         scaling = 1
  #       )
  #     ) %>% dplyr::select("RDA1", "RDA2") %>% st_as_sf(coords = c("RDA1", "RDA2"))
  #     
  #     ## calculate costs based on rda density & distance to center of each class
  #     costs  <-
  #       abind::abind(
  #         as.matrix(st_extract(dens, rda_xy) %>% st_as_sf() %>% st_drop_geometry()),
  #         as.matrix(st_extract(dists, rda_xy) %>% st_as_sf() %>% st_drop_geometry()),
  #         along = 3
  #       ) %>% apply(., 1:2, sum, na.rm = T)
  #     
  #     ## no costs for within class movements
  #     costs[cbind(1:nrow(costs), unlist(sapply(classFlow[,y-1], function(x) which(cls_num==x))))] <- 0
  #     
  #     newProb <- psaFlow %>% filter(Age == (psaFlow$Age %>% unique())[y-1])
  #     
  #     res <- lpSolve::lp.transport(
  #          cost.mat  = costs[1:100000,],
  #          direction = "min",
  #          row.signs = rep("==", nrow(costs[1:100000,])),
  #          row.rhs   = rep(1, nrow(costs[1:100000,])),
  #          col.signs = rep("==", ncol(costs[1:100000,])),
  #          col.rhs   = c(300, 300, 100, 100, 100, 100, 0)*100) # round(nrow(costs[1:1000,]) * unlist(c(newProb[-c(1,2)]))/100,0))
  #     
  #     
  #   }
  #   
  #   
  #     
  #   # }
  # }
  
}
  
  