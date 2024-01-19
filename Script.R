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


for(rr in 1:nrow(psa_metadata)) {
  rr <- 32
  
  cat("\n")
  print(paste0("region grid: ", psa_metadata$Dataset_ID[rr]," (",rr,"/",nrow(psa_metadata),")"))
  
  
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
  }
  
  #### 1.3 Create vegetation object for PSA
  {
    # load data, filter pollen source areas in region and/or calibration subset
    if(!file.exists(glue::glue("{dir_out}/summaryResults/psaVeg.rda")) & check_files) {
      psaVeg <- make_psaVeg(psa, radius, calibration_buffer, resolution, proj = "laea",
                            veg_folder = file.path(data_dir, "data", "vegetation_cover"), 
                            interpolate = TRUE, 
                            fc = 1/500, dt = 500, k = 5, cluster = NULL)
      save(psaVeg, file = glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
      
      plot(psaVeg)
      ggsave(glue::glue("{dir_out}/summaryResults/Map.png"), width = 10, height = 10, units = "cm")
      
    } else load(glue::glue("{dir_out}/summaryResults/psaVeg.rda"))
    
    ## Plot often not usefull since data only contains surface information for recent time slice
    # plotVegSite(psaVeg, ID = unique(psaVeg@vegetation$Dataset_ID)[6], n_taxa = 10, nrow = 2, ncol = 5)
  }
  
  
  ###################################
  ### Landcover                   ###
  ###################################
  
  #### 2.1 Modern LandCover share
  {
    ### modern time threshold for calibration
    modern_age_threshold <- 1000 
    
    modernVegetation <- psaVeg@vegetation %>% filter(Interp, Age <= modern_age_threshold) %>%
      dplyr::select(-c(Interp, Age)) %>%
      group_by(Dataset_ID) %>%
      summarise(across(everything(), function(x) mean(x, na.rm = TRUE))) %>% filter(!apply(., 1, function(x) all(is.nan(x[-1]))))
    
    class_defs <- readxl::read_excel(glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx")) 
  } ## end 2.1
  
  #### 2.2 CCI LandCover map of PSA
  {
    roi = sf_as_ee(psaVeg@regionMap %>% st_transform(4326) %>% st_shift_longitude())
    dataset <- ee$Image(glue::glue("users/slisovski/LandCoverChange/LandCover_CCI_{resolution}"))
    
    if(!file.exists(glue::glue("{dir_out}/summaryResults/dataset_rast.tiff")) & check_files) {
      dataset_rast <- getLCCfromGEE(dataset, psaVeg@regionMap, y_max = 500, buffer = 0.1)
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
  } ## end 2.2
  
  #### 2.3 Modern LandCover shares in PSA
  {
    ## remove classes for calibration that have very little pixels in study area
    class_defs[class_defs$Include_Class_In_Calibration & is.na(class_defs$Merge_To_Class) & 
                 class_defs$Class_Code %in% regLCov$Class_Code[regLCov$Percent<0.25], 'Include_Class_In_Calibration'] <- FALSE
    
    lcovShares <- modernLC(psaVeg, psaIDs = unique(modernVegetation$Dataset_ID), class_defs = class_defs, resolution = resolution)
    
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
      psaLCover <- make_psaLCover(ID = current_psa$Dataset_ID, psaVeg, lcovShares, veg_lc_model)
      save(psaLCover, file = glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
    } else load(glue::glue("{dir_out}/summaryResults/psaLCover.rda"))
    
    # knitr::kable(psaLCover@landcover_ts[1:10,])
    
  } ## end 2.4
  
  
  ###################################
  ### RDAs                        ###
  ###################################
  
  #### 3.1 Pixel selection
  {
    pxl_extract_PSA        <- 800       # max. points for RDA to search within the PSA   
    pxl_extract_PSA_buffer <- 50        # min. points for RDA to search in buffer rings if they are not found in the PSA 
    
    max_buffer      <- 1000           # maximum buffer for search
    method          <- "getInfo"     # Options "ee_as_sf", "getInfo"
    
    if(!file_exists(glue::glue("{dir_out}/psaCrds/psaCrds_{current_psa$Dataset_ID}.rda")) & check_files) {
      
      psa      <- psaVeg@psas %>% filter(Dataset_ID==current_psa$Dataset_ID) %>% st_buffer(.$'Pollen_Source_Radius [m]')
      bufferS  <- seq(0, max_buffer*1000, length = 10)
      
      nPxl    <- pixelsBuffer(psa, psaLCover@resolution, class_defs, bufferS, cutoff_PSA = pxl_extract_PSA, cutoff_buffer = pxl_extract_PSA_buffer) %>%
                  group_split(lcov) %>% lapply(., function(x) {
                        x %>% mutate(count = x$count - c(0, x$count[-nrow(x)]),
                                     cumsum = cumsum(count),
                                     aim = c(pxl_extract_PSA, rep(pxl_extract_PSA_buffer, nrow(.)-1)),
                                     filter = cumsum >= aim) %>% slice(1:min(which(filter))) %>%
                              rowwise() %>% mutate(extrN = min(cumsum, aim)) %>% ungroup()
                      }) %>% Reduce("rbind", .)
                      
      psaCrds <- lapply(unique(nPxl$lcov), function(class) {
        
        classTab <- nPxl %>% filter(lcov==class)
        
        lapply(1:nrow(classTab), function(rad) {
          
          if(classTab$count[rad]>0) {
            
            if(rad==1) {
              bbox_ee  <- psa %>% sf_as_ee()
            } else {
              bbox_ee <- psa %>% st_buffer(classTab$radius[rad]) %>% st_difference(psa %>% st_buffer(classTab$radius[rad-1])) %>% st_geometry() %>% sf_as_ee()
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
            
            # Map$addLayer(dts) +
            #   Map$addLayer(bbox_ee)
                  
            table = dts$sample(
              region = bbox_ee,
              scale = (dts$projection()$nominalScale()$getInfo()),
              geometries = TRUE
            )$limit(classTab$extrN[rad])$getInfo()
            
            lapply(table[[4]], function(f) {
                  tibble(Dataset_ID = current_psa$Dataset_ID,
                           lcov       = as.numeric(as.character(class)),
                           dist       = NA,
                           lon        = f$geometry$coordinates[1],
                           lat        = f$geometry$coordinates[2])
                  }) %>% Reduce("rbind",.) %>% st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
              mutate(dist = (st_distance(., psa %>% st_centroid() %>% st_geometry() %>% st_transform(4326))/1000)[,1]) %>%
              suppressWarnings()
            
          } else NULL
          
        }) %>% Reduce("rbind", .)
        
      }) %>% Reduce("rbind",.)  
      
      psaCrds_summary <- psaCrds %>% st_drop_geometry() %>% group_by(lcov) %>% summarise(pixel_found = n())
      
      # save pixel summary
      openxlsx::write.xlsx(psaCrds_summary, file = glue::glue("{dir_out}/summaryResults/psaCrds_summary_{ID}.xlsx"), asTable = TRUE)
      
      ### save results
      save(psaCrds, file = glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda"))
        
    } 
  
  }## end 3.1
  
  
  #### 3.2 Add climate variables & z-transformation parameters
  {
      load(glue::glue("{dir_out}/psaCrds/psaCrds_{ID}.rda"))
      
      if(!file.exists(glue::glue("{dir_out}/psaEnv/psaEnv_{ID}.rda")) & check_files) {
        
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
        
        ##### NUR KLASSEN DIE AUCH BEWEGT WERDEN SOLLEN!!!!
        
        ## transformation parameters
        allCim <- append(lapply(clim_variables$Var_Name, function(x) {
          ind  <- which(sapply(clim_dir, function(y) grepl(glue::glue("_{x}_"), y, fixed = T)))
          c((read_stars(clim_dir[ind]) %>% setNames(x) %>% 
            st_crop(psaVeg@psas %>% filter(Dataset_ID==current_psa$Dataset_ID) %>% st_buffer(.$'Pollen_Source_Radius [m]') %>% st_transform(4326)) %>%
            st_as_stars())[[1]]) %>% suppressMessages()
          }), list((read_stars(clim_dir[grepl("ETOPO", clim_dir)]) %>% st_set_crs(4326) %>%
                 st_crop(psaVeg@psas %>% filter(Dataset_ID==current_psa$Dataset_ID) %>% st_buffer(.$'Pollen_Source_Radius [m]') %>% st_transform(4326)) %>%
                 st_as_stars())[[1]]) %>% suppressMessages())  
        
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
        
        
        save(psaEnv, file = glue::glue("{dir_out}/psaEnv/psaEnv_{ID}.rda"))
        save(psaEnv, file = glue::glue("{dir_out}/psaEnv/z_tranTab_{ID}.rda"))
      }
    } ## end 3.2
    
  #### 3.3   
  
}
  
  