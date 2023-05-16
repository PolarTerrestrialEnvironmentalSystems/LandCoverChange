#' Simulation of Landcover change over time, based on sedimentary pollen record time series
#'
#' ....
#'
#' @name LandCoverChange-package
#' @docType package
#' @author Peter Ewald and Simeon Lisovski.
NULL

## Data import
##'
##' The function \code{psaMap} loads all vegetation cover files from the specified folder, builds pollen source areas (PSA), i.e. circular polygons around the site's centres with its pollen source radius. 
##' The vegetation cover is filtered to all PSAs within the two maps provided (for computing past land cover and for the calibration of the model for vegetation cover to land cover translation.)
##' @title psaMap
##' @param veg_folder character, absolute path to the folder where vegetation cover files were copied to
##' @param regMap sfc_MULTIPOLYGON, map for that the land cover change shall be computed
##' @param calibMap sfc_MULTIPOLYGON, map from which sites are selected to build the translation model from vegetation cover to landcover time series
##' @param save_folder character, absolute path to the directory where the results list and if plt == T the plots are saved to
##' @param plt boolean, whether to plot the pollen source areas, the number of samples per site, and the oldest age of the samples per site within the maps passed
##' @return a list object
##' \item{\code{data.frame}}{veg_cover: vegetation cover in % per taxon with Dataset_ID, Age, Intersects_regMap, Intersects_calibMap columns}
##' \item{\code{...}}{...: ...}
##' @importFrom data.table fread
##' @importFrom dplyr select
##' @importFrom dplyr full_join
##' @importFrom dplyr pull
##' @importFrom dplyr group_by
##' @importFrom sf st_as_sf
##' @importFrom sf st_intersects
##' @importFrom sf st_union
##' @importFrom sf st_crs
##' @importFrom sf st_set_crs
##' @importFrom sf st_transform
##' @importFrom sf st_make_valid
##' @importFrom sf st_buffer
##' @export
psaMap <- function(veg_folder, regMap, calibMap, save_folder = NULL, plt = T) {
  
  # packages needed: data.table, sf, dplyr,rnaturalearth
  
  sf_use_s2(F)
  return_list = list()
  
  # load vegetation and meta data
  {
    cat("\nLoad meta and vegetation data\n")
    # read vegetation cover and meta data
    file_names <- list.files(veg_folder)
    veg_cover <- lapply(file_names, 
                        function(x) {
                          veg_cov = data.table::fread(file.path(veg_folder, x))
                        })
    names(veg_cover) <- file_names
    
    # extract meta data
    meta_cols_veg = c("Dataset_ID", "Longitude", "Latitude", "Pollen_Source_Radius [m]")
    meta_data = unique(do.call(rbind,
                               lapply(1:length(veg_cover),
                                      function(x){
                                        veg_cover[[x]] %>% dplyr::select(all_of(meta_cols_veg)) %>% 
                                          cbind(File = names(veg_cover)[x])
                                      })))
    centres = meta_data %>% sf::st_as_sf(coords = c("Longitude", "Latitude")) %>% sf::st_set_crs(4326)

    # extract vegetation cover (Dataset_ID, Age, and taxa percentages)
    veg_cover = lapply(1:length(veg_cover),
                       function(x){
                         veg_cov = veg_cover[[x]] %>%
                           dplyr::select(all_of(c("Dataset_ID", 
                                           grep("mean", colnames(veg_cover[[x]]), value = T))))
                         colnames(veg_cov) = gsub(" [mean of cover in %]", "", 
                                                  colnames(veg_cov), fixed = T)
                         colnames(veg_cov) = gsub("_mean [yrs BP]", "", 
                                                  colnames(veg_cov), fixed = T)
                         return(veg_cov)
                       })
    for (ii in 1:length(veg_cover)){
      if (ii == 1){
        veg_cover_full = veg_cover[[ii]]
      } else{
        veg_cover_full = suppressMessages(
          veg_cover_full %>% dplyr::full_join(veg_cover[[ii]])
        )
      }
    } 
    veg_cover = veg_cover_full  %>% as.data.frame(); rm(veg_cover_full) 
    
    return_list$veg_cover = veg_cover
    return_list[["meta_data"]] = list()
    return_list$meta_data[["centres"]] = centres
  }
  
  # create map with pollen source areas based on provided maps
  {
    cat("\nCreate map with pollen source areas")
    
    # project sites to target geometry from maps
    if (!(st_crs(regMap) == st_crs(calibMap))){
      stop("The same projection for the region and calibration map must be provided.")
    }
    centres = centres %>% st_transform(st_crs(regMap))
    
    # filter sites in any regMap or calibrationMap
    all_maps = st_union(st_union(regMap, calibMap$geometry))
    intersects = st_intersects(centres$geometry, all_maps, sparse = F)
    centres = centres[intersects,]
    
    # build psas, enlarge maps (such that all psas fit in; seperately for calibMap and regMap)
    psas = centres %>% st_buffer(.$`Pollen_Source_Radius [m]`)  
    all_calibMaps = st_union(calibMap)
    worldMap  <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
      st_shift_longitude() %>% st_transform(st_crs(regMap)) %>% st_make_valid()
    intersects = data.frame(Dataset_ID = psas$Dataset_ID,
                            Intersects_regMap = st_intersects(psas$geometry, regMap, sparse = F),
                            Intersects_calibMap = st_intersects(psas$geometry, all_calibMaps, sparse = F))
    psas = psas %>% left_join(intersects, by = "Dataset_ID")
    centres = centres %>% left_join(intersects, by = "Dataset_ID")
    regMap = st_bbox(regMap %>% st_union(psas %>% filter(Intersects_regMap) %>% pull(geometry))) %>% 
      st_as_sfc() %>% st_set_crs(st_crs(regMap))
    regMap = suppressWarnings(
      worldMap %>% st_intersection(regMap) %>% dplyr::select(geometry) %>% st_buffer(50) %>% st_union()
    )
    calibMap = suppressWarnings(
      st_bbox(all_calibMaps %>% st_union(psas %>% filter(Intersects_calibMap) %>% pull(geometry))) %>% 
      st_as_sfc() %>% st_set_crs(st_crs(calibMap))
    )
    calibMap = suppressWarnings(
      worldMap %>% st_intersection(calibMap) %>%  
        dplyr::select("continent") %>% rename(groups = continent) %>% group_by(groups) %>% 
        summarize(geometry = st_union(st_buffer(geometry, 20)))
    )
    
    # return enlarged maps
    return_list$meta_data[["regMap"]] = regMap
    return_list$meta_data[["calibMap"]] = calibMap
    return_list$meta_data[["psas"]] = psas
    return_list$meta_data[["centres"]] = centres
  }
  
  # filter vegetation data based on selected pollen source areas, drop taxa not present, and normalize rows to 1
  {
    cat("\nFilter vegetation data within maps")
    
    # filter sites if in any map
    veg_cover = veg_cover %>% 
      left_join(intersects, by = "Dataset_ID")
    meta_cols_veg = c("Dataset_ID", "Age", c("Intersects_regMap", "Intersects_calibMap"))
    veg_cover = veg_cover %>% filter(Intersects_regMap | Intersects_calibMap) 
    
    # remove samples with no age and replace any NA percentage with 0
    veg_cover = veg_cover %>% filter(!is.na(Age)) %>% 
      mutate(across(where(is.numeric), ~replace_na(.x, 0)))
    
    # reorder columns
    veg_cover = veg_cover[,c(meta_cols_veg, 
                             colnames(veg_cover)[!(colnames(veg_cover) %in% meta_cols_veg)])]
    
    # drop taxa not present at any site
    taxa_cols_veg = colnames(veg_cover)[!(colnames(veg_cover) %in% meta_cols_veg)]
    taxa_cols_veg = names(which(!apply(veg_cover[,taxa_cols_veg], 2, function(u) all(u == 0))))
    veg_cover = veg_cover[, c(meta_cols_veg, taxa_cols_veg)]
    
    # normalize all taxa to 1
    veg_cover = data.frame(veg_cover[,meta_cols_veg],
                            t(apply(veg_cover[,taxa_cols_veg], 1, function(u){ u/sum(u) }))) 
    
    # return vegetation cover
    return_list$veg_cover = veg_cover
  }
  
  # plot sites and vegetation oldest /nr. of vegetation samples
  if (plt){
    
    # plot all PSAs on regionMap and calibMap
    {
      cat("\nPlot Pollen Source areas within the region map")
      if (!is.null(save_folder)){
        png(file.path(save_folder, paste0("pollen_source_areas_regMap.png")),
            width = 1200, height = 800)
      }
      plot(regMap, border = "grey10", lwd = 0.5, 
           main = paste0(nrow(psas[psas$Intersects_regMap,]), " Pollen Source Areas and site centres within region map"))
      plot(psas$geometry[psas$Intersects_regMap], col = "darkred", border = NA, add = T)
      plot(centres$geometry[psas$Intersects_regMap], pch = 3, size = .5, col = "white", add = T)
      if (!is.null(save_folder)){ dev.off() }
      
      cat("\nPlot Pollen Source areas within the calibration map")
      if (!is.null(save_folder)){
        png(file.path(save_folder, paste0("pollen_source_areas_calibMap.png")),
            width = 1200, height = 800)
      }
      plot(calibMap$geometry, border = "grey10", lwd = 0.5, 
           main = paste0(nrow(psas[psas$Intersects_calibMap,]), " Pollen Source Areas and site centres within calibration map"))
      plot(psas$geometry[psas$Intersects_calibMap], add = T, col = "darkred", border = NA)
      plot(centres$geometry[psas$Intersects_calibMap], pch = 3, size = .5, col = "white", add = T)
      if (!is.null(save_folder)){ dev.off() }
    }
    
    # plot oldest samples and number of samples per site
    {
      cat("\nPlot number of and age of oldest pollen samples per site")
      # prepare plot data frame
      plt_df = centres %>% 
        left_join(veg_cover %>% group_by(Dataset_ID) %>% summarize(Oldest_Sample = max(Age)), by = "Dataset_ID") %>%
        left_join(veg_cover %>% group_by(Dataset_ID) %>% summarize(Nr_Samples = length(Age)), by = "Dataset_ID") %>% 
        mutate(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2])
      
      # plot number of samples within calibration map
      plt = ggplot() +
        geom_sf(data = calibMap$geometry, fill='darkgreen', alpha = 0.25) +
        geom_point(data = plt_df %>% filter(Intersects_calibMap), 
                   aes(x = x, y = y, 
                       colour = Nr_Samples, 
                       size = Nr_Samples)) + 
        scale_colour_distiller(palette = "YlOrBr", direction = 1) +
        labs(x = "Longitude", y = "Latitude", title=paste0("Number of pollen samples per site within calibration map"),
             colour = "Number of pollen samples",
             size = "Number of pollen samples") +
        theme_minimal()
      if (!is.null(save_folder)){
        png(file.path(save_folder, paste0("number_pollen_samples_calibMap.png")),
            width = 1200, height = 800)
        print(plt)
        dev.off()
      } else{
        print(plt)
      }
      
      # plot number of samples within region map
      plt = ggplot() +
        geom_sf(data = regMap, fill='darkgreen', alpha = 0.25) +
        geom_point(data = plt_df %>% filter(Intersects_regMap), 
                   aes(x = x, y = y, 
                       colour = Nr_Samples, 
                       size = Nr_Samples)) + 
        scale_colour_distiller(palette = "YlOrBr", direction = 1) +
        labs(x = "Longitude", y = "Latitude", title=paste0("Number of pollen samples per site within region map"),
             colour = "Number of pollen samples",
             size = "Number of pollen samples") +
        theme_minimal()
      if (!is.null(save_folder)){
        png(file.path(save_folder, paste0("number_pollen_samples_regMap.png")),
            width = 1200, height = 800)
        print(plt)
        dev.off()
      } else{
        print(plt)
      }
     
      # plot oldest age of samples within calibration map
      plt = ggplot() +
        geom_sf(data = calibMap$geometry, fill='darkgreen', alpha = 0.25) +
        geom_point(data = plt_df %>% filter(Intersects_calibMap), 
                   aes(x = x, y = y, 
                       colour = Oldest_Sample, 
                       size = Oldest_Sample)) + 
        scale_colour_distiller(palette = "YlOrBr", direction = 1) +
        labs(x = "Longitude", y = "Latitude", title=paste0("Age of oldest pollen sample per site within calibration map"),
             colour = "Age of oldest pollen sample",
             size = "Age of oldest pollen sample") +
        theme_minimal()
      if (!is.null(save_folder)){
        png(file.path(save_folder, paste0("oldest_pollen_samples_calibMap.png")),
            width = 1200, height = 800)
        print(plt)
        dev.off()
      } else{
        print(plt)
      }
      
      # plot oldest age of samples within region map
      plt = ggplot() +
        geom_sf(data = regMap, fill='darkgreen', alpha = 0.25) +
        geom_point(data = plt_df %>% filter(Intersects_regMap), 
                   aes(x = x, y = y, 
                       colour = Oldest_Sample, 
                       size = Oldest_Sample)) + 
        scale_colour_distiller(palette = "YlOrBr", direction = 1) +
        labs(x = "Longitude", y = "Latitude", title=paste0("Age of oldest pollen sample per site within region map"),
             colour = "Age of oldest pollen sample",
             size = "Age of oldest pollen sample") +
        theme_minimal()
      if (!is.null(save_folder)){
        png(file.path(save_folder, paste0("oldest_pollen_samples_regMap.png")),
            width = 1200, height = 800)
        print(plt)
        dev.off()
      } else{
        print(plt)
      }
    }
  }
  
  # save or return results
  if (is.null(save_folder)){ 
    return(return_list) 
  } else{
    fl_name =  file.path(save_folder, "psaMap.rda")
    cat("\nSave to file", fl_name)
    save(return_list, file = fl_name)
  }
}

## interpolate vegetation cover
##'
##' ...
##' @title modernShares
##' @param veg_cover
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
modernVegShares = function(veg_cover, 
                           time_span_modern,
                           meta_cols_veg =  c("Dataset_ID", "Age", 
                                              "Intersects_regMap", "Intersects_calibMap"),
                           save_folder = NULL){
  
  taxa_cols_veg = colnames(veg_cover %>% select(-all_of(meta_cols_veg)))
  veg_modern = veg_cover %>% filter(Intersects_calibMap) %>% filter(Age <= time_span_modern)
  veg_modern = as.data.frame(veg_modern[,c("Dataset_ID", taxa_cols_veg)] %>%
                               group_by(Dataset_ID) %>%
                               summarise(across(everything(), mean, na.rm=T)))
  row.names(veg_modern) = NULL
  
  if (!is.null(save_folder)){
    save(veg_modern, file = file.path(save_folder, "modernVegShares.rda"))
  } else{
    return(veg_modern)
  }
}

## interpolate vegetation cover of one site
##'
##' ...
##' @title modernShares
##' @param veg_cover
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
interpolate_ages_per_site = function(df, 
                                     #ID, 
                                     meta_cols,
                                     fc = 1/500, dt = 500,
                                     k = 5, int.method = "linear",
                                     appliedFilter = "gauss"){
  library(corit)
  
  taxa_cols = colnames(df)[!(colnames(df) %in% meta_cols)] 
  df_taxa = df[, taxa_cols]
  
  # sort out all NA columns if existing
  df_taxa = df_taxa %>% select(which(colSums(is.na(.)) == 0))
  taxa_cols = colnames(df_taxa)[!(colnames(df_taxa) %in% meta_cols)] 
  
  # in case all sample are within one time intervall average the vegetation cover to age dt
  # to keep 1 time slice at least. from interpolation there would results 0 time slices.
  if (df$Age[nrow(df)] <= dt){
    df = data.frame(Dataset_ID = unique(df$Dataset_ID), 
                    Age = dt, 
                    t(apply(df %>% select(-Dataset_ID, -Age), 2, mean)),
                    stringsAsFactors = F)
  } else{
    
    # interpolate to 500 years grid
    df_taxa_interp = do.call(cbind, 
                             lapply(
                               taxa_cols,
                               function(u){
                                 col = zoo(df_taxa[,u], order.by = df$Age)
                                 col_interp = InterpolationMethod(X = col, fc=fc, dt=dt,
                                                                  timser.length = tail(df$Age, n=1), 
                                                                  int.method = int.method, 
                                                                  appliedFilter = appliedFilter, 
                                                                  k=k)
                                 ages = index(col_interp)
                                 df_ret = data.frame(Age = ages, col_interp)
                                 colnames(df_ret) = c("Age", u); row.names(df_ret) = NULL
                                 
                                 return(df_ret)
                                 
                               }))
    df_taxa_interp$Dataset_ID = unique(df$Dataset_ID)
    df = df_taxa_interp[, c("Dataset_ID", "Age", taxa_cols)]
  }
  
  
  return(df)
}

## interpolate vegetation cover
##'
##' can have one or multiple IDs in columns Dataset_ID. If multiple, Interpolation is done seperately and the result is joined.
##' @title modernShares
##' @param veg_cover
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
 

interpVegTS = function(veg_cover, save_folder = NULL){
  
  # packages needed: parallel, zoo, dplyr
  
  cat("\n\nInterpolate past vegetation time series")
  IDs = unique(veg_cover$Dataset_ID)
  
  veg_cover_interp = parallel::mclapply(IDs,
                                        function(x){
                                          veg_site = veg_cover %>% filter(Dataset_ID == x) %>% arrange(Age)
                                          veg_site_interp = suppressWarnings(
                                            interpolate_ages_per_site(df = veg_site %>%
                                                                        select(-Intersects_regMap, 
                                                                               -Intersects_calibMap),
                                                                      meta_cols = c("Dataset_ID",
                                                                                    "Age"),
                                                                      fc = fc, dt = dt, k = k)
                                          )
                                          
                                          veg_site_interp = suppressMessages(
                                            veg_site_interp %>% 
                                              left_join(unique(veg_site %>% 
                                                                 select(all_of(c("Dataset_ID", "Intersects_regMap", "Intersects_calibMap")))))
                                          )
                                          
                                          return(veg_site_interp)
                                        }, mc.cores = parallel_max)
  veg_cover_interp = suppressMessages(Reduce(full_join, veg_cover_interp))
  
  if (!is.null(save_folder)){
    save(veg_cover_interp, file = file.path(main_dir, "interim_results", "interpVegTS", "interpVegTS.rda"))
  } else{
    return(veg_cover_interp)
  }
}


# ////////////////////////////////
# helper functions (not exported)

create_subfolders = function(x, folder_path, indent=1) {
  
  if (!file.exists(folder_path)){ dir.create(folder_path) }
  
  if ((is.list(x)) & (length(x) >= 1)){
    for (ii in 1:(length(x))){
      folders = x[[ii]]
      if (!all(sapply(folders, class) == "list")){
        stop("Subfolders are only accepted as named lists.")
      }
      folder_name = names(x)[ii]
      stop_at_last = F
      if (is.null(folder_name)){
        folder_name = names(folders)[1]
        stop_at_last = T
      } 
      folder_name = gsub(" ", "", folder_name)
      
      cat("\n",paste0(rep("\t", indent, collapse=T)), folder_name)
      final_path = file.path(folder_path, gsub(" ", "", folder_name))
      if (!file.exists(final_path)){
        dir.create(final_path)
      }
      
      if (length(folders) >= 1){
        if (!stop_at_last){
          create_subfolders(x = folders, folder_path = file.path(folder_path, folder_name), 
                            indent = indent+1)
        }
      } 
    }
    
    if (indent == 1) { cat("\n\n") }
  }
}

load_rda = function(file_name){
  # loads an rda file, and returns it
  load(file_name)
  get(ls()[ls() != "file_name"]) 
}


