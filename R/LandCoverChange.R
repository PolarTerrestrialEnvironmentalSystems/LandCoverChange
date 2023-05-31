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
    
    # drop taxa not present at any site and if existing an "Indeterminable" column (taxa in pollen cores not determined)
    taxa_cols_veg = colnames(veg_cover)[!(colnames(veg_cover) %in% meta_cols_veg)]
    taxa_cols_veg = names(which(!apply(veg_cover[,taxa_cols_veg], 2, function(u) all(u == 0))))
    taxa_cols_veg = taxa_cols_veg[taxa_cols_veg != "Indeterminable"]
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
  veg_modern = veg_cover %>% 
    filter(Age <= time_span_modern)
  veg_modern = as.data.frame(veg_modern[,c("Dataset_ID", taxa_cols_veg)] %>%
                               group_by(Dataset_ID) %>%
                               summarise(across(everything(), mean, na.rm=T))) %>% 
    left_join(veg_cover %>% select(all_of(meta_cols_veg)) %>% select(-Age), by = "Dataset_ID")
  row.names(veg_modern) = NULL
  veg_modern = veg_modern[, c(meta_cols_veg[meta_cols_veg != "Age"], colnames(veg_modern %>% select(-any_of(meta_cols_veg))))]
  
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

## plot time series of interpolated and original vegetation cover
##'
##' ...
##' @title plotVegSite
##' @param veg_cover
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
plotVegSite = function(vegCover1, name1, 
                  vegCover2 = NULL,  name2 = NULL, 
                  ID = NULL, meta_cols = c("Dataset_ID", "Age",
                                           "Intersects_regMap", "Intersects_calibMap"),
                  save_folder = NULL,
                  n_taxa = NULL){
  if (is.null(ID)){
    df_plot = suppressMessages(
      full_join(cbind(vegCover1, Label = name1),
                cbind(vegCover2, Label = name2))
    )
    
    ID = paste0(unique(data$Dataset_ID), collapse = "_")
  } else{
    df_plot = suppressMessages(
      full_join(cbind(vegCover1 %>% filter(Dataset_ID == ID), 
                      Label = name1),
                cbind(vegCover2 %>% filter(Dataset_ID == ID), 
                      Label = name2))
    )
  }
  
  cols = colnames(df_plot %>% select(-any_of(c(meta_cols, "Label"))))
  
  # select n_taxa taxa
  if (is.null(n_taxa)){ n_taxa = length(cols) }
  cols = names(sort(apply(df_plot[,cols], 2, mean, na.rm = T), decreasing = T))[1:n_taxa]
  df_plot = df_plot[, c(meta_cols, "Label", cols)]
  
  df_plot$Label = factor(df_plot$Label, levels = c(name1, name2))
  plt = suppressMessages(
    plotTS(data = df_plot, 
                     cols = cols,
                     age_name = "Age",
                     round_ages = F, age_steps = dt, lines_var = "Label",
                     color = "indianred", alpha = 0.7,
                     title = "",
                     share_to_percent = F,
                     subtitle = paste0("fc = ", fc, ", dt = ", dt,
                                       " k = ", k))
  )
  
  if (nrow(df_plot) == 1){
    plt = plt + 
      labs(subtitle = paste0("Core has only samples within present and ", dt, ".",
                             " Those have been averaged and set to year ", 
                             dt, ". These are currently not visible in this plot,",
                             "actually there would be points."))
  }
  
  # save or return plt for further modification
  if (!is.null(save_folder)){
    clear()
    png(filename = file.path(save_folder, paste0("site_", ID, )), 
        width=1200, height=800)
    print(plt)
    dev.off()
  } else{
    return(plt)
  }
}

## plot time series of interpolated and original vegetation cover
##'
##' ...
##' @title plotVegSite
##' @param veg_cover
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
plotTS = function(data, cols, age_name, round_ages = T, age_steps = 1000,
                            color, lines_var = NULL, lwd = 1, alpha = 1, share_to_percent = F,
                            title = "", subtitle = ""){
  
  library(ggplot2); library(dplyr); library(tidyr); library(scrutiny)    
  
  # prepare data
  {
    if (!is.null(lines_var)){
      data = data[, c(age_name, cols, lines_var)]
    } else{
      data = data[, c(age_name, cols)]
    }
    
    data = data %>% pivot_longer(cols = all_of(cols), 
                                 names_to = "Plant",
                                 values_to = "Percentage")
    if (share_to_percent){ data$Percentage = data$Percentage * 100 }
    
    if (!is.null(lines_var)){
      if (!is.factor(data[[lines_var]])){
        data[[lines_var]] = as.factor(data[[lines_var]])
      }
    }
  }
  
  # round axis limits
  {
    if (round_ages){
      scale_max_x = round_ceiling(max(data[[age_name]], na.rm = T), -3)
      scale_min_x = round_floor(min(data[[age_name]], na.rm = T), -3)
      step_size_x = age_steps
      values_x = seq(scale_min_x, scale_max_x, step_size_x)
    } else {
      values_x = unique(data[[age_name]])
    }
    
    scale_max_y = round_ceiling(max(data$Percentage, na.rm = T), -1)
    scale_min_y = round_floor(min(data$Percentage, na.rm = T), -1)
    step_size_y = 10
    values_y = seq(scale_min_y, scale_max_y, step_size_y)
  }
  
  plt = ggplot(data, aes(x = !!as.name(age_name), y = Percentage)) + 
    list(
      if (!is.null(lines_var)){
        geom_line(aes(colour = !!as.name(lines_var)), linewidth = lwd,
                  alpha = alpha)
      } else{
        geom_line(colour = color, linewidth = lwd,
                  alpha = alpha)
      }
    ) + 
    facet_wrap(~Plant, ncol = length(cols),
               scales = "fixed") +
    scale_x_reverse(breaks = values_x) +
    scale_y_continuous(breaks = values_y) +
    coord_flip() + 
    theme_minimal() + 
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.grid.major.x = element_line(colour="grey80", 
                                        linewidth=0.075,
                                        linetype="solid"), 
      panel.grid.major.y = element_line(colour="grey80", 
                                        linewidth=0.075,
                                        linetype="solid"),
      panel.background = element_rect(color=NA, fill="grey97"), 
      axis.line = element_line(colour = "grey20"), 
      strip.text.x = element_text(face = "bold", 
                                  size = 12, 
                                  angle = 0, 
                                  vjust = 0.5, 
                                  hjust = 0.5), 
      strip.background = element_rect(color=NA, fill="grey97"),
      strip.text.y = element_text(angle = 0), 
      panel.border = element_blank(), 
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      plot.title = element_text(size = 16,
                                hjust = 0),
      plot.subtitle = element_text(size = 14),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.margin = unit(rep(1.25, 4), "cm")) +
    labs(x = "Percentage (%)", y = "Age BP (years)",
         title = title, subtitle = subtitle) 
  
  return(plt)
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

## extract modern landcover shares with google earth engine
##'
##' 
##' @title modernLCShares
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
modernLCShares = function(polygons, resolution, class_defs, merge_classes = NULL,
                          save_folder = NULL){
  
  
  # extract shares within polygons with GEE
  {
    poly_ee = sf_as_ee(polygons$geometry)
    
    #dataset = ee$Image('ESA/GLOBCOVER_L4_200901_200912_V2_3')$select('landcover')
    dataset = ee$Image("users/slisovski/LandCoverChange_LCCS")
    # Map$addLayer(dataset, {}, 'Landcover')
    # ee_print(dataset)
    
    # extract number of pixel per class
    class_areas  = ee$Image$pixelArea()$addBands(dataset)$reduceRegions(
      collection = poly_ee,
      reducer    = ee$Reducer$count()$group(groupField = 1),
      scale      = resolution
    )$getInfo()
    
    class_defs_orig = class_defs %>% filter(is.na(Merge_To_Class)) %>% select(Class_Code, Class_Plotlabel) %>% 
      rename(lcov = Class_Code)
    class_defs_orig$lcov = factor(class_defs_orig$lcov)
    
    regLCov = lapply(class_areas$features, function(x) {
      lapply(x[[4]][[1]], function(y) tibble(group = x[[3]], 
                                             lcov = as.factor(as.numeric(y$group)), 
                                             count = y$count
      )) %>% bind_rows() 
    }) %>% do.call("rbind",.) 
    regLCov = regLCov %>% left_join(data.frame(group = as.character(0:(nrow(polygons)-1)),
                                               Dataset_ID = polygons$Dataset_ID))
    
    regLCov = regLCov %>% 
      full_join(tibble(group = as.character(rep(seq(0, length(polygons)-1), each = nrow(class_defs_orig))), 
                       lcov = rep(class_defs_orig$lcov, length(polygons))), 
                by = c("group", "lcov")) %>%
      arrange(as.numeric(as.character(group)), as.numeric(as.character(lcov))) %>% 
      filter(!is.na(Dataset_ID)) %>%
      mutate(count = ifelse(is.na(count), 0, count)) %>%
      group_by(group) %>% mutate(perc = (count/sum(count, na.rm = T))*100) %>% dplyr::select(-count) %>%
      pivot_wider(names_from = lcov, values_from = perc, values_fill = 0) %>% ungroup() %>% select(-group)
    
    
  }
  
  # merge classes and sort columns
  {
    if (!is.null(merge_classes)){
      for (merge_class in merge_classes){
        old_cl = merge_class$old_classes
        new_cl = merge_class$new_class
        
        regLCov = regLCov %>% 
          mutate(new_cl = apply(regLCov %>% select(any_of(as.character(old_cl))), 1, sum)) %>%
          select(-any_of(as.character(old_cl)))
        regLCov[[as.character(new_cl)]] = regLCov[["new_cl"]]
        regLCov = regLCov %>% select(-new_cl)
      }
    }
    
    regLCov = regLCov[, c("Dataset_ID", as.character(sort(as.numeric(colnames(regLCov %>% select(-Dataset_ID))))))]
  }
  
  if (!is.null(save_folder)){
    save(regLCov, file = file.path(save_folder, "modernLCShares.rda"))
  } else{
    return(regLCov)
  }
}

## ...
##'
##' 
##' @title k_fold_cv
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
k_fold_cv = function(..., data, k, x_names, y_names, analysis_function, 
                     set_seed = F, seed_int = NULL, prnt = F, 
                     tune = F, average="num", verb = 1){
  
  dots = list(...)
  
  if (k > nrow(data)){
    stop("k-fold CV: k must be smaller than the number of rows in data.")
  }
  # build folds
  {
    cv_results = list()
    folds = list() 
    fold.size = nrow(data)/k
    remain = 1:nrow(data)
    for (ii in 1:k){
      
      if(set_seed){
        if (!is.null(seed_int)){
          set.seed(seed_int)
        } else {
          stop("k-fold CV: provide seed_int when set_seed = T")
        }
      }
      
      subsample = sample(remain, fold.size, replace = F)
      
      folds[[ii]] = subsample # store indices
      if (ii == k){
        folds[[ii]] = remain
      }
      
      remain = setdiff(remain, subsample)
    }
    cv_results[["folds"]] = list()
  }
  
  for (ii in 1:k){
    if (prnt){ cat("\n\nFold:", ii) }
    inds = folds[[ii]] 
    train = data[-inds, ] 
    test = data[inds, ]
    
    # analysis is a function that fits, tunes, whatever a model on training set,
    # evaluates it on the test set and computes a metric for the current fold
    # it must return at least the model and the metric
    analysis_results = analysis_function(..., df_train = train, df_test = test, 
                                         x_names = x_names, y_names = y_names,
                                         prnt = prnt, tune = tune, set_seed = set_seed,
                                         verb = verb) 
    if (prnt){ cat("\n\n\tMean absolute error on training set of fold", 
                   ii, ":", analysis_results$metric) }
    
    cv_results$folds[[ii]] = analysis_results
  }
  
  # extract average best parameter over folds
  {
    if (tune){
      best_parameters = do.call(rbind,
                                lapply(cv_results$folds,
                                       function(u){
                                         u$tuning$best$parameter
                                       }))
      if (average == "num"){ 
        tuned_param = mean(best_parameters, na.rm = T) 
      } else if (average == "int"){
        tuned_param = as.integer(round(mean(best_parameters, na.rm = T),0))
      } else {
        stop("Tuning in k-fold CV: Choose average from ('num','int')")
      }
      
      # here hard coded which parameter to change in what way. could be generalized.
      top_taxa = (sort(apply(data[,x_names], 2, sd), decreasing = T))
      keep_taxa = names(top_taxa)
      perc_sd_all = sum(top_taxa)
      x_names = keep_taxa[c(1:tuned_param)]
      perc_sd = sum(top_taxa[1:tuned_param], na.rm = T)
      perc_sd_rel = perc_sd / perc_sd_all * 100
      
      cv_results$tuning = list()
      cv_results$tuning[["x_names"]] = x_names
      cv_results$tuning[["perc_sd_rel"]] = perc_sd_rel
    } 
  }
  
  # extract mean and sd of metrics per fold 
  {
    metrics = data.frame(fold = 1:k,
                         metric = do.call(rbind, lapply(cv_results$folds, 
                                                        function(u){ u$metric })),
                         mean_mae = do.call(rbind, lapply(cv_results$folds, 
                                                          function(u){ u$mae })),
                         do.call(rbind, lapply(cv_results$folds, 
                                               function(u){ u$mae_per_class })))
    mean_metric = mean(metrics$metric, na.rm = T)
    cv_results[["mean_metric"]] = mean_metric
    sd_metric = sd(metrics$metric, na.rm = T)
    cv_results[["sd_metric"]] = sd_metric
    if (prnt){ cat("\n\nMean Metric +/- st. deviation over all ", k, " folds:", 
                   mean_metric, "+/-", sd_metric) }
    
    mean_mae_per_class = apply(metrics[,-c(1:3)], 2, mean, na.rm = T)
    sd_mae_per_class = apply(metrics[,-c(1:3)], 2, sd, na.rm = T)
    cv_results[["mean_mae_per_class"]] = mean_mae_per_class
    cv_results[["sd_mae_per_class"]] = sd_mae_per_class
    
    mean_mae = mean(mean_mae_per_class)
    sd_mae = mean(sd_mae_per_class)
    cv_results[["mean_mae"]] = mean_mae
    cv_results[["sd_mae"]] = sd_mae
  }   
  
  # train and predict final model on full set and pass with mean metric from folds
  analysis_results_full = analysis_function(..., df_train = data, df_test = NULL, 
                                            x_names = x_names, y_names = y_names,
                                            prnt = prnt, tune = F, verb = verb)
  cv_results[["model_full_set"]] = analysis_results_full$model
  cv_results[["results_full_set"]] = analysis_results_full
  cv_results[["x_names"]] = analysis_results_full$x_names
  return(cv_results)
}

## ...
##'
##' 
##' @title multi_reg_analysis
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
multi_reg_analysis = function(..., df_train, df_test, 
                              x_names, y_names, prnt = F, tune = F, verb = 1){
  
  ## initialize ====
  {
    return_list = list()
    
    dots = list(...)
    if ("prnt" %in% names(dots)){ prnt = dots$prnt }
    if ("tune" %in% names(dots)){ tune = dots$tune }
    if ("set_seed" %in% names(dots)){ set_seed = dots$set_seed }
    if ("seed_int" %in% names(dots)){ seed_int = dots$seed_int }
  }
  
  # tune nr. of taxa to include (sorted after highest variation) ====
  
  if (tune){
    if (verb >= 1){ cat("\n\tTuning optimal predictors") }
    return_list[["tuning"]] = list()
    x_names_all = x_names
    
    # do simple train-validation split from df_train
    {
      if (set_seed){
        set.seed(seed_int)
      }
      sub_sample = sample(1:nrow(df_train), ceiling(nrow(df_train)*0.8), 
                          replace = F)
      train_inner = df_train[sub_sample,]
      validation_inner = df_train[-sub_sample,]
    }
    
    for (nr_taxa in 2:max_nr_predictors){
      nr_taxa_name = paste0(nr_taxa, "_taxa")
      return_list$tuning[[nr_taxa_name]] = list()
      x_names = x_names_all[c(1:nr_taxa)]
      
      # how much of mass, variation, and correlation with y do we keep in advance
      {
        keep_sd_n = sum(top_taxa_sd[x_names], na.rm = T)
        keep_sd_rel_n = keep_sd_n / sum(keep_sd, na.rm = T) * 100
        return_list$tuning[[nr_taxa_name]][["perc_sd"]] = keep_sd_rel_n
        
        keep_mean_n = sum(top_taxa_mean[x_names], na.rm = T)
        keep_mean_rel_n = keep_mean_n / sum(keep_mean, na.rm = T) * 100
        return_list$tuning[[nr_taxa_name]][["perc_mean"]] = keep_mean_rel_n
        
        keep_cor_n = sum(top_taxa_cor[x_names,], na.rm = T)
        keep_cor_rel_n = keep_cor_n / sum(keep_cor, na.rm = T) * 100
        return_list$tuning[[nr_taxa_name]][["perc_cor"]] = keep_cor_rel_n
        
        if (verb == 2){
          cat("\n\t\t\tTry ", nr_taxa, " of (total ", length(x_names_all), ") taxa as predictors:",
              "\n\t\t\t\twith relative continental st. deviation of:", round(keep_sd_rel_n, 4),
              "\n\t\t\t\twith relative continental means of:", round(keep_mean_rel_n, 4),
              "\n\t\t\t\trelative continental abs. correlation with landcover:", round(keep_cor_rel_n, 4),
              "\n\t\t\t\t\t% compared to including all predictors.")
        }
      }
      
      fit_results = compute_multi_reg(df = train_inner, 
                                      x_names = x_names, y_names = y_names)
      model = fit_results$model
      return_list$tuning[[nr_taxa_name]][["model"]] = model
      
      # predict inner test data ====
      {
        pred_test = predict(model, newdata = validation_inner[, x_names])  
        df_pred_test = validation_inner; df_pred_test[,y_names] = pred_test
        df_actual_test = validation_inner
        
        # only y columns for metrics evaluation
        predicted = pred_test
        actual = df_actual_test[,y_names]  
        
        # shifting to > 0 % and rescale such that row sums are 100% (part of the model!)
        {
          # add lowest value per row to get rid of negative probabilites
          predicted = t(apply(predicted, 1, 
                              function(u){
                                if (any(u < 0)){
                                  u = u - min(u, na.rm = T)     
                                } else{
                                  u
                                }
                              }))
          row_sums = apply(predicted, 1, sum)
          row_sums[row_sums == 0] = 1
          predicted = predicted / row_sums 
          
          # update df's with shifted and rescaled probabilities
          df_pred_test[,y_names] = predicted
          predicted_save = cbind("Dataset_ID" = validation_inner$Dataset_ID,
                                 as.data.frame(predicted))
        }
        
        return_list$tuning[[nr_taxa_name]][["df_actual"]] = df_actual_test
        return_list$tuning[[nr_taxa_name]][["df_pred"]] = df_pred_test
      }
      
      # compute metrics for inner test data ====
      {
        mae = multi_reg_metrics(predicted = predicted, actual = actual,
                                per_class = F, metric_name = "MAE")
        mae_per_class = multi_reg_metrics(predicted = predicted, actual = actual,
                                          per_class = T, metric_name = "MAE")
        
        return_list$tuning[[nr_taxa_name]][["mae_test"]] = mae
        return_list$tuning[[nr_taxa_name]][["mae_per_class_test"]] = mae_per_class
        if (verb == 2){ cat("\n\t\t\t\tMAE: ", mae) }
      }
      
      # predict inner training data ====
      {
        pred_train = predict(model, newdata = train_inner[, x_names])  
        df_pred_train = train_inner; df_pred_train[,y_names] = pred_train
        df_actual_train = train_inner
        
        # only y columns for metrics evaluation
        predicted = pred_train
        actual = df_actual_train[,y_names]  
        
        # shifting to > 0 % and rescale such that row sums are 100% (part of the model!)
        {
          # add lowest value per row to get rid of negative probabilites
          predicted = t(apply(predicted, 1, 
                              function(u){
                                if (any(u < 0)){
                                  u = u - min(u, na.rm = T)     
                                } else{
                                  u
                                }
                              }))
          row_sums = apply(predicted, 1, sum)
          row_sums[row_sums == 0] = 1
          predicted = predicted / row_sums 
          
          # update df's with shifted and rescaled probabilities
          df_pred_train[,y_names] = predicted
          predicted_save = cbind("Dataset_ID" = train_inner$Dataset_ID,
                                 as.data.frame(predicted))
        }
        
        return_list$tuning[[nr_taxa_name]][["df_actual"]] = df_actual_train
        return_list$tuning[[nr_taxa_name]][["df_pred"]] = df_pred_train
      }
      
      # compute metrics for inner training data ====
      {
        mae = multi_reg_metrics(predicted = predicted, actual = actual,
                                per_class = F, metric_name = "MAE")
        mae_per_class = multi_reg_metrics(predicted = predicted, actual = actual,
                                          per_class = T, metric_name = "MAE")
        
        # additionally all other metrics
        return_list$tuning[[nr_taxa_name]][["mae_train"]] = mae
        return_list$tuning[[nr_taxa_name]][["mae_per_class_train"]] = mae_per_class
      }
    }
    
    # find model / nr predictors with minimum metric
    {
      # extract and plot training vs test metrics
      {
        metrics = data.frame(mae_test = do.call(rbind, 
                                                lapply(return_list$tuning, 
                                                       function(u){ u$mae_test*100})),
                             mae_train = do.call(rbind, 
                                                 lapply(return_list$tuning, 
                                                        function(u){ u$mae_train*100})))
        metrics$nr_taxa = 2:(nrow(metrics)+1) 
        colnames(metrics)[which(colnames(metrics) == "mae_train")] = "Training Set"
        colnames(metrics)[which(colnames(metrics) == "mae_test")] = "Validation Set"
        return_list$tuning[["metrics"]] = metrics
      }
      
      # find optimal model and extract model etc.
      {
        min_metric = min(metrics[["Validation Set"]], na.rm = T)
        min_ind = which.min(metrics[["Validation Set"]])
        metrics[["Validation Set"]][min_ind] == min_metric
        
        return_list$tuning[["best_model"]] = list()
        return_list$tuning$best_model[["model_ind"]] = min_ind
        return_list$tuning$best_model[["min_metric"]] = min_metric
        best_model = return_list$tuning[[min_ind]]$model
        return_list$tuning$best_model[["model"]] = best_model
        
        best_nr_taxa = metrics$nr_taxa[min_ind]
        return_list$tuning$best_model[["parameter"]] = best_nr_taxa
        best_x_names = keep_taxa[1:best_nr_taxa]
        return_list$tuning$best_model[["best_x_names"]] = best_x_names
        
        return_list$tuning$best_model[["perc_mean"]] = return_list$tuning[[min_ind]]$perc_mean
        return_list$tuning$best_model[["perc_sd"]] = return_list$tuning[[min_ind]]$perc_sd
        return_list$tuning$best_model[["perc_cor"]] = return_list$tuning[[min_ind]]$perc_cor
        
        if (verb >= 1){
          cat("\n\n\t\t Minimum MAE for current fold: ", 
              return_list$tuning$best_model$min_metric, "%", 
              "\n\t\t for including:", best_nr_taxa, "taxa", 
              "\n\t\t top variation of predictors with relative st. dev. of:", 
              round(return_list$tuning$best_model$perc_sd,4), " %",
              "\n\t\t relative prob. mass of:", 
              round(return_list$tuning$best_model$perc_mean, 4), " %", 
              "\n\t\t relative abs correlation with y:", 
              round(return_list$tuning$best_model$perc_cor, 4), "%",
              "\n\t\t\t compared to the sums with all predictors (", 
              "passed to the regression function).")
        }
      }
    }
    
    # set x_names now to the optimal predictor set to train on the full training set
    x_names = best_x_names 
  }
  
  # fit on full set ====
  {
    fit_results = compute_multi_reg(..., df = df_train, 
                                    x_names = x_names, y_names = y_names)
    model = fit_results$model
    return_list[["model"]] = model
    return_list[["formula"]] = fit_results$formula
    A = fit_results$A
    return_list[["A"]] = A
    return_list[["lc_classes"]] = fit_results$lc_classes
    return_list[["lc_class_inds"]] = fit_results$lc_class_inds
  }
  
  # predict values ====
  {
    if (!is.null(df_test)){
      
      # if test data is provided, return metric, actual, and predicted for this test data
      pred_test = predict(model, newdata = df_test[, x_names])  
      df_pred_test = df_test; df_pred_test[,y_names] = pred_test
      df_actual_test = df_test
      
      df_pred_train = NULL; df_actual_train = NULL
      
      # only y columns for metrics evaluation
      predicted = pred_test
      actual = df_actual_test[,y_names]  
      
    } else {
      
      # if no test data provided, return metric, actual, and predicted for training data
      pred_train = predict(model, newdata = df_train)
      df_pred_train = df_train; df_pred_train[,y_names] = pred_train
      df_actual_train = df_train
      
      df_pred_test = NULL; df_actual_test = NULL
      
      predicted = pred_train
      actual = df_actual_train[,y_names]
    }
    
    # shifting to > 0 % and rescale such that row sums are 100% (part of the model!)
    {
      return_list[["rescaled"]] = list()
      return_list[["rescaled"]][["max_row_sums_before"]] = max(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["min_row_sums_before"]] = min(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["max_pred_value_before"]] = max(predicted)*100
      return_list[["rescaled"]][["min_pred_value_before"]] = min(predicted)*100
      
      # add lowest value per row to get rid of negative probabilites
      predicted = t(apply(predicted, 1, 
                          function(u){
                            if (any(u < 0)){
                              u = u - min(u, na.rm = T)     
                            } else{
                              u
                            }
                          }))
      row_sums = apply(predicted, 1, sum)
      row_sums[row_sums == 0] = 1
      predicted = predicted / row_sums 
      
      return_list[["rescaled"]][["max_row_sums_after"]] = max(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["min_row_sums_after"]] = min(as.numeric(apply(predicted, 1, sum)))*100
      return_list[["rescaled"]][["max_pred_value_after"]] = max(predicted)*100
      return_list[["rescaled"]][["min_pred_value_after"]] = min(predicted)*100
      
      # update df's with shifted and rescaled probabilities
      if (!is.null(df_test)){
        df_pred_test[,y_names] = predicted
        predicted_save = cbind("Dataset_ID" = df_test$Dataset_ID,
                               as.data.frame(predicted))
      } else {
        df_pred_train[,y_names] = predicted
        predicted_save = cbind("Dataset_ID"=df_train$Dataset_ID,
                               as.data.frame(predicted))
      }
      
      
    }
    
    if (!is.null(df_test)){
      return_list[["df_actual"]] = df_actual_test
      return_list[["df_pred"]] = df_pred_test
    } else {
      return_list[["df_actual"]] = df_actual_train
      return_list[["df_pred"]] = df_pred_train
    }
  }
  
  # compute metrics ====
  {
    mae = multi_reg_metrics(..., predicted = predicted, actual = actual,
                            per_class = F, metric_name = "MAE")
    mae_per_class = multi_reg_metrics(..., predicted = predicted, actual = actual,
                                      per_class = T, metric_name = "MAE")
    
    # over all metric for cross-fold validation
    return_list[["metric"]] = mae 
    
    # additionally other metrics
    return_list[["mae"]] = mae
    return_list[["mae_per_class"]] = mae_per_class
  }
  
  
  return(return_list)
}

## ...
##'
##' 
##' @title multi_reg_metrics
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
multi_reg_metrics = function(..., predicted, actual,
                             per_class = F, metric_name = "MAE"){
  
  dots = list(...)
  if (('per_class' %in% names(dots))){ per_class = dots$per_class }
  if (('metric_name' %in% names(dots))){ metric_name = dots$metric_name }
  
  if(metric_name == "MAE"){
    resids = predicted - actual
    mae_per_class = apply(abs(resids), 2, 
                          mean, na.rm=T)
    mae = mean(mae_per_class, na.rm=T)
    if (per_class){ 
      metric = mae_per_class
    } else {
      metric = mae
    }
    
  } else if (metric_name == "KLD"){
    keras_available = require(keras)
    if (keras_available){
      kld_per_sample = as.array(keras::metric_kullback_leibler_divergence(actual, predicted))
      kld = mean(kld_per_sample, na.rm=T)
      metric = kld            
    } else { stop("Install Keras for metric KLD")}
  } else if (metric_name == "RMSE"){
    
    # default RMSE (root mean squared error) 
    resids = predicted - actual
    resids_sq = resids**2
    rmse_per_class = apply(resids_sq, 2, 
                           function(u){
                             sqrt(mean(u, na.rm = T))    
                           })
    rmse = mean(rmse_per_class, na.rm = T)
    
    if (per_class){ 
      metric = rmse_per_class
    } else {
      metric = rmse
    }
  }
  
  return(metric)
}

## ...
##'
##' 
##' @title compute_multi_reg
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
compute_multi_reg = function(..., df, x_names, y_names,
                             with_ic = T){
  
  dots = list(...)
  return_list = list()
  
  # seperate data
  {
    x = as.matrix(df[,x_names])
    y = as.matrix(df[,y_names])
    df_model = df[, c(x_names)]
  }
  
  # do regression
  {
    form = as.formula(paste0("y ~ .", c(" + 0", "")[(as.numeric(with_ic)+1)]))   
    reg = lm(form , data = df_model)
    return_list[["model"]] = reg
    return_list[["formula"]] = form
    
    A = reg$coefficients
    lc_classes = colnames(A) 
    lc_class_inds = 1:ncol(A)
    return_list[["A"]] = A
    return_list[["lc_classes"]] = lc_classes
    return_list[["lc_class_inds"]] = lc_class_inds
  }
  
  return(return_list)
}

## ...
##'
##' 
##' @title vegLcModel
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
vegLcModel = function(veg_modern, meta_cols_veg,
                      lc_modern, meta_cols_lc, 
                      prefilter_taxa = F, tune_taxa = F, 
                      with_ic = F,
                      save_folder = NULL){
  
  cols_veg = setdiff(colnames(veg_modern), meta_veg) 
  cols_lc = colnames(lc_modern %>% select(-any_of(meta_lc)))
  
  data  = full_join(veg_modern, lc_modern, by = "Dataset_ID")
  
  # run the analysis cross validated and tune nr. of predictors to include
  if (prefilter_taxa){
    # subset predictors based on mass, variation, and correlation 
    # use these as a basis predictor set for tuning the ideal number of taxa (tune = T) 
    # or skip tuning and use all these predictors (tune = F)
    # keep taxa with relevant probability mass and variation in the data set
    {
      top_taxa_mean = (sort(apply(data[,cols_veg], 2, mean), decreasing = T))
      top_mean = names(top_taxa_mean[top_taxa_mean > min_cont_mean_taxa])
      
      top_taxa_sd = (sort(apply(data[,cols_veg], 2, sd), decreasing = T))
      top_sd = names(top_taxa_sd[top_taxa_sd > min_cont_sd_taxa])
      top_mean_sd = unique(c(top_mean, top_sd))
    }
    
    # add those taxa with high abs correlation with y variables
    {
      top_taxa_cor = suppressWarnings(as.data.frame(do.call(rbind,
                                                            lapply(cols_veg,
                                                                   function(u){
                                                                     sum(abs(cor(cbind(data %>% 
                                                                                         select(as.name(u)),
                                                                                       data[,cols_lc]))[,1]), 
                                                                         na.rm = T) - 1
                                                                   }))))
      row.names(top_taxa_cor) = cols_veg; colnames(top_taxa_cor) = "sum_abs_cor_ys"
      top_taxa_cor = top_taxa_cor %>% arrange(desc(sum_abs_cor_ys))
      top_cor = row.names(top_taxa_cor %>% filter(sum_abs_cor_ys > min_sum_abs_cor_ys))
      top_mean_sd_cor = unique(c(top_mean_sd, top_cor))
    }
    
    # remove predictor with high correlations with other predictors, if existent
    {
      pred_cors = data.frame(suppressWarnings(cov(data[,top_mean_sd_cor])))
      cor_taxa = lapply(1:nrow(pred_cors),
                        function(u){
                          taxon_cor = sort(unlist(abs(pred_cors[u,-c(u,which(is.na(pred_cors[u,])))])),
                                           decreasing = T)
                          taxon_cor = taxon_cor[taxon_cor > max_cor_among_taxa]
                          ret_list = list()
                          ret_list[[colnames(pred_cors)[u]]] = names(taxon_cor)
                          return(ret_list)
                        })
      rem_taxa = c(); dont_rem_taxa = c()
      for (taxon in cor_taxa){
        dont_rem_taxa = c(dont_rem_taxa, names(taxon))
        del_taxa = taxon[[1]]
        del_taxa = del_taxa[del_taxa %in% dont_rem_taxa == F]
        rem_taxa = unique(c(rem_taxa, del_taxa))
      }
      
      top_mean_sd_cor = top_mean_sd_cor[!(top_mean_sd_cor %in% rem_taxa)]
    }
    
    # how much of mass, variation, and correlation with y do we keep (before even starting tuning)
    {
      keep_taxa = top_mean_sd_cor 
      keep_sd = sum(top_taxa_sd[keep_taxa], na.rm = T)
      keep_sd_rel = keep_sd / sum(top_taxa_sd, na.rm = T) * 100
      keep_mean = sum(top_taxa_mean[keep_taxa], na.rm = T)
      keep_mean_rel = keep_mean / sum(top_taxa_mean, na.rm = T) * 100
      keep_cor = sum(top_taxa_cor[keep_taxa, "sum_abs_cor_ys"], na.rm = T)
      keep_cor_rel = keep_cor / sum(top_taxa_cor$sum_abs_cor_ys, na.rm = T) * 100
    }
    
    cols_veg = top_mean_sd_cor
  } 
  
  # build model in cross validated (and potentially tune number of predictors if tune_taxa == T)
  {
    if (k_folds == "loo"){ k_folds = dim(data)[1] }
    cv_results = k_fold_cv(
      # CV parameters
      data = data, k = k_folds, 
      x_names = cols_veg, 
      y_names = cols_lc, 
      analysis_function = multi_reg_analysis, 
      set_seed = F, seed_int = 1234,
      prnt = T, 
      tune = tune_taxa, average = "int",
      
      # fit parameters
      with_ic = F,
      
      # other
      verb = 1
    )
    
    # extract how well prediction of props normalized to 1 and probs between 0 .. 1 worked
    {
      # cat(paste0("\nShow how well the conditions for probability ", 
      #            "distributions (compositions) are met:\n"))
      # norm_errors = do.call(rbind,lapply(cv_results$folds, function(u){
      #   do.call(cbind,u$rescaled)
      # }))
      # round(apply(norm_errors, 2, mean),4)
    }
    
    if (tune_taxa){
      # keep only predictors which was tuned to be the best number 
      # (pred_rel_variation_perc % of the total variation with all predictors)
      cols_veg = cv_results$tuning$cols_veg
      
      # plot average training vs. test curve over folds
      {
        metrics_per_nr_pred = lapply(1:min(length(cols_veg),max_nr_predictors), 
                                     function(u){
                                       ret_df = do.call(rbind, 
                                                        lapply(cv_results$folds,
                                                               function(z){
                                                                 z$tuning$metrics[u,]
                                                               }))
                                       row.names(ret_df) = NULL
                                       return(ret_df)
                                     })
        
        # mean absolute error
        {
          metrics_curve_mae = do.call(rbind,
                                      lapply(2:length(metrics_per_nr_pred),
                                             function(z){
                                               u = metrics_per_nr_pred[[z]]
                                               data.frame(
                                                 n_pred = unique(u$nr_taxa),
                                                 mean_mae_train = mean(u[["Training Set"]], 
                                                                       na.rm = T),
                                                 sd_mae_train = sd(u[["Training Set"]], 
                                                                   na.rm = T),
                                                 mean_mae_val = mean(u[["Validation Set"]], 
                                                                     na.rm = T),
                                                 sd_mae_val = sd(u[["Validation Set"]], 
                                                                 na.rm = T)
                                               )
                                             }))
          metrics_curve_mae$upper_mae_train = metrics_curve_mae$mean_mae_train + 
            metrics_curve_mae$sd_mae_train
          metrics_curve_mae$lower_mae_train = metrics_curve_mae$mean_mae_train - 
            metrics_curve_mae$sd_mae_train
          metrics_curve_mae$upper_mae_val = metrics_curve_mae$mean_mae_val + 
            metrics_curve_mae$sd_mae_val
          metrics_curve_mae$lower_mae_val = metrics_curve_mae$mean_mae_val - 
            metrics_curve_mae$sd_mae_val
          df_plot_train = metrics_curve_mae %>% 
            select((c("n_pred", ends_with("train")))) %>% 
            select(-"sd_mae_train") 
          colnames(df_plot_train) = gsub("_train", "", colnames(df_plot_train))
          df_plot_train$Dataset = "Training Set"
          df_plot_val = metrics_curve_mae %>% 
            select((c("n_pred", ends_with("val")))) %>% 
            select(-"sd_mae_val") 
          colnames(df_plot_val) = gsub("_val", "", colnames(df_plot_val))
          df_plot_val$Dataset = "Validation Set"
          df_plot = rbind(df_plot_train, df_plot_val)
          df_plot$Dataset = factor(df_plot$Dataset, levels = c("Validation Set",
                                                               "Training Set"))
          
          plt = ggplot(df_plot) +
            geom_line(aes(x = n_pred, y = mean_mae, colour = Dataset),
                      size = 0.75) +
            geom_ribbon(aes(x = n_pred, ymin = lower_mae, ymax = upper_mae, fill = Dataset), alpha = 0.4) +
            scale_x_continuous(breaks = c(2, seq(0,100,2))) + 
            scale_color_manual(values=c("indianred", "steelblue")) +
            theme_minimal() + 
            labs(x = "Nr. of predictors included in mulitvariate regression model ( )",
                 y = "Mean absolute error ( )",
                 title = "Model Evaluation on Training and Validation Set",
                 subtitle = paste0("Mean and standard deviations over ", k_folds, 
                                   " folds\nas a function of the number ", 
                                   "of predictors included."),
                 colour = "Dataset")
          print(plt)
        }
        
        # .. other metrics added
      }
    }
    
    # results from final training on full set
    {
      final_model_results = cv_results$results_full_set
      
      A = final_model_results$A
      
      df_pred = final_model_results$df_pred
      df_actual = final_model_results$df_actual
      
      lc_classes = final_model_results$lc_classes
    }
    
    # extract metrics
    {
      # on average over classes
      mean_metric = cv_results$mean_metric; sd_metric = cv_results$sd_metric
      cat("\n\nMultiple Regression without intercept:",
          "\n\tMean +/- SD of metric", round(mean_metric, 4), "+/-", 
          round(sd_metric, 4))
      #metric_full_set = cv_results$results_full_set$metric
      model_full_set = cv_results$model_full_set
      
      # per lc class
      mae_per_class = cv_results$mean_mae_per_class; names(mae_per_class) = cols_lc
      sd_per_class = cv_results$sd_mae_per_class; names(sd_per_class) = cols_lc
    }
  }
  
  # find contribution of taxa per land cover class to landcover shares
  {
    # extract siginificant taxa:
    {
      model_full_set_sum = summary.lm(model_full_set)
      sig_per_class = coef(summary(model_full_set))
      sig_taxa = lapply(sig_per_class, 
                        function(u){
                          inds = which(u[,4]<=sig_level)
                          names(inds)
                        })
      names(sig_taxa) = gsub("Response ", "", names(sig_taxa))
    }
    
    # plot influence of vegetation cover (taxa) on landcover
    {
      # rescale coefficient matrix for plotting (keep pos and negative influences)
      {
        na_inds = which(apply(A, 1, function(u){all(is.na(u))}))
        if (length(na_inds) > 0){ 
          A_rescaled = A[-na_inds,] 
        } else {
          A_rescaled = A
        }
        cols_veg = cols_veg[!(cols_veg %in% names(na_inds))]
        
        # chose top variation taxa
        A_rescaled = A_rescaled[1:min(nrow(A_rescaled), 
                                      show_top_var_taxa_in_matrix), ]
        A_rescaled = A_rescaled * 100 / max(abs(A_rescaled))
      }
      
      # plot taxa to landcover_modern calib translation matrix
      {
        translation_df = t(as.data.frame(A_rescaled))
        translation_df = cbind(data.frame("Landcover_Class" = row.names(translation_df)), 
                               translation_df)
        row.names(translation_df)=NULL
        
        translation_df = pivot_longer(translation_df, names_to = "Taxon", 
                                      cols = all_of(colnames(translation_df)[-1]), 
                                      values_to = "LC_Share")
        # add information if taxa significant for lc class
        {
          translation_df = do.call(rbind, lapply(1:length(lc_classes),
                                                 function(u){
                                                   df1 = translation_df %>% 
                                                     filter((Landcover_Class == 
                                                               lc_classes[u]) & 
                                                              (Taxon %in% sig_taxa[[u]])) %>%
                                                     mutate(Significant = T)
                                                   df2 = translation_df %>% 
                                                     filter((Landcover_Class == 
                                                               lc_classes[u]) &
                                                              !(Taxon %in% sig_taxa[[u]])) %>% 
                                                     mutate(Significant = F)
                                                   return(rbind(df1, df2))
                                                 }))
        }
        
        # veg_modern abundances over all sites
        landcover_region_means = data.frame("Class" = colnames(lc_modern %>% select(-Dataset_ID)),
                                            "regional_Mean" = as.numeric(apply(lc_modern %>% select(-Dataset_ID), 2, 
                                                                               mean, na.rm=T)),
                                            "regional_Variance" = as.numeric(apply(lc_modern %>% select(-Dataset_ID), 2, 
                                                                                   var, na.rm=T))) %>%
          arrange(desc(regional_Mean)) %>%
          filter(regional_Mean >= 1)
        landcover_region_means$Class = factor(landcover_region_means$Class, levels = landcover_region_means$Class)
        landcover_region_means[,"regional_Relative_StD"] = sqrt(landcover_region_means$regional_Variance) /
          landcover_region_means$regional_Mean * 100 # %
        
        lc_defs_plot = class_defs %>% 
          filter(Class_Code %in% levels(landcover_region_means$Class)) %>%
          select(all_of(c("Class_Code", "Class_Plotlabel", "Color_Code"))) %>%
          arrange(Class_Code)
        labels = paste0(lc_defs_plot$Class_Code,":\n",
                        lc_defs_plot$Class_Plotlabel)
        translation_df$Landcover_Class = factor(translation_df$Landcover_Class, 
                                                levels = unique(c(translation_df$Landcover_Class, 
                                                                  as.character(lc_defs_plot$Class_Code))))
        translation_df$Taxon = factor(translation_df$Taxon,
                                      levels=cols_veg)
        translation_df$Coef_Label = as.character(round(translation_df$LC_Share,0))
        plt = ggplot(translation_df, aes(x = Taxon, y = Landcover_Class, 
                                         fill=LC_Share, shape = Significant)) + 
          geom_tile() + 
          geom_point(size=10, colour = "grey75") + 
          scale_shape_manual(values=c(NA, 16)) +
          geom_text(aes(label=Coef_Label)) +
          labs(shape = paste0("Signifcant\n(p < ",sig_level,")"), 
               x = "Taxon", y = "Land cover Type",
               title = paste0("Taxon Contributions to Land Cover Types"),
               subtitle = paste0("Regression coefficients normalized by the highest ", 
                                 "absolute value.\nOnly the ratios between taxa shown are ", 
                                 "meaningfull, not the absolute values.",
                                 "\n\nValues of significant regression coefficients ", 
                                 "for taxa in grey circles."),
               fill = "Resulting\nLC share (%)",
               fill = "Coefficient",
               y = "Land cover Type", x = "Taxon") +
          scale_y_discrete(labels = labels) + 
          scale_fill_gradient2(low = "steelblue", mid="grey95", high="#d1495b") +
          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                panel.grid.major.x = element_blank(), 
                panel.grid.major.y = element_blank(),
                panel.grid.minor=element_line(colour="grey20")) 
        
        plt
      }
    }
  }
  
  return_list = list(model_full_set = model_full_set, cols_veg = cols_veg,
                     veg_lc_mat = A,
                     mean_metric = mean_metric, sd_metric = sd_metric,
                     mae_per_class = mae_per_class, sd_per_class = sd_per_class)
  
  if (!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, "vegLcModel.rda"))
  } else{
    return(return_list)
  }
}

## ...
##'
##' 
##' @title predictLcExplOnly
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
predictLcExplOnly = function(model, veg_predictors, meta_cols){
  
  # predict the landcover that is explained by the vegetation cover
  {
    cols_veg = colnames(veg_predictors %>% select(-all_of(meta_cols)))
    x = veg_predictors[, cols_veg]
    y_pred = cbind(veg_predictors[, meta_cols],
                   predict(model, x)) %>% 
      filter(rowSums(.)!=0)
  }
  
  # shifting to minimum 0 % and rescaling to 100% row sums is part of the model
  {
    lc = y_pred %>% select(-all_of(meta_cols))
    
    # add lowest value per row to get rid of negative probabilites
    lc = t(apply(lc, 1, 
                 function(u){
                   if (any(u < 0, na.rm = T)){
                     u = u - min(u, na.rm = T)     
                   } else{
                     u
                   }
                 }))
    
    row_sums = apply(lc, 1, sum, na.rm = T)
    row_sums[row_sums == 0] = 1
    lc = lc / row_sums
    row_sums = apply(lc, 1, sum, na.rm = T)
    
    y_pred = cbind(y_pred[, meta_cols], lc)
  }
  
  x = data.frame(Dataset_ID = veg_predictors$Dataset_ID,
                 Age = veg_predictors$Age,
                 x)
  
  return_list = list(x = x, y_pred = y_pred)
  return(return_list)
}


## ...
##'
##' 
##' @title predictLcAllClasses
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
predictLcAllClasses = function(lc_veg_only, lc_all_classes, verb = 0){
  
  # normalize extracted land cover shares in PSAs or percentages to 1
  lc_all_classes[-1] = lc_all_classes[,-1] / apply(lc_all_classes[, -1], 1, sum) 
  
  IDs = unique(lc_veg_only$Dataset_ID)
  for (ID in IDs){
    
    # filter modern with all classes
    {
      y_all_site = lc_all_classes %>% filter(Dataset_ID == ID) %>% mutate(Origin = "Actual", Age = -1000)  # place holder age for the modern period
      y_all_site = y_all_site[,c("Dataset_ID", "Age", "Origin", 
                                 colnames(y_all_site %>% select(-all_of(c("Dataset_ID", "Age", "Origin")))))]
    }
    
    # filter only classes explained by vegetation for actual extracted (from land cover map) and predicted
    {
      y_veg_pred_site = lc_veg_only %>% filter(Dataset_ID == ID) %>% mutate(Origin = "Predicted")
      y_veg_actual_site = y_all_site %>% select(all_of(c("Dataset_ID", "Age", "Origin", as.character(include_calib_codes)))) 
      site_df = rbind(y_veg_actual_site, y_veg_pred_site)
      site_df = site_df[, c("Dataset_ID", "Age", "Origin", as.character(include_calib_codes))]
    }
    
    # renormalize shares
    {
      # recalculate probs with all classes included in calibration and prediction
      # and add excluded (either with constant probs or set to zero)
      # 1st go back to actual share in psa of vegetation classes ('includeded in calibration')
      {
        #apply(site_df[,-c(1:3)], 1, sum)
        perc_calib = sum(y_all_site[, as.character(include_calib_codes)])
        site_df[,-c(1:3)] = site_df[,-c(1:3)] * perc_calib /
          (apply(site_df[,-c(1:3)], 1, sum))
        apply(site_df[,-c(1:3)], 1, sum)
        
        y_all_site = (y_all_site %>% select(-Dataset_ID, -Age, -Origin))
        sum(y_all_site)
      }
      
      # 2nd share of excluded classes
      {
        excl_classes_present = as.character(exclude_calib_codes)[which(as.character(exclude_calib_codes) %in% colnames(lc_all_classes))]
        y_exc_site = y_all_site[,excl_classes_present]
        perc_exc = sum(y_exc_site)
        abs(1 - perc_calib - perc_exc) <= 10**(-3)
      }
      
      # 3rd now differentiate excluded classes, i.e. set defined classes to zero
      # others keep constant and renormalize all
      {
        const_cols_present = as.character(keep_const_codes)[which(as.character(keep_const_codes) %in% 
                                                                    colnames(y_all_site))]
        y_exc_const_site = y_all_site[,const_cols_present] 
        perc_keep_constant = sum(y_exc_const_site) 
        
        zero_cols_present = as.character(set_zero_codes)[which(as.character(set_zero_codes) %in% 
                                                                 colnames(y_all_site))]
        y_exc_zero_site = y_all_site[,zero_cols_present]
        perc_set_zero = sum(y_exc_zero_site) 
        abs(perc_exc - (perc_keep_constant + perc_set_zero)) <= 10**(-3)
        y_exc_zero_site[,names(y_exc_zero_site)] = 0
      }
      
      # 4th renormalize predicted including the set zero columns 
      # i.e. multiply with amplification factor
      {
        y_pred_new_site = cbind(site_df %>% filter(Origin == "Predicted"), 
                                y_exc_zero_site)
        all(abs(apply(y_pred_new_site[,-c(1:3)], 1, sum) - perc_calib) <= 10**(-3))
        amp_factor = ((perc_calib + perc_set_zero) / perc_calib)
        y_pred_new_site[,-c(1:3)] = y_pred_new_site[,-c(1:3)] *
          (amp_factor)
        
        all(apply(y_pred_new_site[,-c(1:3)],1,sum) + perc_keep_constant/100 - 1 <= 10**(-3))
        
        if (verb > 0){
          cat(paste0("\n\tSite ",ID,": \n\t\tIncreased share of predicted lc types",
                     "\n\t\tdue to setting others to zero in the past",
                     "\n\t\tby a factor of ", 
                     round(amp_factor,4)))
        }
      }
      
      # 5th add columns kept constant in share
      {
        y_pred_new_site = cbind(y_pred_new_site,
                                y_exc_const_site)
        all_cols_present = as.character(all_lc_codes)[which(as.character(all_lc_codes) %in% 
                                                              colnames(y_pred_new_site))]
        y_pred_new_site = y_pred_new_site[,c("Dataset_ID", "Age", all_cols_present)]
        all(abs(apply(y_pred_new_site[,-c(1:2)], 1, sum) - 1) <= 10**(-3))
      }
      
      # do some checks of the landcover shares
      {
        all(abs(apply(y_pred_new_site[,-c(1:2)],1,sum) - 
                  (perc_calib + perc_set_zero + perc_keep_constant)) <= 10**(-3))
        all(abs(apply(y_pred_new_site[,const_cols_present],1,sum) - 
                  perc_keep_constant) <= 10**(-3))
        all(abs(apply(y_pred_new_site[,zero_cols_present],1,sum)) <= 10**(-3))
        all(abs(apply(y_pred_new_site[,as.character(include_calib_codes)], 1, sum) - 
                  (perc_calib + perc_set_zero)) <= 10**(-3))
      }
      
      # add the modern time slice 
      # (that has not been modified by setting classes to zero and amplifying the other classes)
      {
        y_pred_new_site = 
          rbind(y_all_site %>% mutate(Dataset_ID = ID, Age = -1000, .before = 1),
                y_pred_new_site)
        
        # and check set to zero, keep constant, and amplification is correct
        #y_pred_new_site[, as.character(const_cols_present)] # all slices (actual and predicted) should have the same values
        #y_pred_new_site[, as.character(zero_cols_present)] # only the modern actual slice should have values
        apply((y_pred_new_site[, as.character(include_calib_codes)]), 1, sum)[2] / 
          apply((y_pred_new_site[, as.character(include_calib_codes)]), 1, sum)[1] == amp_factor
      }
      
      # add to data frame with all sites
      if (ID == IDs[1]){
        y_pred_new_full = y_pred_new_site
      } else {
        y_pred_new_full = rbind(y_pred_new_full, y_pred_new_site)
      }
    }
  }
  
  return_list = list(y_pred_all_classes = y_pred_new_full)
  return(return_list)
}

## ...
##'
##' 
##' @title predictLC
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
predictLC = function(model, veg_predictors, meta_cols, lc_modern,
                     save_folder = NULL){
  
  # predict only the land cover classes predicted by the vegetation cover (normalized to 1)
  results = predictLcExplOnly(model = model, 
                              veg_predictors = veg_predictors, 
                              meta_cols = meta_cols)
  lc_ts_veg_only = results$y_pred
  veg_cover_predictors = results$x # if filtered and/or tuned the predictors are less (columns) than the original input to model building
  rm(results)
  
  
  # calculate real land cover shares within in pollen source areas (as for model building only vegetational classes were normalized to 1)
  # and set defined classes to zero in the past (e.g. croplands) or keep constant (e.g. water bodies)
  results = predictLcAllClasses(lc_veg_only = lc_ts_veg_only, 
                                lc_all_classes = lc_modern,
                                verb = 0)
  lc_ts_all = results$y_pred_all_classes
  rm(results)
  
  
  return_list = list(lc_ts_all = lc_ts_all)
  
  if (!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, "predictLC.rda"))
  } else{
    return(return_list)
  }
}

## ...
##'
##' 
##' @title extractDataPSA
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
extractDataPSA = function(poly, landcover, elevation, climate,
                          save_folder, del_classes = NULL){
  
  ID = poly$Dataset_ID
  
  # check if the pollen source area crosses the date line
  lons = as.data.frame(st_coordinates(poly))$X
  if (sign(min(lons)) != sign(max(lons))){
    # no fix implemented yet, but I think the easiest solution is to shift longitudes (centered in the PSA center, i.e shift the dateline do not tranform the data), 
    # extract the data (still in wgs84) and shift the longitudes back. I have code that as I did this for Chenzhi, ask me if needed.
    # basically the st_shift_longitude function is limited to shift 180° and there is no inverse function. I generalized it.
    # for now skip the site
    cat("\tSkip this site, repair extraction of data for sites crossing the date line")
    next
  }
  
  # extract landcover
  {
    lc_site = suppressMessages((landcover[poly,,] %>% st_as_stars())[]) # keep this object for dimensions below
    lc_site = suppressWarnings(
      lc_site %>% 
        st_as_sf(crs = 4326) %>% 
        st_centroid()
    )
    coords = as.data.frame(st_coordinates(lc_site))
    df_site = as.data.frame(lc_site %>% st_drop_geometry())
    df_site$cell_ind = paste0(ID, "_", 1:nrow(df_site))
    df_site$x = coords$X; df_site$y = coords$Y; rm(coords)
    df_site = df_site[,c("cell_ind", "x", "y", "Landcover")]
  }
  
  # extract env. variables with the resolution of the land cover product")
  {
    df_site = cbind(df_site,
                    do.call(cbind, 
                            lapply(c(list(elevation), climate),
                                   function(u) {
                                     #cat("\n\t\tfrom", names(u))
                                     var_site = suppressMessages(
                                       as.data.frame(
                                         st_extract(u, lc_site$geometry) %>%
                                           st_drop_geometry()
                                       )) 
                                     #cat("\tdim: ", dim(var_site))
                                     return(var_site)
                                   })))
    df_site$geometry = lc_site$geometry
  }        
  
  colnames(df_site) = c("cell_ind", "x", "y", "Landcover", "Elevation", names(climate), "geometry")
  
  # edit NAs
  {
    df_site = df_site %>% mutate(across(all_of(c("gdd0", "gdd5", "gdd10")), 
                                        ~replace_na(.x, 0)))
    apply(df_site %>% select(where(is.numeric)), 2, function(u){ sum(is.na(u)) })
    
    # drop remaining rows with any NA-value
    df_site = df_site %>% filter(rowSums(is.na(.)) == 0)
  }
  
  # merge lc classes
  {
    cat("\n\t\tMerge to class:")
    for (merge_class in classes_to_merge){
      class_inds = which(df_site$Landcover %in% merge_class$old_classes)
      df_site$Landcover[class_inds] = merge_class$new_class
      cat("\n\t\t", merge_class$new_class, "-", merge_class$new_label,
          "-", length(class_inds), "cells.")
    }
  }
  
  # set classes removed completely (del_classes, see section 'lc_modern_shares' above) to NA, i.e. here, remove the rows
  {
    if (!is.null(del_classes)){
      df_site = df_site %>% filter(!(Landcover %in% del_classes))
    }
  }
  
  save(df_site, file = file.path(save_folder, paste0("extractDataPSA_", ID, ".rda")))
}


## ...
##'
##' 
##' @title extractDataPSAs
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
extractDataPSAs = function(polys, lc_env_folder, clim_flnames, lc_flname, elev_flname,
                           save_folder = NULL, del_classes = NULL){
  
  library(stars); library(sf); sf_use_s2(F)
  library(dplyr)
  
  # load land cover and environmental data sets
  {
    landcover = read_stars(file.path(lc_env_folder, lc_flname)) %>%
      st_set_crs(4326) %>% 
      setNames("Landcover")
    
    elevation = read_stars(file.path(lc_env_folder, elev_flname)) %>% 
      st_set_crs(4326) %>% setNames("Elevation")
    
    # climate
    {
      climate = lapply(clim_flnames, 
                       function(u) {
                         clim_var = read_stars(file.path(lc_env_folder, u)) %>%
                           st_set_crs(4326) 
                         clim_var = setNames(clim_var, u)
                         return(clim_var)
                       })
      
      names(climate) = clim_vars
      climate = climate[clim_vars]
      }
    
    env_vars = c("Elevation", clim_vars)
  }
  
  IDs = polys$Dataset_ID
  cat("\nExtract landcover and environmental data:")
  for (ID in IDs){
    cat("\n\n\tID: ", ID)
    extractDataPSA(poly = polys %>% filter(Dataset_ID == ID),
                   landcover = landcover, elevation = elevation, climate = climate,
                   save_folder = save_folder, del_classes = del_classes)
  }
}

## ...
##'
##' 
##' @title downsampleForRDA
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
downsampleForRDA = function(lc_env_df,
                            save_folder = NULL){
  
  downsample = function(df, class_col){
    n = min(table(df[[class_col]]))
    idx = unlist(tapply(1:nrow(df),df[[class_col]],sample,n))
    df[idx,]
  }
  
  lc_stats = lc_env_df %>% group_by(Landcover) %>% 
    summarize(Nr_Cells = length(Landcover)) %>%
    mutate(Share_Cells = Nr_Cells / sum(Nr_Cells)) %>% 
    arrange(Landcover)
  
  lc_env_df$Landcover = factor(lc_env_df$Landcover, levels = sort(unique(lc_env_df$Landcover)))
  lc_env_downs = downsample(df = lc_env_df, class_col = "Landcover")
  
  lc_stats_downs = lc_env_downs %>% 
    group_by(Landcover) %>% 
    summarize(Nr_Cells = length(Landcover)) %>%
    mutate(Share_Cells = Nr_Cells / sum(Nr_Cells))
  
  return_list = list(lc_env_downs = lc_env_downs, 
                     lc_stats_orig = lc_stats, lc_stats_downs = lc_stats_downs)
  
  if (!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, paste0("downsampleForRDA_", ID, ".rda")))
  } else{
    return(return_list)
  }
}

## ...
##'
##' 
##' @title compute_rda
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
compute_rda = function(df, env_cols, significance_tests = F,  
                       plt_model, width, height){
  
  library(vegan)
  
  # one-hot-encode land cover type
  {
    lc_type = df %>% select(cell_ind, "Landcover") 
    lc_type$Landcover = as.numeric(as.character(lc_type$Landcover))
    lc_type = lc_type %>% as_tibble() %>% 
      pivot_wider(names_from = Landcover,
                  values_from = Landcover) %>%
      unnest(all_of(colnames(.)))
    colnames(lc_type) = c("cell_ind", paste0("LC", colnames(lc_type)[-1]))
    
    lc_type = cbind(lc_type$cell_ind,
                    as.data.frame(t(apply(lc_type[,-1], 1,
                                          function(u){
                                            u[!is.na(u)] = 1
                                            u[is.na(u)] = 0
                                            u
                                          }))))
    lc_type = as.matrix(lc_type[,-1])
  }
  
  # ‘lc’: orthogonal linear combinations of the explanatory variable
  # ‘wa’: more robust to noise in the environmental variables but are a step between constrained towards unconstrained.
  #     Most of the time, the default is “wa”. Because these are the most robust to noise in the data.
  # Scaling for species and site scores: 
  #     Either site (1)  
  #     or     species (2)
  #     scores are scaled by eigenvalues, and the other set of scores is left unscaled, or with 3 both are scaled symmetrically by square root of eigenvalues.
  
  rda_scale = F
  
  #formula = df[, env_cols] ~ df[, "Landcover"], 
  formula = as.formula(paste0("df[,env_cols] ~ ", paste0("LC", levels(df$Landcover), collapse = " + ")))
  rda_model = rda(formula = formula, 
                  data = as.data.frame(lc_type), #df[, c(env_cols,"Landcover")],
                  scale = rda_scale)
  
  # for comparison do a pca
  {
    pca_model = prcomp(df[,env_cols], center = F, scale. = F, tol = 10**(-2))
    pca_variance = pca_model$sdev / sum(pca_model$sdev)
    pca_A = pca_model$rotation
    pca_pred = pca_model$x
  }
  
  # extract model results
  {
    R2_adj = RsquareAdj(rda_model)$adj.r.squared
    eigen_vals_all = eigenvals(rda_model, "all")
    eigen_vals_constr = eigenvals(rda_model, "constrained")
    eigen_vals_unconstr = eigenvals(rda_model, "unconstrained")
    variance_df = data.frame("Axis" = c(names(eigen_vals_constr), names(eigen_vals_unconstr)), 
                             Variance = c(eigen_vals_constr, eigen_vals_unconstr),
                             Relative_Variance = c(eigen_vals_constr / sum(eigen_vals_all), 
                                                   eigen_vals_unconstr / sum(eigen_vals_all))); row.names(variance_df) = NULL
    
    # see how prediction is done and how the different scalings change euclidic distances in env space
    #vars_transf_rds_wa = as.data.frame(rda_model$CCA$wa) ==  predict(rda_model, new_data = df[, env_cols], type = "wa")
    vars_transf_rds_wa_sc_1 = as.data.frame(scores(rda_model, choices = c(1,2), display = "wa", scaling = 1))
    #vars_transf_rds_wa_sc_2 = as.data.frame(scores(rda_model, choices = c(1,2), display = "wa", scaling = 2))
    
    # check how scaling types differ in scaling over RDA axes
    #apply(round(vars_transf_rds_wa_sc_1 / vars_transf_rds_wa[,1:2], 6), 2, unique)
    #apply(round(vars_transf_rds_wa_sc_2 / vars_transf_rds_wa[,1:2], 6), 2, unique)
    
    transl_lc_rds = as.data.frame(coef(rda_model))
    row.names(transl_lc_rds) = gsub("df[, y]", "LC", row.names(transl_lc_rds), fixed = T)
    
    # matrix that translates environmental data to RDA axes
    transl_env_rds = as.data.frame(predict(rda_model, new_data = df[, env_cols], type = "sp"))
  }
  
  # check collinearity and whether to add all predictors
  {
    lc_type_vifs = sqrt(vif.cca(rda_model))
    #cat("\n\t", all(sqrt(vif.cca(rda_model)) < 2, na.rm = T)) # rule of thumb for collinearity too high
  }
  
  # fitted and residuals on env. variables 
  {
    residuals = as.data.frame(vegan:::residuals.cca(rda_model))
    predicted = as.data.frame(predict(rda_model, new_data = df[, env_cols], 
                                      type = "response")) # scaling has no effect for env. variables, only when predicting RDA values
    actual = df[, env_cols]
    #cat("\n\t", all(round(actual - (residuals + predicted), 6) == 0)) 
    
    # predicted values can also be directly assessed from the model via fitted function
    #cat("\n\t", all(round(fitted(rda_model, type = "response", scaling = 1) - predicted, 6) == 0))
  }
  
  # significance tests
  if (significance_tests){
    
    # perform significance tests
    {
      sig_global = anova.cca(rda_model)
      sig_axes = anova.cca(rda_model, by = "axis")
      sig_axes[which(sig_axes$`Pr(>F)` <= p_val),]
      sig_terms = anova.cca(rda_model, by = 'terms', step = 100)
    }
    
    # build data frame with axes' significance
    {
      sig_df = (as.data.frame(sig_terms) %>% select(all_of("Pr(>F)")))
      sig_df$P_LC_Axes = row.names(sig_df)
      colnames(sig_df) = c("P", "Axes"); row.names(sig_df) = NULL
      sig_df = sig_df[-nrow(sig_df), ]; sig_df = sig_df[,c("Axes", "P")]
      
      sig_df = rbind(sig_df, 
                     as.data.frame(sig_axes) %>% 
                       mutate("Axes" = row.names(as.data.frame(sig_axes)),
                              P = !!as.name("Pr(>F)"),
                              .keep = "none"))
      sig_df = sig_df[-nrow(sig_df),]; row.names(sig_df) = NULL
      
      sig_df = rbind(data.frame(Axes = "Model",
                                P = sig_global$`Pr(>F)`[1]), sig_df)
      sig_df = sig_df %>% mutate(Significant = (P <= p_val))
    } 
  } else{
    sig_df = NULL
  }
  
  # plotting
  if (plt_model){
    # ‘wa’ plotting 
    # weighted sums of env. variables
    {
      # scaling of scores
      {
        env_scores_sc_1 = scores(rda_model, choices=1:2, scaling = 1, display="sp")
        env_scores_sc_2 = scores(rda_model, choices=1:2, scaling = 2, display="sp") 
        
        # these are just rescaled by a different factor (constant factor per RDA axis)
        transl_env_rds[,1:2] / env_scores_sc_1
        transl_env_rds[,1:2] / env_scores_sc_2
      }
      
      # scaling 1 - distance biplot
      # euclidic distances between objects, but
      # only vectors of response variables (env. vars) and 
      # explanatory variables (lc type) reflect linear correlation
      # i.e. do not interpret correlation between env. variables
      {
        plot(rda_model, scaling=1, main="Triplot env. vars ~ lc - scaling 1 - wa scores")
        arrows(0,0,env_scores_sc_1[,1], env_scores_sc_1[,2], length=0, lty=1, col='darkred')
        text(x = env_scores_sc_1[,1], y = env_scores_sc_1[,2], labels = row.names(env_scores_sc_1),
             col = "darkred")
      }
      
      # scaling 2 - correlation biplot
      # no euclidic distances  between objects, but
      # all vectors reflect linear correlation
      # and length of env. vars arrows reflext importance for RDA (species scores are scaled by eigenvalues)
      {
        plot(rda_model, main="Triplot env. vars ~ lc - scaling 2 - wa scores")
        arrows(0,0,env_scores_sc_2[,1], env_scores_sc_2[,2], length=0, lty=1, col='darkred')
        text(x = env_scores_sc_2[,1], y = env_scores_sc_2[,2], labels = row.names(env_scores_sc_2),
             col = "darkred")
      }
    }
  }
  
  return(list(rda_model = rda_model, 
              R2_adj = R2_adj, 
              transl_lc_rds = transl_lc_rds, 
              vars_transf_rds_wa_sc_1 = vars_transf_rds_wa_sc_1,
              variance_df = variance_df, significance_df = sig_df, 
              pca_model = list(model = pca_model,
                               rel_variances = pca_variance, 
                               A = pca_A, 
                               predicted = pca_pred)))
}

## ...
##'
##' 
##' @title rdaModel
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
rdaModel = function(lc_env_rda, meta_cols, 
                    significance_tests = F, p_val = 0.05,
                    plt_model = F, width = 1200, height = 800,
                    save_folder = NULL){
  
  env_cols = colnames(lc_env_rda %>% select(-any_of(c("Landcover", meta_cols))))
  
  lc_env_rda$Landcover = factor(lc_env_rda$Landcover,
                                levels = sort(unique(lc_env_rda$Landcover)))
  
  # plot correlations among environmental data
  {
    cor_mat = cor(lc_env_rda[, env_cols]) 
    
    library(ggcorrplot)
    plt = ggcorrplot(cor_mat,
                     hc.order = F,
                     type = "upper",
                     lab = T, 
                     colors = c("#6D9EC1", "white", "#E46726")) + 
      labs(title = paste0("Correlation of Env. Variables for Site ", ID)) + 
      theme_minimal()
    plt
  }  
  
  # z-standardization
  {
    means_per_var = apply(lc_env_rda %>% select(all_of(env_cols)), 2, mean)
    sd_per_var = apply(lc_env_rda %>% select(all_of(env_cols)), 2, sd)
    
    lc_env_rda = cbind(lc_env_rda[, meta_cols],
                       do.call(cbind,
                               lapply(env_cols,
                                      function(u){
                                        (lc_env_rda[, u] - means_per_var[u]) / sd_per_var[u]
                                      })),
                       Landcover = lc_env_rda$Landcover)
    colnames(lc_env_rda) = c(meta_cols, env_cols, "Landcover")
    
    # check
    #round(apply(lc_env_rda %>% select(all_of(env_cols)), 2, mean, na.rm = T), 2)
    #round(apply(lc_env_rda %>% select(all_of(env_cols)), 2, sd, na.rm = T), 2)
  }
  
  # compute RDA model
  rda_results = compute_rda(df = lc_env_rda,
                            env_cols = env_cols,
                            significance_tests = significance_tests,
                            plt_model = plt_model,
                            width = width, height = height)
  
  if (!is.null(save_folder)){
    save(rda_results, file = file.path(save_folder, paste0("rdaModel_", unique(lc_env_rda$Dataset_ID), ".rda")))
  } else{
    return(rda_results)
  }
}


## ...
##'
##' 
##' @title plot_env_space
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
plot_env_space = function(df, class_defs, rda_model,
                          variance_df = NULL, 
                          sig_df = NULL, p_val = NULL,
                          title = "Environmental Space", RDA1 = "RDA1", RDA2="RDA2",
                          xlims = NULL, ylims = NULL){
  
  # rda results
  {
    coords = df %>% select("RDA1", "RDA2")
    
    model_prop_1 = as.numeric(variance_df %>% filter(Axis == "RDA1") %>% select(Relative_Variance)) * 100
    model_prop_2 = as.numeric(variance_df %>% filter(Axis == "RDA2") %>% select(Relative_Variance)) * 100
    model_prop_unit = "% of total variance"
    constr_var = paste0(round(sum(eigenvals(rda_model, "constrained")) / 
                              sum(eigenvals(rda_model, "all")) * 100, 1), 
                        model_prop_unit)
    unconstr_var = paste0(round(sum(eigenvals(rda_model, "unconstrained")) / 
                                sum(eigenvals(rda_model, "all")) * 100, 1), 
                          model_prop_unit)
  }
  
  # class definitions for this site
  {
    class_defs_site = class_defs %>% filter(Class_Code %in% (df$Landcover))
    
    # these defintions are used for plotting etc. throughout the script
    class_codes = class_defs_site$Class_Code
    class_names = paste0("LC", class_codes)
    class_labels = paste0("LC",class_defs_site$Class_Code, " ", 
                          gsub(", ", ",\n", class_defs_site$Class_Plotlabel))
    class_cols = class_defs_site$Color_Code
  }
  
  # prepare plot frame
  {
    data = data.frame("x"=coords[["RDA1"]], "y"=coords[["RDA2"]],
                      "Class"=df$Landcover, 
                      stringsAsFactors = F)
    data$Class = factor(data$Class, levels = class_codes)
    n_samples = nrow(data)
    rand_inds = sample(1:nrow(data),nrow(data))
    data = data[rand_inds,]; row.names(data) = NULL

    n_samples_per_class = sapply(levels(data$Class), 
                                 function(u) length(which(data$Class == u)))
    variances_pcs = variance_df[grep("PC", variance_df$Axis),]
    variances_rds = variance_df[grep("RD", variance_df$Axis),]
  }
  
  # plot
  {
    plt = ggplot(data, aes(x, y)) +
      geom_point(aes(x, y, colour = Class), alpha=0.75,
                 shape=20, size = 1.75) +
      scale_color_manual(values=class_cols, labels = class_labels) +
      labs(x = paste0(gsub("_"," ","RDA1"), " (", 
                      round(model_prop_1,2), model_prop_unit,")", sep=""),
           y = paste0(gsub("_"," ","RDA2"), " (", 
                      round(model_prop_2,2),  model_prop_unit,")", sep=""),
           title=title) + 
      labs(subtitle = paste0("\nNr. of samples for model building: ", n_samples, "\n",
                             "Percentage of samples per class:\n\t", 
                             paste0("LC", names(n_samples_per_class), " ",
                                    paste0(round(n_samples_per_class / n_samples * 100, 1), "%"),
                                    collapse = " | "), "\n",
                             "\nConstrained Variance: ", constr_var, "\n",
                             "R^2 and adjusted R^2 from Multivariate Regression: ", 
                             round(sum(variances_rds$Variance) / sum(variance_df$Variance), 4), " ",  
                             round(rda_R2_adj, 4), "\n",
                             "Variance per environmental axis:\n\t", 
                             paste0(variances_rds$Axis, " ", 
                                    round(variances_rds$Relative_Variance*100, 2), "%", collapse = " | "), "\n","\nUnconstrained Variance: ", unconstr_var, "\n",
                             "Variance per unconstrained axis:\n\t", 
                             paste0(variances_pcs$Axis, " ", 
                                    round(variances_pcs$Relative_Variance*100, 2), "%",
                                    collapse = " | "), "\n",
                             c("", paste0(
                                 paste0("\nSignificance: p <= ", p_val),
                                 paste0("\nModel Significance: ", sig_df$P[which(sig_df$Axes == "Model")],
                                        c("", " (sig)")[as.numeric(sig_df$P[which(sig_df$Axes == "Model")] <= p_val) + 1]),
                                 paste0("\nConstraining Variables Significance:"),
                                 paste0("\n\t",paste0(sig_df[grep("LC", sig_df$Axes), "Axes"]," ",
                                                      sig_df[grep("LC", sig_df$Axes), "P"], 
                                                      c("", " (sig)")[as.numeric(sig_df[grep("LC", sig_df$Axes),
                                                                                        "Significant"]) + 1], 
                                                      collapse = " | ")),
                                 paste0("\nEnvironmental Axes Significance:"), 
                                 paste0("\n\t", paste0(sig_df[grep("RD", sig_df$Axes), "Axes"]," ",
                                                       sig_df[grep("RD", sig_df$Axes), "P"], 
                                                       c("", " (sig)")[as.numeric(sig_df[grep("RD", sig_df$Axes), 
                                                                                         "Significant"]) + 1],
                                                       collapse = " | ")))
                               )[as.integer(!is.null(sig_df))+1]
                             ),
           colour = "Land Cover\nType") +
      theme_minimal() +     
      theme(legend.title = element_text(size=16), #change legend title font size
            legend.text = element_text(size=14), #change legend text font size
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 18),
            plot.subtitle = element_text(size = 12)) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) +
      list(
        if (!is.null(xlims)){
          suppressMessages(xlim(xlims))
        },
        if (!is.null(ylims)){
          suppressMessages(ylim(ylims))
        }
      )
    
    return(plt)
  }
}


## ...
##'
##' 
##' @title density_overlap_2D
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
density_overlap_2D = function(scores_1, scores_2, 
                              name_1, name_2,
                              n_env_cells_per_dim = 100,
                              xlims, ylims,
                              plt_densities = F){
  library(MASS)
  library(dplyr); select = dplyr::select
  
  scores_1_in = scores_1; scores_2_in = scores_2
  scores_1 = scores_1 %>% select(x, y)
  scores_2 = scores_2 %>% select(x, y)
  
  # compute density estimates with gaussian bivariate kernel
  {
    dens_1 = kde2d(x = scores_1[, 1], y = scores_1[, 2], 
                   h = apply(scores_1, 2, bandwidth.nrd), 
                   n = n_env_cells_per_dim,
                   lims = c(xlims, ylims))
    
    dens_2 = kde2d(x = scores_2[, 1], y = scores_2[, 2], 
                   h = apply(scores_2, 2, bandwidth.nrd), 
                   n = n_env_cells_per_dim,
                   lims = c(xlims, ylims))
    
    x = dens_1$x; y = dens_1$y
    z1 = as.vector(t(dens_1$z)); z2 = as.vector(t(dens_2$z))
    
    coords_1 = cbind(as.data.frame(expand_grid(x = x, y = y)), 
                     z = as.vector(z1))
    coords_2 = cbind(as.data.frame(expand_grid(x = x, y = y)), 
                     z = as.vector(z2))
    
    # remove infinitesimal values resulting from kernel smoothing
    {
      #coords_1$z[coords_1$z < (max(coords_1$z, na.rm = T) / 1000)] = 0 
      #coords_2$z[coords_2$z < (max(coords_2$z, na.rm = T) / 1000)] = 0 
    }
    
    # remove infinite values and NAs
    {
      coords_1$z[is.na(coords_1$z) | (is.infinite(coords_1$z))] = 0
      coords_2$z[is.na(coords_2$z) | (is.infinite(coords_2$z))] = 0
    }
  }           
  
  # make raster version to return
  {
    dens_1_norm = dens_1
    dens_1_norm$z = dens_1_norm$z / sum(dens_1_norm$z)
    dens_1_raster = raster::raster(dens_1_norm)
    
    dens_2_norm = dens_2
    dens_2_norm$z = dens_2_norm$z / sum(dens_2_norm$z)
    dens_2_raster = raster::raster(dens_2_norm)
  }
  
  if (!(all(round(coords_1$z, 6) == 0) & 
        all(round(coords_2$z, 6) == 0))){
    
    # plot the env. space and the two subspaces
    {
      if (plt_densities){
        # lc 1
        {
          image(dens_1_norm, 
                main = paste0("Normalized denisity for ",
                              name_1))
          contour(dens_1_norm, add = T)
          
          # plot the original scores to verify the correct position
          # scores_1_sf = st_as_sf(scores_1, coords = c(1,2))
          # plot(scores_1_sf, add = T, col = "red", cex = 0.1, pch = 16)
          # 
          # just to check that the matrix entries really match the coordinates
          # coords_1_sf = st_as_sf(coords_1[which(coords_1$z != 0), ], coords = c(1,2))
          # plot(coords_1_sf, add = T, col = "grey50", cex = 0.1, pch = 16)
        }
        
        # lc 2
        {
          dens_2_norm = dens_2
          dens_2_norm$z = dens_2_norm$z / sum(dens_2_norm$z)
          
          image(dens_2_norm, 
                main = paste0("Normalized denisity for ",
                              name_2))
          contour(dens_2_norm, add = T)
          
          # plot the original scores to verify the correct position
          #scores_2_sf = st_as_sf(scores_2, coords = c(1,2)) %>% st_set_crs(4326)
          #plot(scores_2_sf, add = T, col = "darkgreen", cex = 0.1, pch = 16)
          
          # just to check that the matrix entries really match the coordinates
          #coords_2_sf = st_as_sf(coords_2[which(coords_2$z != 0), ], coords = c(1,2))
          #plot(coords_2_sf, add = T, col = "grey50", cex = 0.1, pch = 16)                        
        }
      }
    }
    
    # compute density overlap
    {
      # rescale by the sum, i.e. make the density 2-point-probability
      coords_1$z = coords_1$z / sum(coords_1$z)
      coords_2$z = coords_2$z / sum(coords_2$z)
      
      D = 1 - (0.5 * (sum(abs(coords_1$z - coords_2$z))))			
      I = 1 - (0.5 * (sqrt(sum((sqrt(coords_1$z) - sqrt(coords_2$z))**2))))
    }
  } else{
    cat("\n\tAll density values for classes", name_1, "and", name_2, " are 0.")
    D = NA
    I = NA
  }
  
  return(list(x = x, y = y, 
              scores_1 = scores_1, scores_2 = scores_2,
              dens_1_raster = dens_1_raster, dens_2_raster = dens_2_raster,
              coords_1 = coords_1, coords_2 = coords_2,
              overlap_D = D, overlap_I = I))
}

## ...
##'
##' 
##' @title plot_distances
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
plot_distances = function(distance_matrix, 
                          classes, class_colors,
                          save_folder=NULL,
                          width=1200, height=800,
                          scale_to_max = T){
  
  library(ggtext)
  
  if (length(grep("LC", colnames(distance_matrix))) == 0){
    colnames(distance_matrix) = paste0("LC", colnames(distance_matrix))
    row.names(distance_matrix) = paste0("LC", row.names(distance_matrix))
  }
  if (length(grep("LC", classes)) == 0){
    classes = paste0("LC", classes)
  }
  
  distance_matrix = as.data.frame(distance_matrix)
  
  distance_matrix$Class_2 = row.names(distance_matrix)
  row.names(distance_matrix) = NULL
  
  dist_table = distance_matrix %>% 
    pivot_longer(cols = starts_with("LC"), names_to = "Class_1")
  colnames(dist_table)[3] = "Distance"
  
  dist_table$Class_1 = factor(dist_table$Class_1, levels=rev(classes))
  dist_table$Class_2 = factor(dist_table$Class_2, levels=rev(classes))
  
  # scale to max for plotting
  if (scale_to_max){ dist_table$Distance = dist_table$Distance / max(dist_table$Distance, na.rm=T) }
  
  plt = ggplot(dist_table, aes(Class_1, Class_2)) + 
    labs(fill="Distance", 
         x="Land cover type", y="Land cover type") + 
    geom_tile(aes(fill = Distance)) + 
    scale_fill_gradient(low="white", high="grey30") +
    geom_text(aes(label = round(Distance,2)), 
              color="indianred3", size=4.5,
              fontface="bold"
    ) +
    theme_minimal() + 
    theme(axis.text.x = element_markdown(angle = 45, hjust = 1,
                                         #fill = rev(class_colors),
                                         #color="white",
                                         size=14,
                                         #face="bold"
    ),
    axis.text.y = element_markdown(angle = 0, hjust = 1,
                                   #fill= rev(class_colors),
                                   #color="white",
                                   size=14,
                                   #face="bold"
    ),
    legend.title = element_text(size=16), 
    legend.text = element_text(size=14),
    axis.title = element_text(size = 16),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 16),
    panel.grid.major = element_line(size = 0.5,
                                    color = "grey20"),
    panel.grid.minor = element_line(size = 0.5,
                                    linetype = 1))
  
  if (!is.null(save_folder)){
    png(save_folder, width=width, height=height)
    print(plt)
    dev.off()
  } else {
    return(plt)
  }
}

## ...
##'
##' 
##' @title distancesAndOverlaps
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
distancesAndOverlaps = function(ID, env_df, class_defs,
                                overlap_metric = "D", n_env_cells_per_dim = 100,
                                plot_niches = F, plot_overlap = F, plot_dists = F, width = 1200, height = 800,
                                save_folder = NULL,
                                verb = 0){
  
  return_list = list()
  
  # find min / max boundaries of all landcover classes in x (=RDA1) and y (=RDA2) direction
  # such that the density estimation fuction that is applied per class has information on how to extend the densities to the env. space of all classes
  {
    coords = env_df %>% select("RDA1", "RDA2")
    xlims = c(min(coords[,1]), max(coords[, 1]))
    ylims = c(min(coords[,2]), max(coords[, 2]))
  }
  
  # class definitions for this site
  {
    class_defs_site = class_defs %>% filter(Class_Code %in% (env_df$Landcover))
    
    # these defintions are used for plotting etc. throughout the script
    class_codes = class_defs_site$Class_Code
    class_names = paste0("LC", class_codes)
    class_labels = paste0("LC",class_defs_site$Class_Code, " ", 
                          gsub(", ", ",\n", class_defs_site$Class_Plotlabel))
    class_cols = class_defs_site$Color_Code
  }
  
  # compute overlap per unique pair of land cover classes
  {
    per_pair = list()
    
    # all combination of landcover classes
    {
      # das kann sicherlich einfacher machen, aber ich wusste nicht wie. Ziel ist es einfach alle Kombinationen der landcover cover Klassen bilden
      # ohne die wo lc1 == lc2 und ohne die Vertauschungen von Kombinationen (als lc2, lc1 nicht aufnehmen, wenn lc1, lc2 schon drin ist)
      pairs = as.data.frame(expand_grid(class_codes, class_codes))
      colnames(pairs) = c("lc1", "lc2")
      pairs = pairs %>% filter(lc1 != lc2)
      pairs_new = data.frame(lc1 = 1, lc2 = 1)
      for (class1 in unique(pairs$lc1)){
        for (class2 in unique(pairs$lc2)){
          if (nrow(pairs_new %>% filter((lc2 == class1) & (lc1 == class2))) == 0){
            pairs_new = rbind(pairs_new, data.frame(lc1 = class1, lc2 = class2))
          }
        }
      }
      pairs = pairs_new[-1,]; pairs = pairs %>% filter(lc1 != lc2); row.names(pairs) = NULL
      
      per_pair[["pairs"]] = pairs
    }
    
    # compute niches and overlaps per pair of land cover classes
    {
      
      if(verb >0){ cat("\n\tCompute density in env. space (niche) per class:") }
      
      per_pair$results = list()
      for (ind in 1:nrow(pairs)){
        
        # prepare pair's results list
        {
          class1 = pairs$lc1[ind]; class2 = pairs$lc2[ind]
          per_pair$results[[ind]] = list()
          per_pair$results[[ind]]$pair = pairs[ind,]
          if(verb >0){ cat("\n\t\tPair: ", paste0(class1, " - ", class2)) }
        }
        
        # extract scores
        {
          scores_1 = env_df[env_df$Landcover == class1, c("cell_ind", "RDA1", "RDA2", "Landcover")]
          scores_2 = env_df[env_df$Landcover == class2, c("cell_ind", "RDA1", "RDA2", "Landcover")]
          
          scores_1_coords = env_df[env_df$Landcover == class1, c("cell_ind", "x", "y", "Landcover")]
          scores_2_coords = env_df[env_df$Landcover == class2, c("cell_ind", "x", "y", "Landcover")]
          scores_coords = rbind(scores_1_coords, scores_2_coords)
        }
        
        # compute niches and overlap
        {
          overlap_obj = density_overlap_2D(scores_1 =  scores_1 %>% rename(x = RDA1, y = RDA2), 
                                           scores_2 = scores_2 %>% rename(x = RDA1, y = RDA2),
                                           name_1 = paste0("LC", class1), 
                                           name_2 = paste0("LC", class2),
                                           n_env_cells_per_dim = n_env_cells_per_dim,
                                           xlims = xlims, ylims = ylims,
                                           plt_densities = F # this gives a plot of the densities of each class too, but later there will be better plots
          )
          per_pair$results[[ind]]$overlap_obj = overlap_obj
          
          x = overlap_obj$x; y = overlap_obj$y
          dens_1 = cbind(overlap_obj$coords_1, Landcover = class1)
          dens_2 = cbind(overlap_obj$coords_2, Landcover = class2)
          
          dens_1_raster = overlap_obj$dens_1_raster
          dens_2_raster = overlap_obj$dens_2_raster
          
          per_pair$results[[ind]]$dens_1_raster = dens_1_raster
          per_pair$results[[ind]]$dens_2_raster = dens_2_raster
        }
        
        # niches data frame and density estimation
        {
          niches_df = rbind(dens_1, dens_2)
          niches_df$Landcover = factor(niches_df$Landcover, levels = sort(unique(niches_df$Landcover)))
          
          scores = rbind(scores_1, scores_2)
          scores = suppressMessages(left_join(scores, scores_coords))
          scores$Landcover = factor(scores$Landcover, levels = sort(unique(niches_df$Landcover)))
          
          per_pair$results[[ind]]$niches = niches_df
          per_pair$results[[ind]]$scores = scores
        }
        
        # overlap, i.e. where to sample the switching cells from = overlap regions in env. space
        {
          overlap_df = data.frame(x = dens_1$x, y = dens_1$y,
                                  z = (dens_1$z * dens_2$z))
          
          # remove infinitesimal values
          #overlap_df$z[overlap_df$z <= max(overlap_df$z)/1000] = 0
          
          # renormalize overlap density to 1
          overlap_df$z = overlap_df$z / sum(overlap_df$z)
          
          # replace NA with zero
          overlap_df$z[which(is.na(overlap_df$z))] = 0
          
          per_pair$results[[ind]]$overlap_df = overlap_df
        }
        
        # overlap as raster
        {
          # is identical to what happens in the next step (overlap_grid_raster)
          # overlap_raster = dens_1_raster * dens_2_raster
          # overlap_raster = overlap_raster / raster::cellStats(overlap_raster, stat = "sum")
          # 
          # per_pair$results[[ind]]$overlap_raster = overlap_raster
          
        }
        
        # match RDA scores and cells on the gridded env. space to plot the overlap on the map
        # e.g. for connecting the density or joint densities in gridded env. space with spatial cooridnates on a map
        # and to later sample points in env. space that switch landcover classes in part 3
        {
          # create grid matrix with overlap density values of the full environmental space
          {
            overlap_grid_matrix = list(x = unique(overlap_df$x),
                                       y = unique(overlap_df$y),
                                       z = matrix(overlap_df$z, nrow = length(unique(overlap_df$y)), 
                                                  ncol = length(unique(overlap_df$x)), byrow = T))
            overlap_grid_raster = raster::raster(overlap_grid_matrix)
            
            per_pair$results[[ind]]$overlap_grid_matrix = overlap_grid_matrix
            per_pair$results[[ind]]$overlap_grid_raster = overlap_grid_raster
            
          }
        }
      }
    }
    
    return_list[["per_pair"]] = per_pair
  }
  
  # extract results (niches and overlap) from different pairs and bind together
  {
    if(verb >0){ cat("\n\n\tExtract data from all pairs.") }
    return_list[["all_pairs"]] = list() 
    return_list$all_pairs[["Pairs"]] = pairs
    
    dens_rasters = lapply(return_list$per_pair$results, 
                          function(u){
                            list(lc1 = u$pair$lc1,
                                 lc2 = u$pair$lc2,
                                 dens_1_raster = u$dens_1_raster,
                                 dens_2_raster = u$dens_2_raster
                            )
                          })
    return_list$all_pairs[["dens_rasters"]] = dens_rasters
    
    niches_df = do.call(rbind, 
                        lapply(return_list$per_pair$results, 
                               function(u){
                                 cbind(u$niches,
                                       lc1 = u$pair$lc1,
                                       lc2 = u$pair$lc2,
                                       D = u$overlap_obj$overlap_D,
                                       I = u$overlap_obj$overlap_I)
                               }))
    niches_df = niches_df %>% mutate(Pair = paste0("LC", lc1, " - LC", lc2))
    return_list$all_pairs[["niches_df"]] = niches_df
    
    scores = do.call(rbind, 
                     lapply(return_list$per_pair$results, 
                            function(u){
                              cbind(u$scores,
                                    lc1 = u$pair$lc1,
                                    lc2 = u$pair$lc2,
                                    D = u$overlap_obj$overlap_D,
                                    I = u$overlap_obj$overlap_I)
                            }))
    scores = scores %>% mutate(Pair = paste0("LC", lc1, " - LC", lc2))
    return_list$all_pairs[["scores_df"]] = scores
    
    overlap_df = do.call(rbind, 
                         lapply(return_list$per_pair$results, 
                                function(u){
                                  cbind(u$overlap_df,
                                        lc1 = u$pair$lc1,
                                        lc2 = u$pair$lc2,
                                        D = u$overlap_obj$overlap_D,
                                        I = u$overlap_obj$overlap_I)
                                }))
    overlap_df = overlap_df %>% mutate(Pair = paste0("LC", lc1, " - LC", lc2))
    return_list$all_pairs[["overlap_df"]] = overlap_df
    
    overlap_grid_matrix = lapply(return_list$per_pair$results,
                                 function(u){
                                   ret_list = u$overlap_grid_matrix
                                 })
    return_list$all_pairs[["overlap_grid_matrix"]] = overlap_grid_matrix
    
    overlap_grid_raster = lapply(return_list$per_pair$results,
                                 function(u){
                                   ret_list = u$overlap_grid_raster
                                 })
    return_list$all_pairs[["overlap_grid_raster"]] = overlap_grid_raster
    
    
    Is = do.call(rbind, 
                 lapply(return_list$per_pair$results, 
                        function(u){ u$overlap_obj$overlap_I }))
    return_list$all_pairs[["Overlap_metrics_I"]] = Is
    
    Ds = do.call(rbind, 
                 lapply(return_list$per_pair$results, 
                        function(u){ u$overlap_obj$overlap_D }))
    return_list$all_pairs[["Overlap_metrics_D"]] = Ds
  }
  
  # create distance matrix from overlap metrics per pair
  {
    if(verb >0){ cat("\n\n\tCreate Distance matrix.") }
    return_list$distance_matrix = list()
    
    dist_mat = matrix(nrow = length(class_codes), ncol = length(class_codes))
    colnames(dist_mat) = row.names(dist_mat) = class_codes
    for (ind in 1:length(return_list$per_pair$results)){
      class1 = as.character(return_list$per_pair$results[[ind]]$pair$lc1)
      class2 = as.character(return_list$per_pair$results[[ind]]$pair$lc2)
      
      # what gets minimized in the optimal transport problem is the product of probability mass, p,  and distance, d.
      # extreme examples:
      # 1. no overlap in env. space means D = 0 (= 0% area overlapping), thus the distance to move the probability mass (share of class A) = d = 1, max distance as max dissimilar classes. 
      # 2. shortest distance possible is a flow class A -> class A (actually no "flow", cells do not change landcover type and stay within their class) 
      #    and this means max similarity -> minimal distance. i.e. dist matrix entry d must be zero 0 (where overlap D is 100% = 1).
      if (overlap_metric == "D"){
        dist_mat[class1, class2] = dist_mat[class2, class1] = 1 - return_list$per_pair$results[[ind]]$overlap_obj$overlap_D   
      } else if (overlap_metric == "I"){
        dist_mat[class1, class2] = dist_mat[class2, class1] = 1 - return_list$per_pair$results[[ind]]$overlap_obj$overlap_I   
      }
      
    }
    
    # each class has d=0 distance, i.e. 100% similarity to itself (D = 1)
    diag(dist_mat) = 0
    
    # normalize the distance matrix by the highest value, such that all sites have distances between 0 and 1
    dist_mat_norm = dist_mat / max(dist_mat, na.rm = T)
    
    return_list$distance_matrix[["dist_mat_orig"]] = dist_mat
    return_list$distance_matrix[["dist_mat_norm"]] = dist_mat_norm
  }
  
  library(metR)
  
  if (plot_dists){
    
    if(verb >0){ cat("\n\n\tPlot distance matrix.") }
    if (!is.null(save_folder)){
      plt_dists = plot_distances(distance_matrix = as.data.frame(dist_mat_norm),
                                 classes = class_codes, 
                                 class_colors = class_cols,
                                 save_folder = file.path(save_folder, paste0("distance_matrix_site_", ID, ".png")),
                                 width = width, height = height,
                                 scale_to_max = T)
    } else {
      plot_distances(distance_matrix = as.data.frame(dist_mat_norm),
                     classes = class_codes, 
                     class_colors = class_cols,
                     save_folder = NULL,
                     width = width, height = height,
                     scale_to_max = T)
    }
    
  }
  
  if (plot_niches){
    
    cat("\n\n\tPlot niches per Pair:")
    niches_df$Landcover = as.numeric(as.character(niches_df$Landcover))
    
    for (ind in 1:nrow(pairs)){
      class1 = pairs[ind, 1]; class2 = pairs[ind, 2]
      class_defs_pair = class_defs_site[which(class_defs_site$Class_Code %in% c(class1, class2)),]
      if(verb >0){  cat("\n\t\tPair:", paste0(class1, " - ", class2)) }
      
      # plot niches seperately as raster
      {
        them = theme(plot.title = element_text(size = 20),
                     plot.subtitle = element_text(size=18),
                     legend.position="right",
                     legend.title = element_text(size=16),
                     legend.text = element_text(size=16),
                     legend.key.size = unit(1,"cm"),
                     axis.text = element_text(size = 16),
                     axis.title = element_text(size = 18),
                     strip.text.x = element_text(size = 16, color = "grey20", face = "bold"),
                     strip.background = element_rect(fill="white"))
        # with points
        plt_1 = ggplot() +
          geom_raster(data = raster::as.data.frame(dens_rasters[[ind]]$dens_1_raster, 
                                                   xy = T), 
                      aes(x = x, y = y, fill = layer)) + 
          scale_fill_gradient(low = "white", high = class_defs_pair$Color_Code[1]) + 
          coord_fixed() + 
          theme_minimal() + 
          labs(title = paste0("Site", ID, " - Density of Class LC", class1),
               fill = "Density") +
          them
        
        
        plt_2 = ggplot() +
          geom_raster(data = raster::as.data.frame(dens_rasters[[ind]]$dens_2_raster, 
                                                   xy = T), 
                      aes(x = x, y = y, fill = layer)) + 
          scale_fill_gradient(low = "white", high = class_defs_pair$Color_Code[2]) + 
          coord_fixed() + 
          theme_minimal() + 
          labs(title = paste0("Site", ID, " - Density of Class LC", class2),
               fill = "Density") + 
          them
        
        # with points added
        plot_df_scores = scores %>% filter((lc1 == class1) & (lc2 == class2)) %>% 
          left_join(env_df %>% select(any_of(c("cell_ind", "Origin_Landcover"))), by = "cell_ind")
        plot_df_scores$Origin_Landcover = factor(plot_df_scores$Origin_Landcover, 
                                                 levels = c("Within PSA",
                                                            "Outside PSA"))
        plt_3 = ggplot() +
          geom_raster(data = raster::as.data.frame(dens_rasters[[ind]]$dens_1_raster, 
                                                   xy = T), 
                      aes(x = x, y = y, fill = layer)) + 
          scale_fill_gradient(low = "white", high = class_defs_pair$Color_Code[1]) + 
          geom_point(data = plot_df_scores %>% filter(Landcover == class1), 
                     mapping = aes(x = RDA1, y = RDA2, colour = Origin_Landcover),
                     size = 1.25, alpha = 0.3) + 
          coord_fixed() + 
          theme_minimal() + 
          labs(title = paste0("Site", ID, " - Density of Class LC", class1),
               fill = "Density",
               Colour = "Origin of Landcover Cell") +
          guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) +
          them
        
        plt_4 = ggplot() +
          geom_raster(data = raster::as.data.frame(dens_rasters[[ind]]$dens_2_raster, 
                                                   xy = T), 
                      aes(x = x, y = y, fill = layer)) + 
          scale_fill_gradient(low = "white", high = class_defs_pair$Color_Code[2]) + 
          geom_point(data = plot_df_scores %>% filter(Landcover == class2), 
                     mapping = aes(x = RDA1, y = RDA2, colour = Origin_Landcover),
                     size = 1.25, alpha = 0.3) + 
          coord_fixed() + 
          theme_minimal() + 
          labs(title = paste0("Site", ID, " - Density of Class LC", class2),
               fill = "Density",
               Colour = "Origin of Landcover Cell") +
          guides(colour = guide_legend(override.aes = list(alpha = 1, size=4))) +
          them
        
        if (!is.null(save_folder)){
          png(filename = file.path(save_folder, "niches",
                                   paste0("niches_classes_", class1, "_", class2, "_site_", ID, ".png")), 
              width = 2*width, height = 2*height)
          gridExtra::grid.arrange(plt_1, plt_2, plt_3, plt_4, nrow = 2)
          dev.off()
        } else{
          gridExtra::grid.arrange(plt_1, plt_2, plt_3, plt_4, nrow = 2)
        }
      }
    }
  }
  
  if (plot_overlap){
    
    if(verb >0){ cat("\n\n\tPlot niche overlap per Pair:") }
    
    # plot overlap in gridded env. space as raster
    for (ind in 1:nrow(pairs)){
      class1 = pairs[ind, 1]; class2 = pairs[ind, 2]
      class_defs_pair = class_defs_site[which(class_defs_site$Class_Code %in% c(class1, class2)),]
      if(verb >0){ cat("\n\t\tPair:", paste0(class1, " - ", class2)) }
      
      plot_df_scores = scores %>% filter((lc1 == class1) & (lc2 == class2)) %>% 
        left_join(env_df %>% select(any_of(c("cell_ind", "Origin_Landcover"))), by = "cell_ind")
      plot_df_scores$Origin_Landcover = factor(plot_df_scores$Origin_Landcover, 
                                               levels = c("Within PSA",
                                                          "Outside PSA"))
      plot_df_scores$Landcover = factor(plot_df_scores$Landcover,
                                        levels = sort(unique(plot_df_scores$Landcover)))
      plot_df = overlap_df %>% filter((lc1 == class1) & (lc2 == class2)) 
      
      them =  theme(plot.title = element_text(size = 20),
                    plot.subtitle = element_text(size=18),
                    legend.position="right",
                    legend.title = element_text(size=16),
                    legend.text = element_text(size=16),
                    legend.key.size = unit(1,"cm"),
                    axis.text = element_text(size = 16),
                    axis.title = element_text(size = 18),
                    strip.text.x = element_text(size = 16, color = "grey20", face = "bold"),
                    strip.background = element_rect(fill="white"))
      
      # without points
      plt_1 = ggplot() +
        geom_raster(data = raster::as.data.frame(overlap_grid_raster[[ind]], xy = T), 
                    aes(x = x, y = y, fill = layer)) + 
        geom_point(data = plot_df, aes(x = x, y = y),
                   size = 0.5, colour = "grey60", alpha = 0.3) +
        scale_fill_gradient(low = "white", high = "darkred") + 
        coord_fixed() +  
        theme_minimal() + 
        labs(title = paste0("Site", ID, " - Overlap of Classes LC", class1, " and LC", class2),
             subtitle = paste0("Overlap, D = ", round(Ds[ind], 4)),
             fill = "Density") +
        them +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size=4)))
      
      # with points of class1
      plt_2 = ggplot() +
        geom_raster(data = raster::as.data.frame(overlap_grid_raster[[ind]], xy = T), 
                    aes(x = x, y = y, fill = layer)) + 
        geom_point(data = plot_df, aes(x = x, y = y),
                   size = 0.5, colour = "grey60", alpha = 0.3) +
        scale_fill_gradient(low = "white", high = "darkred") + 
        geom_point(data = plot_df_scores %>% filter(Landcover == class1), 
                   mapping = aes(x = RDA1, y = RDA2, colour = Origin_Landcover),
                   size = 1.25, alpha = 0.5) + 
        coord_fixed() +  
        theme_minimal() + 
        labs(title = paste0("Site", ID, " - Overlap of Classes LC", class1, " and LC", class2,
                            "\nwith points of LC", class1),
             subtitle = paste0("Overlap, D = ", round(Ds[ind], 4)),
             fill = "Density") +
        them +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size=4)),
               fill = guide_legend(override.aes = list(alpha = 1, size=4)))
      
      # with points of class2
      plt_3 = ggplot() +
        geom_raster(data = raster::as.data.frame(overlap_grid_raster[[ind]], xy = T), 
                    aes(x = x, y = y, fill = layer)) + 
        geom_point(data = plot_df, aes(x = x, y = y),
                   size = 0.5, colour = "grey60", alpha = 0.3) +
        scale_fill_gradient(low = "white", high = "darkred") + 
        geom_point(data = plot_df_scores %>% filter(Landcover == class2), 
                   mapping = aes(x = RDA1, y = RDA2, colour = Origin_Landcover),
                   size = 1.25, alpha = 0.5) + 
        coord_fixed() +  
        theme_minimal() + 
        labs(title = paste0("Site", ID, " - Overlap of Classes LC", class1, " and LC", class2,
                            "\nwith points of LC", class2),
             subtitle = paste0("Overlap, D = ", round(Ds[ind], 4)),
             fill = "Density") +
        them +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size=4)),
               fill = guide_legend(override.aes = list(alpha = 1, size=4)))
      
      if (!is.null(save_folder)){
        png(filename = file.path(save_folder, "overlap",
                                 paste0("overlap_classes_", class1, "_", class2, "_site_", ID, ".png")), 
            width = 2*width, height = 2*height)
        gridExtra::grid.arrange(plt_1, plt_2, plt_3, nrow = 3)
        dev.off()
      } else{
        gridExtra::grid.arrange(plt_1, plt_2, plt_3, nrow = 3)
      }
      
    }
  }
  
  if (!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, "distancesAndOverlaps.rda"))
  } else{
    return(return_list)
  }
}



## ...
##'
##' 
##' @title regGrid
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
regGrid = function(map, cell_length, 
                   save_folder = NULL, width = 1200, height = 800){
  
  library(sf); sf_use_s2(F)
  library(stars)
  
  # land and water cells as polygons
  grid_polys = st_as_sf(st_as_stars(map, dx = cell_length, dy = cell_length))
  grid_polys$cell_ind = 1:nrow(grid_polys)
  
  # land and water cells as points
  grid_points = suppressWarnings(st_centroid(grid_polys))
  
  # land cells as points
  grid_points_land = grid_points %>% filter(values == 1) %>% select(-values)
  
  # land cells as polygons
  grid_polys_land = grid_polys %>% filter(values == 1) %>% select(-values)
  
  # plot water and land cells, only for land cells the land cover will extracted
  if (!is.null(save_folder) & (nrow(grid_polys) <= 10**5)){ 
    png(file.path(save_folder, paste0("land_water_mask_res_", as.integer(cell_length), ".png")), 
        width = width, height = height)
    plot(regMap)
    plot(grid_polys$geometry, add = T, border = "yellow", col = NA)
    plot(grid_polys_land$geometry, add = T, border = "darkgreen", col = NA)
    dev.off()
  } 
  
  return_list = list(grid_polys = grid_polys, grid_points = grid_points,
                     grid_polys_land = grid_polys_land, grid_points_land = grid_points_land)
  
  if(!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, paste0("regGrid_res_", as.integer(cell_length), "_m.rda")))
  } else{
    return(return_list)
  }
}



## ...
##'
##' 
##' @title extractDataGrid
##' @param 
##' save_folder
##' @return ...
##' \item{\code{...}}{...: ...}
##' @importFrom dplyr filter
##' @export
extractDataGrid = function(grid,
                           lc_env_folder, clim_flnames, lc_flname, elev_flname,
                           all_lc_codes, merge_classes, class_defs,
                           plt_aggregated = F, width = 1200, height = 800,
                           save_folder = NULL){
  
  cat("\nExtract land cover and environmental data:")
  # load land cover and environmental data sets
  {
    landcover = read_stars(file.path(lc_env_folder, lc_flname)) %>%
      st_set_crs(4326) %>% setNames("Landcover")
    
    elevation = read_stars(file.path(lc_env_folder, elev_flname)) %>% 
      st_set_crs(4326) %>% setNames("Elevation")
    
    # climate
    {
      climate = lapply(clim_flnames, 
                       function(u) {
                         clim_ind = which(clim_flnames == u) 
                         clim_var = read_stars(file.path(lc_env_folder, u)) %>%
                           st_set_crs(4326) 
                         clim_var = setNames(clim_var, clim_vars[clim_ind])
                         return(clim_var)
                       })
      
      names(climate) = clim_vars
      }
    
    env_vars = c("Elevation", clim_vars)
  }
  
  # project grid to wgs84
  {
    grid = grid %>% st_transform(4326)
  }
  
  geometry_type = as.character(st_geometry_type(grid$geometry[1])) 
  if (geometry_type == "POLYGON"){
    # extract Landcover
    {
      cat("\n\tVariable:", "Landcover")
      
      sample_region = 1:nrow(grid) # sample instead for developing: sample_region = sample(1:nrow(grid), min(nrow(grid), 10))
      
      # extract the landcover compostion in each cell
      lc_region = do.call(rbind,
                          parallel::mclapply(sample_region, 
                                             function(x) {
                                               lc_dist = suppressMessages(
                                                 (landcover[grid[x,]$geometry,,] %>% 
                                                    st_as_stars()) %>% as_tibble() %>% 
                                                   filter(!is.na(Landcover)) %>%
                                                   count(Landcover) %>% 
                                                   mutate(share = n / sum(n))  %>% 
                                                   select(-n)
                                               )
                                               
                                               # merge lc classes and fill missing with 0%
                                               {
                                                 for (merge_class in classes_to_merge){
                                                   class_inds = which(lc_dist$Landcover %in% merge_class$old_classes)
                                                   del_shares = lc_dist$share[class_inds] 
                                                   lc_dist = rbind(lc_dist,
                                                                   data.frame(Landcover = merge_class$new_class,
                                                                              share = sum(del_shares)))
                                                   if (length(class_inds) > 0){ lc_dist = lc_dist[-class_inds,] }
                                                 }
                                                 
                                                 lc_dist = rbind(lc_dist, 
                                                                 data.frame(Landcover = all_lc_codes[!(all_lc_codes %in% lc_dist$Landcover)],
                                                                            share = 0))
                                                 
                                                 lc_dist = lc_dist %>% arrange(Landcover)
                                                 
                                                 if (round(sum(lc_dist$share), 5) != 1){
                                                   stop("Check aggregation of landcover cells. Sum after merging is not equal 1.")
                                                 }
                                               }
                                               
                                               ret_df = t(data.frame(lc_dist$share))
                                               colnames(ret_df) = lc_dist$Landcover; row.names(ret_df) = NULL
                                               
                                               return(ret_df)
                                             }, mc.cores = 1))
      lc_region = data.frame(grid$cell_ind[sample_region],
                             st_coordinates(st_centroid(grid)$geometry)[sample_region,1:2]) %>% 
        cbind(lc_region)
      colnames(lc_region) = c("cell_ind", "x", "y", paste0("LC", all_lc_codes))
      row.names(lc_region) = NULL
      
      region_lc_stats = apply(lc_region %>% select(starts_with("LC")), 2, mean, na.rm = T)
    } # output: lc_region und region_lc_stats
    
    # find dominant landcover type and plot aggregated landcover of the region
    {
      lc_region$Landcover = apply(lc_region%>% select(starts_with("LC")), 1, 
                                  function(u){
                                    aggr_class = as.numeric(gsub("LC", "", names(sort(u, decreasing = T)))[1])
                                  })
      lc_region = st_as_sf(lc_region, coords = c("x", "y"), crs = 4326) %>% st_shift_longitude()
      lc_region$Landcover = factor(lc_region$Landcover,
                                   levels = all_lc_codes[all_lc_codes %in%
                                                           unique(lc_region$Landcover)])
      if (plt_aggregated){
        plt = ggplot() +
          geom_sf(data = regMap %>% st_transform(4326) %>% st_shift_longitude(), fill = NA, col = "grey20", linewidth = 0.4) +
          geom_sf(data = lc_region[sample(1:nrow(lc_region), min(nrow(lc_region), 10000)),],
                  aes(colour = Landcover), 
                  size = 1.5, alpha = 0.7) +
          scale_colour_manual(breaks = all_lc_codes[all_lc_codes %in% levels(lc_region$Landcover)],
                              labels = paste0((class_defs %>% filter(Class_Code %in% levels(lc_region$Landcover)))$Variable_Name, " ",
                                              (class_defs %>% filter(Class_Code %in% levels(lc_region$Landcover)))$Class_Plotlabel),
                              values = (class_defs %>% filter(Class_Code %in% levels(lc_region$Landcover)))$Color_Code) +
          labs(x = "Longitude", y = "Latitude",
               fill = "Landcover Type", title = paste0("Sample with ", min(nrow(lc_region), 5000), 
                                                       " cells of aggregated and merged Land Cover Types")) +
          theme_minimal() + 
          guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
        
        if (!is.null(save_folder)){
          png(file.path(save_folder, paste0("landcover_regGrid_res_", as.integer(cell_length), "_m.png")),
              width = width, height = height)
          print(plt)
          dev.off()
        } else{
          print(plt)
        }
      }
    } # output: save lc_region from here as sf object or the one before?
    
    # extract elevation
    {
      
      cat("\n\tVariable:", "Elevation")
      
      # extract the median elevation value for each grid cell
      elev_region = parallel::mclapply(sample_region,  
                                       function(x) {
                                         suppressMessages(
                                           median((elevation[grid[x,],,] %>% 
                                                     st_as_stars()) %>% as_tibble() %>% 
                                                    filter(!is.na(Elevation)) %>% 
                                                    pull(Elevation), na.rm = T)
                                         )
                                       }, mc.cores = 1) %>% unlist()
      
      elev_region = cbind(data.frame(cell_ind = grid$cell_ind[sample_region],
                                     elev_region))
    }
    
    # extract bio-clim variables
    {
      for (clim_ind in 1:length(climate)){
        
        clim_name = names(climate[[clim_ind]])
        clim_var = climate[[clim_ind]]
        cat("\n\tVariable:", clim_name)
        
        # extract the median elevation value for each grid cell
        var_region = parallel::mclapply(sample_region, 
                                        function(x) {
                                          suppressMessages(
                                            median((clim_var[grid[x,],,] %>% 
                                                      st_as_stars()) %>% as_tibble() %>% 
                                                     filter(!is.na(!!as.name(clim_name))) %>% 
                                                     pull(!!as.name(clim_name)), 
                                                   na.rm = T)
                                          )
                                        }, mc.cores = 1) %>% unlist()
        
        var_region = cbind(data.frame(grid$cell_ind[sample_region],
                                      var_region))
        colnames(var_region) = c("cell_ind", clim_name)
        
        if (clim_ind == 1){
          vars_region = var_region
        } else{
          vars_region = left_join(vars_region, var_region, by = "cell_ind")
        }
      }
    }
    
    # bind all variables together
    {
      coords = st_coordinates(lc_region$geometry)[,1:2]
      df_all = st_as_sf(cbind(cell_ind = lc_region$cell_ind,
                              x = coords[,1],
                              y = coords[,2],
                              Landcover = lc_region$Landcover,
                              Elevation = elev_region  %>% select(-cell_ind), 
                              vars_region %>% select(-cell_ind),
                              geometry = lc_region$geometry
      ), crs = 4326)
      colnames(df_all) = c("cell_ind", "x", "y", "Landcover", env_vars, "geometry")
    }
    
    
    # edit NAs
    {
      # set growing degree days NAs to zero
      df_all = df_all %>% mutate(across(all_of(c("gdd0", "gdd5", "gdd10")), 
                                        ~replace_na(.x, 0)))
      
      # drop remaining rows with any NA-value
      df_all = df_all %>% filter(rowSums(is.na(.)) == 0)
    }
    
    # sort out classes not included in the calibration for the full region
    {
      
      df_all = df_all %>% filter(Landcover %in% include_calib_codes)
      
      lc_stats_final = as.data.frame(
        df_all %>% st_drop_geometry() %>% group_by(Landcover) %>% 
          summarize(Nr_Cells = length(Landcover),
                    Share_Cells =  length(Landcover) / nrow(df_all) * 100)
      )
    }
    
    
    # plot aggregated env. variables
    if (plt_aggregated){
      for (plot_var in env_vars){
        plt = ggplot() +
          geom_sf(data = regMap, fill = NA, col = "grey20", linewidth = 0.4) +
          geom_sf(data = st_as_sf(df_all)[sample(1:nrow(df_all), min(nrow(df_all),10000)),],
                  aes(colour = !!as.name(plot_var)), 
                  size = 2, alpha = 0.8)  +
          labs(x = "Longitude", y = "Latitude",
               fill = plot_var) +
          theme_minimal() + 
          guides(colour = guide_legend(override.aes = list(size=5, alpha = 1)))
        
        
        if (!is.null(save_folder)){
          png(file.path(save_folder, paste0("aggregated_", plot_var, "_regGrid_res_", as.integer(cell_length), "_m.png")),
              width = width, height = height)
          print(plt)
          dev.off()
        } else{
          print(plt)
        }
      }
    }
    
    # correlation between env. variables
    {
      cor_mat = cor(df_all[, env_vars] %>% st_drop_geometry())
      
      library(ggcorrplot)
      row.names(cor_mat) = c("Elevation", clim_labels_breaks)
      colnames(cor_mat) = c("Elevation", clim_labels_breaks)
      plt = ggcorrplot(cor_mat,
                       hc.order = F,
                       type = "upper",
                       lab = T, 
                       colors = c("#6D9EC1", "white", "#E46726")) + 
        labs(title = "Correlation of Environmental Variables in regMap") + 
        theme_plot
      
      if (!is.null(save_folder)){
        png(file.path(save_folder, "correlation_env_vars.png"), width = width, height = height)
        print(plt)
        dev.off()
      } else{
        print(plt)
      }
    }
    
  } else if (geometry_type == "POINT"){
    cat("\n not implemented for point geometries yet.")
    
    # extract exactly the value at each grid point's centre, e.g.:
    # elev_region = suppressMessages(
    #     st_extract(elevation, grid) %>% st_as_stars()
    # )
  } 
  
  return_list = list(lc_env_region = df_all, 
                     lc_stats_extracted = region_lc_stats, lc_stats_dominant = lc_stats_final)
  if (!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, "extractDataGrid.rda"))
  } else{
    return(df_all)
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

