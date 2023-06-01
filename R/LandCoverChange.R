#' Simulation of Landcover change over time, based on sedimentary pollen record time series
#'
#' ....
#'
#' @name LandCoverChange-package
#' @docType package
#' @author Peter Ewald and Simeon Lisovski.
NULL


setClass(
  "psaVeg",
   slots = c(regMap     = "sfc_MULTIPOLYGON",
             calibMap   = "sf",
             psas       = "sf",
             vegetation = "data.frame",
             vegetation_interp = "data.frame")
)

setMethod("plot", "psaVeg", function(x) {
  
  ps_poly <- psaVeg@psas %>% left_join(psaVeg@vegetation %>% group_by(Dataset_ID) %>% summarise(oldestAge = max(Age)), by = 'Dataset_ID')
  
  ggplot() +
    geom_sf(psaVeg@calibMap, mapping = aes(geometry = geometry, fill = groups)) +
    scale_fill_brewer(palette = "Set3", name = "Calibration groups") +
    ggnewscale::new_scale_fill() +
    geom_sf(psaVeg@regMap, mapping = aes(geometry = geometry), fill = adjustcolor("white", alpha.f = 0.3)) +
    geom_sf(psaVeg@regMap %>% st_bbox() %>% st_as_sfc() , mapping = aes(geometry = geometry), linewidth = 0.6, fill = NA) +
    geom_sf(ps_poly, mapping = aes(geometry = geometry), alpha = 0.2) +
    # scale_fill_gradientn(colours = c("darkred", "aliceblue")) +
    theme_light()
})


## Data import
##'
##' The function \code{make_psaVeg} loads all vegetation cover files from the specified folder, builds pollen source areas (PSA), i.e. circular polygons around the site's centres with its pollen source radius. 
##' The vegetation cover is filtered to all PSAs within the two maps provided (for computing past land cover and for the calibration of the model for vegetation cover to land cover translation.)
##' @title make_psaVeg
##' @param veg_folder character, absolute path to the folder where vegetation cover files were copied to
##' @param regMap sfc_MULTIPOLYGON, map for that the land cover change shall be computed
##' @param calibMap sfc_MULTIPOLYGON, map from which sites are selected to build the translation model from vegetation cover to landcover time series
##' @param save_folder character, absolute path to the directory where the results list and if plt == T the plots are saved to
##' @param plt boolean, whether to plot the pollen source areas, the number of samples per site, and the oldest age of the samples per site within the maps passed
##' @return a psaVeg object
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
make_psaVeg <- function(veg_folder, regMap, calibMap, silent = TRUE) {
  
  # packages needed: data.table, sf, dplyr, rnaturalearth
  
  if(sf_use_s2()) sf_use_s2(F)
  return_list = list()
  
  # load vegetation and meta data
  {
    if(!silent) cat("\nLoad meta and vegetation data\n")
    
    # read vegetation cover and meta data
    file_names <- list.files(veg_folder)
   
    veg_cover <- lapply(file_names, 
                        function(x) {
                          veg_cov = data.table::fread(file.path(veg_folder, x), verbose = !silent, showProgress = !silent)
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

  }
  
  
  ## Create map with pollen source areas
  {
    if(!silent) cat("\nCreate map with pollen source areas")
    
    # project sites to target geometry from maps
    if (!(st_crs(regMap) == st_crs(calibMap))){
      stop("The same projection for the region and calibration map must be provided.")
    }
    
    # filter sites in any regMap or calibrationMap
    psas = centres %>% st_transform(st_crs(regMap)) %>%
      filter(c(st_intersects(., st_union(st_union(regMap, calibMap$geometry)), sparse = FALSE))) %>%
      mutate(calibMap = calibMap$groups[unlist(apply(st_intersects(., calibMap, sparse = FALSE), 1, function(x) if(any(x)) which(x) else NA))],
             regMap   = c(st_intersects(., regMap, sparse = FALSE))) %>%
      st_buffer(.$`Pollen_Source_Radius [m]`) %>% relocate(geometry, .after = last_col())
  }

  
  # filter vegetation data based on selected pollen source areas, drop taxa not present, and normalize rows to 1
  {
    
    if(!silent) cat("\nFilter vegetation data within maps")
    
    # extract vegetation cover (Dataset_ID, Age, and taxa percentages)
    veg_cover_filtered = lapply(veg_cover,
                                function(x){
                                  veg_cov_out = x %>% filter(Dataset_ID %in% psas$Dataset_ID)
                                  if(nrow(veg_cov_out)>0) {
                                    veg_cov_out <- veg_cov_out %>% dplyr::select(all_of(c("Dataset_ID",
                                                                                          grep("mean", colnames(x), value = T))))
                                    colnames(veg_cov_out) = gsub(" [mean of cover in %]", "",
                                                                 colnames(veg_cov_out), fixed = T)
                                    colnames(veg_cov_out) = gsub("_mean [yrs BP]", "",
                                                                 colnames(veg_cov_out), fixed = T)
                                    veg_cov_out
                                  } else NULL
                                })
    
    suppressMessages({
    veg_cover_filtered <- veg_cover_filtered[!sapply(veg_cover_filtered, is.null)] %>% 
      Reduce("full_join", .)  %>%
      filter(!is.na(Age)) %>% 
      mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>% as.data.frame()
    })
    
    # reorder columns
    veg_cover_filtered = veg_cover_filtered[,c(meta_cols_veg[meta_cols_veg %in% colnames(veg_cover_filtered)], 
                                               colnames(veg_cover_filtered)[!(colnames(veg_cover_filtered) %in% meta_cols_veg)])] %>%
                         dplyr::select(-which(apply(veg_cover_filtered, 2, sum)==0))
    
    veg_cover_filtered[,-c(1,2)] <- t(apply(veg_cover_filtered[,-c(1,2)],  1, function(u) u/sum(u)))
    
  }
  
  
  new(
    "psaVeg",
    regMap = regMap,
    calibMap = calibMap,
    psas = psas,
    vegetation = veg_cover_filtered
  )
  
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
plotVegSite = function(psaVeg, 
                  ID = NULL, 
                  meta_cols = c("Dataset_ID", "Age"),
                  n_taxa = NULL) {
  
  if (is.null(ID)) {ID <- psaVeg@vegetation$Dataset_ID[1]}

    
    df_plot = suppressMessages(
      full_join(cbind(psaVeg@vegetation        %>% filter(Dataset_ID == ID), 
                      Label = "Original Vegetation Cover"),
                cbind(psaVeg@vegetation_interp %>% filter(Dataset_ID == ID), 
                      Label = "Interpolated Vegetation Cover"))
    ) %>% relocate(Label, .after = Age)
    
  
  # select n_taxa taxa
  cols    = names(sort(apply(df_plot[,-c(1:which(names(df_plot)=="Label"))], 2, mean, na.rm = T), decreasing = T))
  if (is.null(n_taxa)) { n_taxa = length(cols) }
  df_plot = df_plot %>% dplyr::select(names(df_plot)[1:which(names(df_plot)=="Label")], cols[1:n_taxa])
  
  df_plot$Label = factor(df_plot$Label, levels = unique(df_plot$Label))
  
  plt = suppressMessages(
    plotTS(data = df_plot, 
                     cols = cols,
                     age_name = "Age",
                     round_ages = F, 
                     age_steps = dt, lines_var = "Label",
                     color = "indianred", alpha = 0.7,
                     title = "",
                     share_to_percent = F,
                     subtitle = paste0("fc = ", fc, ", dt = ", dt,
                                       " k = ", k))
  )
  
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
interpVegTS = function(vegetation, cluster = NULL, fc = 1/500, dt = 500, k = 5, silent = TRUE){

  if(!silent) cat("\n\nInterpolate past vegetation time series")
  
  if(any(class(cluster)=='cluster')) {
    
    clusterEvalQ(cl, {
      library(parallel)
      library(zoo)
      library(dplyr)
      library(corit)
    })
    
    clusterExport(cl, varlist = list('vegetation', 'fc', 'dt', 'k', 'interpolate_ages_per_site'))
    
    suppressMessages({
    (parallel::parLapply(cl, unique(veg_cover$Dataset_ID),
                                          function(x){
                                            veg_site = vegetation %>% filter(Dataset_ID == x) %>% arrange(Age)
                                            suppressWarnings(
                                              interpolate_ages_per_site(df = veg_site,
                                                                        meta_cols = c("Dataset_ID",
                                                                                      "Age"),
                                                                        fc = fc, dt = dt, k = k))})) %>% Reduce("full_join", .)})
  } else {
    suppressMessages({
      (lapply(unique(veg_cover$Dataset_ID),
                                              function(x){
                                                veg_site = vegetation %>% filter(Dataset_ID == x) %>% arrange(Age)
                                                suppressWarnings(
                                                  interpolate_ages_per_site(df = veg_site,
                                                                            meta_cols = c("Dataset_ID",
                                                                                          "Age"),
                                                                            fc = fc, dt = dt, k = k))})) %>% Reduce("full_join", .)})
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
      cat(paste0("\nShow how well the conditions for probability ", 
                 "distributions (compositions) are met:\n"))
      norm_errors = do.call(rbind,lapply(cv_results$folds, function(u){
        do.call(cbind,u$rescaled)
      }))
      round(apply(norm_errors, 2, mean),4)
    }
    
    if (tune_taxa){
      # keep only predictors with highest st. dev. from tuning 
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
            theme_plot + 
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
      #model_full_set = cv_results$model_full_set
      
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
          translation_df = do.call(rbind, lapply(lc_class_inds,
                                                 function(u){
                                                   df1 = translation_df %>% 
                                                     filter((Landcover_Class == 
                                                               lc_classes[u]) & 
                                                              (Taxon %in% sig_taxa[[u]]))
                                                   df1["Significant"] = T
                                                   df2 = translation_df %>% 
                                                     filter((Landcover_Class == 
                                                               lc_classes[u]) &
                                                              !(Taxon %in% sig_taxa[[u]]))
                                                   df2["Significant"] = F
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
  
  return_list = list(model_full_set = model_full_set, 
                     mean_metric = mean_metric, sd_metric = sd_metric,
                     mae_per_class = mae_per_class, sd_per_class = sd_per_class)
  
  if (!is.null(save_folder)){
    save(return_list, file = file.path(save_folder, "vegLcModel.rda"))
  } else{
    return(return_list)
  }
}



# ////////////////////////////////
# helper functions (not exported)

create_subfolders = function(folder_path, dirs = NULL, default = TRUE, silent = TRUE, indent = 1) {
  
  if(default) {
    x = list("settings" = list(),
             "data" = list(),
             "interim_results" = list(psaMap = list(),
                                      interpVegTS = list(),
                                      modernVegShares = list(),
                                      modernLCShares = list(),
                                      vegLcModel = list(),
                                      vegLcPrediction = list()),
             "results" = list())
  } else {
    if(is.null(dirs) | class(dirs)!="list") {
      stop("Subfolders are only accepted as named lists.")
    } else x = dirs
  }
  
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
      
      if(!silent) cat("\n",paste0(rep("\t", indent, collapse=T)), folder_name)
      
      final_path = file.path(folder_path, gsub(" ", "", folder_name))
      
      if (!file.exists(final_path)){
        dir.create(final_path)
      }
      
      if (length(folders) >= 1){
        
        if (!stop_at_last){
          
          create_subfolders(folder_path = file.path(folder_path, folder_name), dirs = folders, default = FALSE,
                            indent = indent+1)
          
        }
      } 
    }
    
    if (!silent & indent == 1) { cat("\n\n") }
  }
  
}

load_rda = function(file_name){
  # loads an rda file, and returns it
  load(file_name)
  get(ls()[ls() != "file_name"]) 
}

