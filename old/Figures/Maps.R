library(sf)
sf_use_s2(FALSE)
library(stars)
library(tidyverse)

prj_dir  <- "/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024"
rda_dir  <- "/Volumes/projects/bioing/data/LandCoverChange/TestRuns_PSA_Script_RDAtest"
data_dir <- "/Volumes/projects/bioing/data/LandCoverChange/LandCoverChangeProject_data_Run_Feb2024"

fls     <- list.files("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/", recursive = T, 
                  pattern = "changeMap.rda", full.names = T)

flsTab  <- tibble(ID = as.numeric(sapply(strsplit(fls, "/"), function(x) strsplit(x[[9]], "_")[[1]][[2]])),
                  folder = sapply(strsplit(fls, "/"), function(x) x[[9]]),
                  path = fls)


class_defs <- readxl::read_excel(glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx")) 


for(i in 2:nrow(flsTab)) {
  
  load(glue::glue("{prj_dir}/{flsTab$folder[i]}/psaChange/psaPixelFlow.rda"))
  load(glue::glue("{rda_dir}/{flsTab$folder[i]}/summaryResults/initRast.rda"))
  
  LcovOut <- lapply(1:ncol(psaPixelFlow$flow), function(x) {
    rast <- initRast[1] %>% setNames(names(psaPixelFlow$flow)[x])
    rast[[1]][psaPixelFlow$crds$pxlID] <- as.numeric(unlist(psaPixelFlow$flow[,x]))
    rast
  }) %>% do.call("c", .) %>% merge() %>% setNames("LandcoverChange")

  mergeTab    <- class_defs %>% mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
      dplyr::select(Class_Code, lcov, Class_Plotlabel, Predicted_In_Past, Color_Code, Transfer)

  dim <- as.numeric(dim(LcovOut)[3])
  
  maps <- ggplot() +
    geom_stars(data = LcovOut, mapping = aes(fill = as.factor(LandcoverChange))) +
    facet_wrap(~attributes, ncol = 6, nrow= ceiling(dim/6)) +
    scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                      breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                      labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
    theme_void() +
    theme(legend.position = "top")
  
  ggsave(glue::glue("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/{flsTab$folder[i]}/summaryResults/changeMap.png"), 
         plot = maps, width = 6*8, height = 6 + ceiling(dim/4)*6, units = "cm", limitsize = FALSE)
  
  
  ### Animation
  library(magrittr)
  library(magick)
  
  dir.create(glue::glue("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/{flsTab$folder[i]}/tmp"))
 
   invisible(parallel::mclapply(1:dim, function(x) {
    
    stars_tmp <- LcovOut[,,,x]
    stars_tmp[[1]][1:length(mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov))] <- mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov)    
    
    anim <- ggplot() +
      geom_stars(data = stars_tmp, mapping = aes(fill = as.factor(LandcoverChange))) +
      scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                        breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                        labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
      ggtitle(glue::glue("PSA {flsTab$ID[i]}: {sapply(strsplit(st_get_dimension_values(LcovOut, 'attributes'), '_'), function(y) y[[2]])[x]}")) +
      theme_void()
    
    ggsave(glue::glue("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/{flsTab$folder[i]}/tmp/anim_{sprintf('%03d', x)}.png"), plot = anim,
            width = 19, height = 12, units = "cm", limitsize = FALSE)
   }, mc.cores = 5))
  
  list.files(path=glue::glue("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/{flsTab$folder[i]}/tmp"), pattern = '*.png', full.names = TRUE) %>% 
    image_read() %>%
    image_join() %>%
    image_animate(fps = 4) %>%
    image_write(glue::glue("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/{flsTab$folder[i]}/summaryResults/anim.gif"))
  
  unlink(glue::glue("/Volumes/projects/bioing/data/LandCoverChange/Run_Feb2024/{flsTab$folder[i]}/tmp"), recursive = TRUE)  
}

