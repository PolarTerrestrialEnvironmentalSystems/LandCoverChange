library(tidyverse)

eval_dir    <- "/Volumes/projects/bioing/data/LandCoverChange/Evaluation/"
meta_dir    <- "/Volumes/projects/bioing/data/LandCoverChange/LandCoverChange_Data/data/psa_metadata"
data_dir    <- "/Volumes/projects/bioing/data/LandCoverChange/LandCoverChange_Data"
project_dir <- "/Volumes/projects/bioing/data/LandCoverChange/LandCoverChange_Preprocessed"

psa_metadata  <- read_csv(glue::glue("{meta_dir}/PSA_locations_northern_hemisphere.csv"), show_col_types = F)


preprocess_files <- tibble(path = list.files(project_dir, recursive = T))


output_tab    <- psa_metadata %>% dplyr::select(Dataset_ID, Longitude, Latitude, Continent, Pollen_Source_Radius) %>%
  mutate(Dataset_ID = glue::glue("PSA_{Dataset_ID}")) %>%
  left_join(preprocess_files %>% 
              mutate(Dataset_ID = sapply(strsplit(path, "_"), function(x) paste(x[[1]], x[[2]], sep = "_")),
                     file       = sapply(strsplit(path, "/"), function(y) unlist(strsplit(y[[length(y)]], "[.]"))[1]),
                     exists     = 1) %>% pivot_wider(id_cols = Dataset_ID, names_from = file, values_from = exists), by = "Dataset_ID")


cols <- c("psaVeg", "modernVegetation", "dataset_rast", "lcovShares", "veg_lc_model", "psaLCover", "psaCrds", "psaEnv", "env_stars", "z_tranTab", "psaOvlp", "rda_summary", "map_init", "initRast", "map_init_withHuman")


write_csv(output_tab %>% dplyr::select(c(names(output_tab)[1:5], cols)), glue::glue("{eval_dir}/PreprocessingEval_22052024.csv"))
