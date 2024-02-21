#################################
### Template: LandCoverChange ###
### 22.10.2023                ###
###  Linux 5/1                ###
#################################

suppressMessages({
library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(readxl)
library(sf); sf_use_s2(FALSE)
library(stars)
library(ggplot2)
library(patchwork)
library(transport)
library(lpSolve)
library(vegan)
})

#############
### Setup ###
#############

### Project folder
project_dir <- "/bioing/data/LandCoverChange/Run_Feb2024"
rda_dir     <- "/bioing/data/LandCoverChange/TestRuns_PSA_Script_RDAtest"
data_dir    <- "/bioing/data/LandCoverChange/LandCoverChangeProject_data_Run_Feb2024"


### LogFile
{
  logfile <- tibble(Name = list.files(rda_dir)) %>% mutate(PSA = grepl("PSA", Name, fixed = TRUE)) %>%
    filter(PSA) %>% mutate(ID = sapply(strsplit(Name, "_"), function(x) as.numeric(x[[2]]))) %>%
    dplyr::select(-PSA) %>%
    left_join(tibble(Name = list.files(project_dir)) %>% 
    mutate(changeMap = grepl("PSA", Name, fixed = TRUE)) %>%
    filter(changeMap) %>% mutate(ID = sapply(strsplit(Name, "_"), function(x) as.numeric(x[[2]]))) %>%
    dplyr::select(ID, changeMap), by = join_by(ID==ID)) %>%
    filter(is.na(changeMap))
}

### Parallel loop
invisible(parallel::mcmapply(1:nrow(logfile), function(run) {

  ID   <- logfile$ID[run]
  name <- logfile$Name[run]
  dir_out <- glue::glue("{project_dir}/{name}")

  ### Create folder and copy files
  {
    dir.create(dir_out)
    invisible(lapply(c("psaChange", "summaryResults"), function(x) dir.create(glue::glue("{dir_out}/{x}"))))
    
    ## Move files
    fls <- list.files(glue::glue("{rda_dir}/{logfile$Name[run]}"), recursive = T)
    invisible(lapply(c("map_init.png", "map_init_withHuman", "Map.png", "rda_summary.png"), function(x) {
      if(any(grepl(x, fls, fixed = TRUE))) {
        file.copy(glue::glue("{rda_dir}/{logfile$Name[run]}/{fls[which(grepl(x, fls, fixed = TRUE))]}"),
                  glue::glue("{project_dir}/{logfile$Name[run]}/{fls[which(grepl(x, fls, fixed = TRUE))]}"))
      }
    }))
}

  runOut <- suppressWarnings(suppressMessages(tryCatch({

    #### class definitions
    {
      ## Class defs
      class_defs <- readxl::read_excel(glue::glue("{data_dir}/settings/landcover_class_definitions_mod.xlsx"))

      ## Native resolution = 309.2208 meter, possible options: 300, 3000, 5000, 10000
      resolution   <- 300

      ## Pixel per class for RDA
      pxl_extract_PSA        <- 800    # max. points for RDA to search within the PSA
      pxl_extract_PSA_buffer <- 50     # min. points for RDA to search in buffer rings if they are not found in the PSA

      ## Maximal spatial buffer around PSA
      max_buffer             <- 1000        # maximum buffer for search

      ## Flow cost parameters
      envs  <- c(0, 100) ## no scaling NULL
      dists <- c(0, 100) ## no scaling NULL
      distF <- 0.05       ## no factor 1

      ## Human classes will be transfered with (TRUE) or without (FALSE) scaling in relation to lancover distribution in PSA
      scaleDensity = T
    }

    ###################################
    ### Pixel flow                  ###
    ###################################

    #### 4.2 Pixel Flow
    {
      if(!file.exists(glue::glue("{dir_out}/psaChange/psaPixelFlow.rda"))) {

        load(glue::glue("{rda_dir}/{name}/summaryResults/initRast.rda"))
        load(glue::glue("{rda_dir}/{name}/summaryResults/rdaOut.rda"))
        load(glue::glue("{rda_dir}/{name}/summaryResults/psaLCover.rda"))
        load(glue::glue("{rda_dir}/{name}/psaOvlp/psaOvlp.rda"))
        load(glue::glue("{rda_dir}/{name}/summaryResults/psaLCover.rda"))
        load(glue::glue("{rda_dir}/{name}/psaEnv/env_stars.rda"))

        ### Landcover classes
        lcovs     <- st_dimensions(psaOvlp@densities)$attributes$values
        lcovs_num <- as.numeric(gsub("LC_", "", lcovs))

        ### Flow table
        psaFlow <- psaLCover@landcover_ts %>% filter(Dataset_ID==ID, Age>500) %>%
          group_split(Age) %>%
          Reduce("rbind", .)

        ### stars as sf & environment
        init_crds <- initRast[1] %>% st_as_sf(na.rm = F) %>% st_centroid() %>%
          rowid_to_column(var = 'pxlID') %>% relocate("pxlID", .before = "Landcover") %>%
          filter(Landcover %in% as.numeric(names(psaFlow)[-c(1:2)])) %>% suppressWarnings()

        init_env <- st_extract(env_stars, init_crds) %>% st_as_sf() %>% st_drop_geometry() %>%
          apply(., 2, function(x) { x[is.na(x)] <- median(x, na.rm = T); x }) %>% suppressMessages()

        init_rda <- as.data.frame(
          predict(
            psaOvlp@rda_mod[[1]],
            init_env %>% as.data.frame(),
            type = "wa",
            scaling = 1
          )) %>% dplyr::select("RDA1", "RDA2") %>% st_as_sf(coords = c("RDA1", "RDA2"))

        # ggplot() +
        #   geom_sf(data = init_rda) +
        #   geom_sf(data = psaOvlp@densities %>% st_bbox() %>% st_as_sfc(), fill = NA)

        init_ovlp <- init_rda %>%
          st_extract(lapply(1:dim(psaOvlp@densities)[3], function(x) {
            split(psaOvlp@densities)[x,] %>% setNames("densities")
          }) %>% do.call("c", .) %>% setNames(lcovs) %>% merge(), .) %>% st_as_sf() %>% st_drop_geometry()

        distCentre <- lapply(1:dim(psaOvlp@densities)[3], function(x) {
          tmp <- split(psaOvlp@densities)[x]
          min <- tmp %>% st_as_sf() %>% pull(1) %>% which.max()
          tmp %>% st_as_sf() %>% slice(min) %>% st_centroid() %>% st_geometry() %>% suppressWarnings()
        })

        init_dist <- distCentre %>% lapply(., function(x) {
          as.numeric(st_distance(init_rda, x))
        }) %>% Reduce("cbind", .) %>% as_tibble() %>% setNames(names(init_ovlp)) %>% suppressWarnings()

        ## costs (with scaling)
        if(!is.null(envs)) {
          envs_scale <- max(envs) - scales::rescale(init_ovlp %>% as.matrix(), envs)
        }  else envs_scale <- max(init_ovlp %>% as.matrix(), na.rm = T) - init_ovlp %>% as.matrix()
        if(!is.null(dists)) {
          dist_scale <- scales::rescale(init_dist %>% as.matrix(), dists)
        }  else dist_scale <- init_dist %>% as.matrix()
        costs <- abind::abind(envs_scale, dist_scale*distF, along = 3) %>%
          apply(., 1:2, sum, na.rm = T)

        classFlow <- matrix(nrow = nrow(init_env), ncol = length(psaFlow$Age %>% unique())+1)
        classFlow[,1] <- init_crds$Landcover

        #### start FLOW ####
        for(y in 2:ncol(classFlow)) {

          # cat("\b\b\b\b\b\b")
          # cat(sprintf("%6d",y))
          # flush.console()

          old_p <- ((tibble(lcov = classFlow[,y-1]) %>% group_by(lcov) %>% summarize(p = n()) %>%
                       right_join(tibble(lcov = lcovs_num), by = "lcov", ) %>% mutate(p = ifelse(is.na(p), 0, p)) %>%
                       arrange(lcov) %>% pull(p))/nrow(classFlow))*100
          new_p <- as.numeric(unlist(psaFlow %>% filter(Age == (psaFlow$Age %>% unique())[y-1]))[-c(1,2)])

          costM   <- 1 - psaOvlp@overlapp[,,2]

          flow    <- transport(old_p, new_p, costM, fullreturn = FALSE) %>% as_tibble() %>%
            mutate(from = lcovs_num[from], to = lcovs_num[to]) %>% filter(from!=to) %>%
            mutate(ident = apply(., 1, function(x) glue::glue("{min(as.numeric(x[1:2]))}_{max(as.numeric(x[1:2]))}"))) %>%
            group_split(ident) %>% lapply(., function(x) {
              if(nrow(x)>1) {
                diff <- x %>% pull(mass) %>% diff() %>% abs()
                x[which.max(x %>% pull(mass)),] %>% mutate(mass = diff)
              } else x
            }) %>% Reduce("rbind",.) %>% dplyr::select(-ident) %>%
            mutate(nrPxl = floor(nrow(classFlow)*mass/100)) %>% arrange(from) %>% group_split(from)

          index_change <- parallel::mclapply(flow, function(x) {

            if(nrow(x)==1) {

              tibble(index = 1:nrow(costs)) %>% bind_cols(costs) %>% filter(classFlow[,y-1]==as.numeric(x[1])) %>%
                dplyr::select('index', lcovs[which(lcovs_num==as.numeric(x[1]))], lcovs[which(lcovs_num==as.numeric(x[2]))]) %>%
                setNames(c("index", "from", "to")) %>% arrange(to) %>% slice(1:as.numeric(x[4])) %>% mutate(newClass = as.numeric(x[2])) %>%
                dplyr::select(index, newClass)

            } else {

              cost_sub <- tibble(index = 1:nrow(costs)) %>% bind_cols(costs) %>% filter(classFlow[,y-1]==unique(x %>% pull(from))) %>%
                dplyr::select(index, lcovs[lcovs_num%in%(x %>% pull(to))]) %>%
                arrange(apply(.[,-1], 1, sum)) %>% slice(1:sum(c(x %>% pull(nrPxl))))


              res <- lpSolve::lp.transport(
                cost.mat  = as.matrix(cost_sub[,-1]),
                direction = "min",
                row.signs = rep("==", nrow(cost_sub)),
                row.rhs   = rep(1, nrow(cost_sub)),
                col.signs = rep("==", ncol(cost_sub[,-1])),
                col.rhs   = c(x %>% pull(nrPxl)))

              cost_sub %>% dplyr::select(index) %>%
                mutate(newClass = lcovs_num[match(names(cost_sub[,-1]), lcovs)][apply(res[["solution"]], 1, function(y) which(y == 1))])
            }

          }, mc.cores = length(flow)) %>% Reduce("rbind", .)

          classFlow[,y] <- classFlow[,y-1]
          classFlow[index_change %>% pull(index),y] <- index_change %>% pull(newClass)

        }
        ### end   FLOW ####
          
        ## test plot
        test_out <- lapply(1:ncol(classFlow), function(x) {
          (as_tibble(classFlow)[,x]) %>%
            table(.) %>% as_tibble() %>% mutate(P = (n/nrow(classFlow))*100) %>%
            dplyr::select(1, 3) %>% setNames(c("LCov","P")) %>%
            full_join(tibble(LCov = as.character(lcovs_num)), by = "LCov") %>%
            arrange(LCov) %>% mutate(P = ifelse(is.na(P), 0, P)) %>%
            mutate(Age = as.numeric(c(-1000, psaFlow$Age)[x]), Type = "pred") %>%
            bind_rows(psaFlow %>% filter(Age == as.numeric(c(-1000, psaFlow$Age)[x])) %>%
                        pivot_longer(-c(Age, Dataset_ID)) %>%
                        dplyr::select(name, value, Age) %>%
                        mutate(type = "orig") %>% setNames(c("LCov", "P", "Age", "Type")))
        }) %>% Reduce("rbind",.)

        pl <- ggplot(test_out, aes(x = P, y = Age, color = Type)) +
          geom_path() +
          facet_wrap(~LCov)
        ggsave(glue::glue("{dir_out}/summaryResults/flow_comparison.png"), plot = pl, width = 30, height = 22, units = "cm")


        psaPixelFlow <- list(crds = init_crds %>% dplyr::select(-Landcover),
                             flow = classFlow %>% as_tibble() %>% setNames(glue::glue("LCov_{c(-1000, psaFlow$Age)}")))
        save(psaPixelFlow, file = glue::glue("{dir_out}/psaChange/psaPixelFlow.rda"))
      } ## end if file.exists
}

    #### 4.3 Output Maps
    {
      
      load(glue::glue("{dir_out}/psaChange/psaPixelFlow.rda"))
      load(glue::glue("{rda_dir}/{name}/summaryResults/initRast.rda"))


      LcovOut <- lapply(1:ncol(psaPixelFlow$flow), function(x) {
        rast <- initRast[1] %>% setNames(names(psaPixelFlow$flow)[x])
        rast[[1]][psaPixelFlow$crds$pxlID] <- as.numeric(unlist(psaPixelFlow$flow[,x]))
        rast
      }) %>% do.call("c", .) %>% merge() %>% setNames("LandcoverChange")

      changeMap    <- class_defs %>% mutate(lcov = ifelse(is.na(Merge_To_Class), Class_Code, Merge_To_Class)) %>%
        dplyr::select(Class_Code, lcov, Class_Plotlabel, Predicted_In_Past, Color_Code)

      save(changeMap, file = glue::glue("{dir_out}/psaChange/changeMap.rda"))
    }
    
    #### Figures
    {
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
      
      ggsave(glue::glue("{project_dir}/{name}/summaryResults/changeMap.png"), 
             plot = maps, width = 6*9, height = 6 + ceiling(dim/4)*6, units = "cm", limitsize = FALSE, bg = "white")
      
      
      dir.create(glue::glue("{project_dir}/{name}/tmp"))
      
      invisible(lapply(1:dim, function(x) {
        
        stars_tmp <- LcovOut[,,,x]
        stars_tmp[[1]][1:length(mergeTab %>% 
                                  filter(!is.na(Color_Code), !duplicated(lcov)) %>% 
                                  pull(lcov))] <- mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov)    
        
        anim <- ggplot() +
          geom_stars(data = stars_tmp, mapping = aes(fill = as.factor(LandcoverChange))) +
          scale_fill_manual(values = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Color_Code),
                            breaks = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(lcov),
                            labels = mergeTab %>% filter(!is.na(Color_Code), !duplicated(lcov)) %>% pull(Class_Plotlabel), name = "Landcover") +
          ggtitle(glue::glue("PSA {ID}: {sapply(strsplit(st_get_dimension_values(LcovOut, 'attributes'), '_'), function(y) y[[2]])[x]}")) +
          theme_void()
        
        ggsave(glue::glue("{project_dir}/{name}/tmp/anim_{sprintf('%03d', x)}.png"), plot = anim,
               width = 19, height = 12, units = "cm", limitsize = FALSE, bg = "white")
      }))
      
    }
    
  }, error = function(e) e)))

  if(any(class(runOut)%in%"error")) {
    png(glue::glue("{dir_out}/summaryResults/error.png"))
    plot(1:10, 1:10, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
    text(5,5,glue::glue("Error Message: \n {runOut$message}"), cex = 2)
    dev.off()
  }
  
}, mc.cores = 50))