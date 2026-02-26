# Script to assess the quality of the taxonomic assignment performed by bioinformatic algorithms at the species level
# Species will receive a grade between 0 and 5, with 5 being the best grade, meaning the species has been qualitatively detected, and 0 being the worst.

# Creation by Romane Rozanski (rrozanski@ethz.ch) 13.07.2022 (dd.mm.yyyy)
# Updated on 21.01.2025 by Romane Rozanski

########################################################################################################################################################################################

### Inputs:

# DATA
# df_assign (File 1): data.frame containing the following information for each species:
#           Species name / Level of sequencing match (between 0-1) / Genus name / Diversification rate (DR) /
#           Range map availability (y/n) / Occurrence data availability (y/n)

# Regional_checklist (File 2): table containing all the species present in the specific studied area (ecoregion) with the following information:
#                             Region name / Marker used / Species name / Genus name / Sequencing of the species (y/n) / Diversification rate (DR)

# list_sp_coord (File 3): list in which each element corresponds to a species and the coordinates (latitude/longitude start and end) of each eDNA filter where it was detected.
#                        If the sampling station is a point and not a transect, the "latitude/longitude end" has the same values as the "latitude/longitude start."

# Neighbor_ecoreg_checklist: an OPTIONAL table containing a list of species present in neighboring ecoregions

# path_range_map: path to the folder containing all the range maps available for the species (File 4)
# path_occurrences: path to the folder containing all the occurrence data available for the species (File 5)


# PARAMETERS
# steps: which steps should be performed in the pipeline (default: all steps c(1:5))

# --------------------- Step 1
# Neighbor_ecoreg: if the user provides a table of neighboring ecoregions, this checklist will check if the species is present in one of the neighboring ones (default = FALSE)


# --------------------- Step 2
# min_threshold_match_sp: Minimum threshold corresponding to the percentage of match between a DNA sequence and the reference sequence from the database, allowing the classification at the species level (default = 97% = 0.97)
# max_threshold_match_sp: Maximum threshold corresponding to the percentage of match between a DNA sequence and the reference sequence from the database, allowing the classification at the species level (default = 100% = 1)
# fitting_sar_function_match_seq: the chosen SAR function to fit the data and give a score between 0-1 for the percentage of sequence match
#                                 The possible functions are: "power", "powerR", "epm1", "epm2", "p1", "p2", "loga", "koba", "monod", "negexpo", "chapman", "weibull3", "asymp", "ratio", "gompertz", "weibull4", "betap", "logistic", "heleg", "linear"
#                                 (default = "linear")
#                                 See SARS package (https://cran.r-project.org/web/packages/sars/sars.pdf) and reference https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.04271 for fitting functions


# --------------------- Step 3
#no parameters


# --------------------- Step 4
# fitting_sar_function_div_rate: the chosen SAR function to fit the data and give a score between 0-1 for the diversification rate (see "fitting_sar_function_match_seq" for details). (Default = "loga")


# --------------------- Step 5
# format_range_map: format of the range map
# project_crs: projection to compute the geographical coordinates and spatial distances between the filter coordinates and the known range map of a species (default = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
# dist_categ_meters: vector containing categories of distance (in meters) to consider to give a score between 0-1 for the distances between the filter coordinates and the known range map/occurrence of a species.
#                    If dist > the last element of the vector, the score will be 0 (too far).
#                    For range maps: if at least one transect overlaps with the range map, the score will be 1.
#                    All distances in between will have a score between 0-1 (both excluded)
# a: value of the major (equatorial) radius of the ellipsoid (default: 6378137, associated with the WGS84 projection). Needed to estimate the shortest distance (m) between two objects
# f: value of the ellipsoid flattening (default: 1/298.257223563, associated with the WGS84 projection). Needed to estimate the shortest distance (m) between two objects
# format_occ_coord: format of the occurrence data
# dist_max_accepted_radius_occ: Only for occurrence data: the maximal radius distance between the occurrence point and the transect for which the score will be 1 (default = 1000 m). If distance > 1000 m, see parameter "dist_categ_meters"

########################################################################################################################################################################################




#############################################################################################################################################################################################

# Library imports
library("raster")
library("geosphere")
library("spatstat.geom")
library("sf") 
library("sars")
library("stringr")

# ----------------------------------------------------------------------------
taxo_assignment_quality <- function(df_assign = NULL, list_sp_coord = list_sp_coord, Regional_checklist = Regional_checklist,
                                    steps = c(1:5),
                                    Neighbor_ecoreg = FALSE, Neighbor_ecoreg_checklist = Neighbor_ecoreg_checklist,
                                    min_threshold_match_sp = 0.97, max_threshold_match_sp = 1, fitting_sar_function_match_seq = "linear",
                                    fitting_sar_function_div_rate = "loga",
                                    path_range_map = "data/range_maps", format_range_map = ".asc", project_crs = "+proj=longlatitude +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",
                                    dist_categ_meters = c(0, 2500, 5000, 10000), a = 6378137, f = 1/298.257223563,  # a, f default values to compute distance in meters in the WGS84 projection
                                    path_occurences = "data/Occurrences", format_occ_coord = ".rds", dist_max_accepted_radius_occ = 1000) {
  
  # Creation of an empty score table
  steps_names <- c("Regional_checklist", "Seq_match", "Gap_analysis_Genus", "Div_rate", "Map/Occ")
  score_table <- as.data.frame(matrix(ncol = length(steps) + 1, nrow = length(df_assign$Species), data = NA))
  colnames(score_table) <- c(steps_names[steps], "TOTAL")
  rownames(score_table) <- unique(df_assign$Species)
  
  
  
  # --------------------------------------------------------------------------
  # STEP 1: Check if the detected species are in the regional checklist or not
  
  if (1 %in% steps) {
    message("STEP 1: Presence/absence from a regional checklist")
    
    Sp_absent_Regional_checklist <- c()
    
    for (x in (1:nrow(score_table))) {
      
      # Species present in the regional checklist
      if (df_assign$Species[x] %in% Regional_checklist$Species) {
        score_table$Regional_checklist[x] <- 1
      } else {
        
        # Option to choose neighboring ecoregions
        if (Neighbor_ecoreg == TRUE) {
          
          # Species present in a neighboring ecoregion
          if (df_assign$Species[x] %in% Neighbor_ecoreg_checklist$Species) {
            score_table$Regional_checklist[x] <- 0.5 
            
            # Species absent
          } else {
            score_table$Regional_checklist[x] <- 0
            Sp_absent_Regional_checklist <- c(Sp_absent_Regional_checklist, df_assign$Species[x])
          }
          
        } # End if Neighbor_ecoreg == TRUE
        
        else {
          # Species absent
          score_table$Regional_checklist[x] <- 0
          Sp_absent_Regional_checklist <- c(Sp_absent_Regional_checklist, df_assign$Species[x])
        }  
        
      }
      
    } # End of for loop
    
    # For species that are not present in the regional checklist:
    if (length(Sp_absent_Regional_checklist) > 0) {
      message("The ", length(Sp_absent_Regional_checklist), " following detected species are not in the regional checklist: ", paste(Sp_absent_Regional_checklist, collapse = ", "),
              ". Check the output file list_potential_sp to see to which regional species these detections could potentially correspond.")
    }
    
    # Creating a list of the species from the same Genus that are present in the regional checklist and could correspond to the 'false detection' to help the user
    list_potential_sp <- list()
    
    list_potential_sp <- lapply(1:length(Sp_absent_Regional_checklist), function(x) {
      Genus_sp <- str_split(Sp_absent_Regional_checklist[x], "_")[[1]][1]
      sp_present_Regional_checklist <- data.frame(Potential_sp = Regional_checklist[which(Regional_checklist$Genus == Genus_sp), "Species"])
      
      # Creating a table with the data from GapeDNA to see if the species are sequenced
      list_potential_sp[[x]] <- merge(sp_present_Regional_checklist, Regional_checklist[, c("Species", "Sequenced")], by.x = 'Potential_sp', by.y = "Species", all.x = TRUE)
    })
    
    names(list_potential_sp) <- Sp_absent_Regional_checklist
    
  } # End if step 1
  
  
  
  # --------------------------------------------------------------------------
  # STEP 2: Score for sequence matching between the sampling data and the database 
  # (default match between 97-100%): application of a function to fit a SAR model (default: linear model)
  
  if (2 %in% steps) {
    message("STEP 2: Sequence matching")
    
    y = c(0, 1) # Score between 0-1: highest score when match = 100%, lowest score when match = min threshold
    x = c(min_threshold_match_sp, max_threshold_match_sp) # Species level attribution if seq_match is within the threshold (default = 97% = 0.97 and 100% = 1)
    tab_seq_match <- data.frame(x = x, y = y)
    
    # Fit a function to get a score between 0-1 based on SAR models
    sar_model <- paste0("sar_", fitting_sar_function_match_seq)
    sar_result <- do.call(sar_model, args = list(tab_seq_match))
    sar_expression <- as.character(sar_result$model$exp) # Compute the formula of the SAR model
    sar_coeff <- sar_result$model$parNames
    
    # Compute the value of each coefficient in the formula
    for (i in (1:length(sar_coeff))) {
      sar_expression <- gsub(sar_coeff[i], sar_result$par[i], sar_expression)
    }
    
    # Compute the score between 0-1 based on the formula of the chosen SAR function
    for (i in (1:nrow(score_table))) {
      A = as.numeric(df_assign$Match_seq[i])
      score_table$Seq_match[i] <- round(eval(parse(text = sar_expression)), digit = 2)
    }
    
  } # End if step 2
  
  
  
  # --------------------------------------------------------------------------
  # STEP 3: Score for the completeness of the database based on GapEdna data per Genus (present in the area)
  
  if (3 %in% steps) {
    message("STEP 3: Gap analysis of the regional reference database")
    
    for (i in (1:nrow(df_assign))) {
      
      if (df_assign$Genus[i] %in% Regional_checklist$Genus) {
        Regional_checklist_Genus <- unique(Regional_checklist[which(Regional_checklist$Genus == df_assign$Genus[i]), c('Sequenced', 'Genus', 'Species')]) 
        # We use "unique" in case the same species appears multiple times due to the fusion of two ecoregions
        seq_yes <- sum(Regional_checklist_Genus$Sequenced == 'Yes')
        seq_no <- sum(Regional_checklist_Genus$Sequenced == 'No')
        score_table$Gap_analysis_Genus[i] <- round(seq_yes / (seq_yes + seq_no), digit = 2)
        
      } else {
        score_table$Gap_analysis_Genus[i] <- 0
      }
    } # End of for loop
    
  } # End if step 3
 
  
 
  # --------------------------------------------------------------------------
  # STEP 4: Score for the diversification rate: application of a function to fit a SAR model (default: logarithmic model)
  if (4 %in% steps) {
    message("STEP 4: Diversification rate")
    
    y_DR = c(1, 0) # Score between 1 and 0 because the smallest DR corresponds to the highest score (less diversification means better conservation of the genetic diversity)
    
    # Looking for the max and min regional diversification rate
    max_div_rate_regio = max(Regional_checklist$DR, na.rm = T)
    min_div_rate_regio = min(Regional_checklist$DR, na.rm = T)
    
    x_DR = c(min_div_rate_regio, max_div_rate_regio) 
    tab_DR <- data.frame(x = x_DR, y = y_DR)
    
    # Fit a function to get a score between 0-1 based on SAR models
    sar_model_DR <- paste0("sar_", fitting_sar_function_div_rate)
    sar_result_DR <- do.call(sar_model_DR, args = list(tab_DR))
    sar_expression_DR <- as.character(sar_result_DR$model$exp) # Compute the formula of the SAR model
    sar_coeff_DR <- sar_result_DR$model$parNames
    
    # Compute the value of each coefficient in the formula
    for (i in (1:length(sar_coeff_DR))) {
      sar_expression_DR <- gsub(sar_coeff_DR[i], sar_result_DR$par[i], sar_expression_DR)
    }
    
    # Compute the score between 1-0 based on the formula of the chosen SAR function
    for (i in (1:nrow(score_table))) {
      A = df_assign$DR[i]
      score_table$Div_rate[i] <- round(eval(parse(text = sar_expression_DR)), digit = 2)
    }
    
  } # End if step 4
  
  
  
  # --------------------------------------------------------------------------
  # STEP 5: Score for the distance with range_map or species occurrence data
  if (5 %in% steps) {
    message("STEP 5: Range map/Occurrence data")
    
    ### Part 1: With range map 
    #############################
    
    message("Range map")
    range_map_sp_vec <- df_assign[which(df_assign$Range_map == "yes"), "Species"]
    
    if (length(range_map_sp_vec) != 0) {
      
      for (i in (1:length(range_map_sp_vec))) {
        message(paste0("i=", i, " = ", range_map_sp_vec[i]))
        
        range_map_sp <- raster(paste0(path_range_map, "/", range_map_sp_vec[i], format_range_map)) #: the range map has a value of 1 on pixels where the species is present and 0 where the species is absent
        coord_pres_sp_filter <- list_sp_coord[[range_map_sp_vec[i]]]
        colour <- hcl.colors(6, palette = "heat", alpha = NULL, rev = TRUE, fixup = TRUE)
        # plot(range_map_sp, col = colour)
        
        # TRANSECT: Creation of spatial lines for transects: (if points: same latitude and longitude start/end, does not change anything)
        Test <- lapply(1:nrow(coord_pres_sp_filter), function(x) {
          Line(cbind(c(as.numeric(coord_pres_sp_filter$Long_start[x]), as.numeric(coord_pres_sp_filter$Long_end[x])),
                     c(as.numeric(coord_pres_sp_filter$Lat_start[x]), as.numeric(coord_pres_sp_filter$Lat_end[x]))))
        })
        
        # Test if the spatial lines are overlapping cells in which the species is present
        Common_cells <- as.data.frame(do.call(rbind, lapply(1:length(Test), function(w) {
          SL1 <- SpatialLines(list(Lines(Test[[w]], ID = "1")), proj4string = CRS(project_crs))
          # plot(SL1, col = "blue", add = T, lwd = 5)
          DF <- do.call(rbind, extract(range_map_sp, SL1, cellnumbers = TRUE)) # Extract the common cells between the range map and the filter transects
        })))
        
        if ((colSums(Common_cells, na.rm = TRUE)[2] > 0) && (!is.na(colSums(Common_cells, na.rm = TRUE)[2]))) { # At least one filter overlaps with the range map
          score_table[which(rownames(score_table) == range_map_sp_vec[i]), "Map/Occ"] <- 1
          
        } else { # No overlap between all the filters and the range map: we have to calculate the geographic distance between the transect and the closest range map presence
          
          # Convert the raster to a polygon (only cells with values = 1 are kept)
          polygon_map <- rasterToPolygons(range_map_sp, fun = function(x) { x == 1 }, n = 4, na.rm = TRUE, digits = 5, dissolve = FALSE)
          polygon_map@proj4string <- CRS(project_crs)
          # plot(polygon_map, col = "red", add = T)
          polygon_map_sf <- as(polygon_map, "sf") # Convert polygon to an sf object
          
          # Convert coordinates of each filter as spatial lines
          list_spatial_lines_tr <- lapply(1:length(Test), function(w) {
            SpatialLines(list(Lines(Test[[w]], ID = "1")), proj4string = CRS(project_crs))
          })
          
          # Computing the minimum distance from the nearest pixel to each transect filter
          min_dist_among_all_tr <- min(do.call(rbind, lapply(1:length(list_spatial_lines_tr), function(j) {
            
            SL1 <- list_spatial_lines_tr[[j]]
            
            # Convert spatial lines into sf objects
            SL1_sf <- as(SL1, "sf")
            
            # Computing the nearest polygon from the SL 
            nearest_polygon <- sf::st_nearest_feature(SL1_sf, polygon_map_sf, check_crs = T)
            
            # Searching for the nearest polygon coordinates based on the previous result
            nearest_poly_coord <- polygon_map@polygons[[nearest_polygon]]@Polygons[[1]]@coords  # Coords includes the coordinates of each side of the polygon 
            
            # Computing the geographical distance (in meters) between the coordinates of each side of the nearest polygon and the coordinates of the transect: check if the starting point or ending point of the transect is the nearest
            dist_start_tr_nearest_poly <- do.call(rbind, lapply(1:nrow(unique(nearest_poly_coord)), function(x) {
              DF <- geosphere::distGeo(SL1@lines[[1]]@Lines[[1]]@coords[1,], nearest_poly_coord[x,], a = a, f = f) # Default values for a and f related to WGS84
            }))
            
            dist_end_tr_nearest_poly <- do.call(rbind, lapply(1:nrow(unique(nearest_poly_coord)), function(x) {
              DF <- geosphere::distGeo(SL1@lines[[1]]@Lines[[1]]@coords[2,], nearest_poly_coord[x,], a = a, f = f) # Default values for a and f related to WGS84
            })) 
            
            # Choosing the minimal distance computed between one side of the nearest polygon and the start/end of the transect
            dist_nearest_poly <- min(min(dist_start_tr_nearest_poly), min(dist_end_tr_nearest_poly))
            
          }))) # End of lapply j
          
          # Compute the maximum distance accepted between the polygon and the transect
          max_dist_accepted <- dist_categ_meters[length(dist_categ_meters)]
          
          if (min_dist_among_all_tr > max_dist_accepted) { # If the nearest polygon is further than the accepted distance
            score_table[which(rownames(score_table) == range_map_sp_vec[i]), "Map/Occ"] <- 0
          } else {
            
            # Categorizing the distances to get a value between 0 and <1 (1 excluded because 1 means the transect is inside the range map)
            pas = 1 / length(dist_categ_meters)  # To get the step to equally divide the dist_categ distance vector 
            score_dist = rev(round(seq(pas, (1 - pas), by = pas), digit = 2)) # Because 1 and 0 are excluded & "rev" (reverse) to get the highest value for the closest range map presence
            
            categ_dist_nearest_poly <- cut(min_dist_among_all_tr, dist_categ_meters, include.lowest = TRUE, dig.lab = 5)
            level_dist <- levels(categ_dist_nearest_poly)
            
            for (x in 1:length(level_dist)) {
              if (categ_dist_nearest_poly == level_dist[x]) {
                score_table[which(rownames(score_table) == range_map_sp_vec[i]), "Map/Occ"] <- score_dist[x]
              }
            } # End for x
            
          } # End else score dist
          
        }  # End else "overlap" 
        
      } # End of for i in length species_range_map
      
    } # End if length != 0
    
    
    
    
    ### Part 2: With occurrence data
    ###############################
    message("Occurrence")
    occ_sp_vec <- df_assign[which(df_assign$Occurrence_data == "yes"), "Species"]
    
    if (length(occ_sp_vec) != 0) {
      
      for (i in (1:length(occ_sp_vec))) {
        message(paste0("i=", i, " = ", occ_sp_vec[i]))
        
        occ_coord_sp <- as.data.frame(read.csv(paste0(path_occurences, "/", occ_sp_vec[i], format_occ_coord), row.names = 1))
        coord_pres_sp_filter <- list_sp_coord[[occ_sp_vec[i]]]
        
        # Transform the coordinates of occurrence into spatial points
        pts_occ <- SpatialPoints(coords = cbind(occ_coord_sp$decimallongitude, occ_coord_sp$decimallatitude), proj4string = CRS(project_crs))
        
        # Transform the filter coordinates into transects spatial lines
        Test <- lapply(1:nrow(coord_pres_sp_filter), function(x) {
          Line(cbind(c(as.numeric(coord_pres_sp_filter$Long_start[x]), as.numeric(coord_pres_sp_filter$Long_end[x])),
                     c(as.numeric(coord_pres_sp_filter$lat_start[x]), as.numeric(coord_pres_sp_filter$Lat_end[x]))))
        })
        
        list_spatial_lines_tr <- lapply(1:length(Test), function(w) {
          SpatialLines(list(Lines(Test[[w]], ID = "1")), proj4string = CRS(project_crs))
        })
        
        # Computing the minimum distance from the nearest point to each transect filter
        min_dist_among_all_tr <- min(do.call(rbind, lapply(1:length(list_spatial_lines_tr), function(j) {
          
          SL1 <- list_spatial_lines_tr[[j]]
          # plot(SL1, col = "black", add = T, lwd = 5)
          
          # Convert spatial lines and points into sf objects
          SL1_sf <- as(SL1, "sf")
          pts_sf <- as(pts_occ, "sf")
          
          # Computing the nearest point from the SL 
          nearest_point <- sf::st_nearest_feature(SL1_sf, pts_sf, check_crs = T)
          
          # Searching for the nearest polygon coordinates based on the previous result
          nearest_point_coord <- pts_occ[nearest_point]@coords
          
          # Computing the geographical distance (in meters) between the coordinates of the nearest point and the coordinates of the transect: check if the starting point or ending point of the transect is the nearest
          dist_start_tr_nearest_point <- geosphere::distGeo(SL1@lines[[1]]@Lines[[1]]@coords[1,], nearest_point_coord, a = a, f = f) # Default values for a and f related to WGS84
          
          dist_end_tr_nearest_point <- geosphere::distGeo(SL1@lines[[1]]@Lines[[1]]@coords[2,], nearest_point_coord, a = a, f = f) # Default values for a and f related to WGS84
          
          # Choosing the minimal distance computed between one side of the nearest polygon and the start/end of the transect
          dist_nearest_point <- min(min(dist_start_tr_nearest_point), min(dist_end_tr_nearest_point))
          
        }))) # End of lapply j
        
        # Compute the maximum distance accepted between the point and the transect
        max_dist_accepted <- dist_categ_meters[length(dist_categ_meters)]
        
        if (min_dist_among_all_tr > max_dist_accepted) { # If the nearest point is further than the accepted distance
          score_table[which(rownames(score_table) == occ_sp_vec[i]), "Map/Occ"] <- 0
        } else {
          
          if (min_dist_among_all_tr < dist_max_accepted_radius_occ) { # 1000 is the default value, can be changed 
            score_table[which(rownames(score_table) == occ_sp_vec[i]), "Map/Occ"] <- 1 # If the transect is in a 1000m area around the nearest point, value = 1
          } else {
            
            # Categorizing the distances to get a value between 0 and <1 (1 excluded because 1 means the transect is inside the range map)
            pas = 1 / length(dist_categ_meters)  # To get the step to equally divide the dist_categ distance vector 
            score_dist = rev(round(seq(pas, (1 - pas), by = pas), digit = 2)) # Because 1 and 0 are excluded & "rev" (reverse) to get the highest value for the closest range
            
            categ_dist_nearest_point <- cut(min_dist_among_all_tr, dist_categ_meters, include.lowest = TRUE, dig.lab = 5)
            level_dist <- levels(categ_dist_nearest_point)[1] 
            
            for (x in 1:length(level_dist)) {
              if (categ_dist_nearest_point == level_dist[x]) {
                score_table[which(rownames(score_table) == occ_sp_vec[i]), "Map/Occ"] <- score_dist[x]
              }
            } # End for x
            
          }
          
        } # End else score dist
        
      }  # End of for i in length species_occurrence
      
    } # End if length != 0
    
  } # End if step 5
  
 
  
  #--------------------------------------------------------------------------
  # STEP 6: Calculation of the total score
  message("STEP 6: TOTAL SCORE")
  score_table$TOTAL <- round(rowSums(score_table[, 1:ncol(score_table) - 1], na.rm = F), digit = 2)
  
  score_table <- score_table[order(score_table$TOTAL, decreasing = T),] # Order the table to get the highest score first
  
  if (1 %in% steps) {
    return(list(score_table = score_table, list_potential_sp = list_potential_sp))
  } else {
    return(score_table)
  }

  
    
} # End of the function taxonomic quality
#################################################################################################################################################################
