# ---
#
# title: "Calculate Bedding Dip and Dip Direction (Functional Version)"
# author: "Original script by Michal Brezny, heavily refactored using Gemini"
# date: "2026-01-15"
# description: >
#   This script calculates bedding dip and dip direction using a functional approach.
#   It supports two workflows:
#   1. From a point layer: Calculates normal vectors, interpolates them, and then
#      computes dip/dip direction.
#   2. From pre-calculated rasters of normal vector components.
#
#   The method is based on Santangelo, M., et al. (2015). A method for the
#   assessment of the influence of bedding on landslide abundance and types.
#   Landslides, 12(2), 295–309.
#
# ---

# 1. SETUP
# ------------------------------------------------------------------------------
# Load required packages.
# If you don't have them, run: install.packages(c("terra", "gstat"))
library(terra)
library(gstat)

# 2. FUNCTION DEFINITIONS
# ------------------------------------------------------------------------------

#' Calculate Normal Vectors from Dip and Dip Direction
#'
#' @param points SpatVector of points.
#' @param dip_col Character, name of the dip column.
#' @param dip_dir_col Character, name of the dip direction column.
#' @return SpatVector with new columns 'nx', 'ny', 'nz'.
calculate_normal_vectors <- function(points, dip_col, dip_dir_col) {
  message("Calculating normal vectors from dip and dip direction...")
  
  # Ensure the specified columns exist
  if (!all(c(dip_col, dip_dir_col) %in% names(points))) {
    stop("Dip and/or Dip Direction columns not found in the point file.")
  }
  
  # Convert degrees to radians
  dip_rad <- points[[dip_col]] * pi / 180
  dip_dir_rad <- points[[dip_dir_col]] * pi / 180
  
  # Calculate normal vector components (pole to the plane)
  points$nx <- sin(dip_rad) * sin(dip_dir_rad)
  points$ny <- sin(dip_rad) * cos(dip_dir_rad)
  points$nz <- cos(dip_rad)
  
  return(points)
}


#' Interpolate Vector Components onto a Raster Grid
#'
#' @param points SpatVector with nx, ny, nz columns.
#' @param template_raster SpatRaster to use for grid definition.
#' @param method Character, "idw" or "kriging".
#' @param save_variance Logical, if TRUE, saves kriging variance maps.
#' @param output_dir Character, path to save intermediate files.
#' @return A SpatRaster stack with interpolated 'rx', 'ry', 'rz' layers.
interpolate_components <- function(points, template_raster, method = "kriging", save_variance = TRUE, output_dir = "intermediate_rasters") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Create interpolation grid
  grid <- rast(ext(template_raster), resolution = res(template_raster), crs = crs(template_raster))
  
  if (method == "idw") {
    message("Interpolating vector components using Inverse Distance Weighting (IDW)...")
    g_nx <- gstat::gstat(formula = nx ~ 1, data = points)
    g_ny <- gstat::gstat(formula = ny ~ 1, data = points)
    g_nz <- gstat::gstat(formula = nz ~ 1, data = points)
    
    rx <- interpolate(grid, g_nx, ext = ext(template_raster))
    ry <- interpolate(grid, g_ny, ext = ext(template_raster))
    rz <- interpolate(grid, g_nz, ext = ext(template_raster))
    
  } else if (method == "kriging") {
    message("Interpolating vector components using Ordinary Kriging...")
    
    # Define a helper for the kriging process
    krige_component <- function(formula, data, grid, component_name) {
      message("...fitting variogram for ", component_name)
      var <- gstat::variogram(formula, data)
      fit <- gstat::fit.variogram(var, model = vgm("Sph"))
      g <- gstat::gstat(formula = formula, model = fit, data = data)
      
      message("...performing kriging for ", component_name)
      kriged_result <- interpolate(grid, g, ext = ext(template_raster), xy = TRUE)
      
      # Save variance if requested
      if (save_variance) {
        message("...saving variance map for ", component_name)
        writeRaster(kriged_result$var1.var, file.path(output_dir, paste0(component_name, "_variance.tif")), overwrite = TRUE)
      }
      return(kriged_result$var1.pred)
    }
    
    rx <- krige_component(nx ~ 1, points, grid, "rx")
    ry <- krige_component(ny ~ 1, points, grid, "ry")
    rz <- krige_component(nz ~ 1, points, grid, "rz")
    
  } else {
    stop("Invalid interpolation_method. Choose 'idw' or 'kriging'.")
  }
  
  # Stack results
  interpolated_rasters <- c(rx, ry, rz)
  names(interpolated_rasters) <- c("nx", "ny", "nz")
  
  # Save intermediate results
  message("...saving interpolated component rasters")
  writeRaster(interpolated_rasters$nx, file.path(output_dir, "rx_interpolated.tif"), overwrite = TRUE)
  writeRaster(interpolated_rasters$ny, file.path(output_dir, "ry_interpolated.tif"), overwrite = TRUE)
  writeRaster(interpolated_rasters$nz, file.path(output_dir, "rz_interpolated.tif"), overwrite = TRUE)
  
  return(interpolated_rasters)
}


#' Calculate Dip and Dip Direction from Vector Components
#'
#' @param component_rasters SpatRaster stack with 'nx', 'ny', 'nz' layers.
#'   Can optionally include pre-calculated 'r' and 'm' layers.
#' @return A SpatRaster stack with 'dip_direction' and 'dip' layers.
calculate_bedding_geometry <- function(component_rasters) {
  message("Calculating final Dip and Dip Direction rasters...")
  
  # Calculate magnitudes if not provided
  if (!"r" %in% names(component_rasters)) {
    component_rasters$r <- sqrt(component_rasters$nx^2 + component_rasters$ny^2 + component_rasters$nz^2)
  }
  if (!"m" %in% names(component_rasters)) {
    component_rasters$m <- sqrt(component_rasters$nx^2 + component_rasters$ny^2)
  }
  
  # Helper function for Dip Direction
  dipdir_calc_fun <- function(nx, ny, nz, m) {
    z <- rep(NA, length(nx))
    i <- which(nz >= 0 & nx >= 0 & ny >= 0); z[i] <- asin(nx[i] / m[i])
    i <- which(nz >= 0 & nx < 0 & ny >= 0);  z[i] <- 2 * pi - asin(-nx[i] / m[i])
    i <- which(nz >= 0 & nx < 0 & ny < 0);  z[i] <- pi + asin(-nx[i] / m[i])
    i <- which(nz >= 0 & nx >= 0 & ny < 0); z[i] <- pi - asin(nx[i] / m[i])
    i <- which(nz < 0 & nx >= 0 & ny >= 0);  z[i] <- pi + asin(nx[i] / m[i])
    i <- which(nz < 0 & nx < 0 & ny >= 0);  z[i] <- pi - asin(-nx[i] / m[i])
    i <- which(nz < 0 & nx < 0 & ny < 0);  z[i] <- asin(-nx[i] / m[i])
    i <- which(nz < 0 & nx >= 0 & ny < 0);  z[i] <- asin(nx[i] / m[i])
    return(z * 180 / pi)
  }
  
  # Helper function for Dip
  dip_calc_fun <- function(nz, r) {
    z <- rep(NA, length(nz))
    ia <- which(nz >= 0); z[ia] <- acos(nz[ia] / r[ia])
    ib <- which(nz < 0); z[ib] <- pi - acos(-nz[ib] / r[ib])
    return(z * 180 / pi)
  }
  
  # Apply calculations
  dip_direction <- lapp(component_rasters[[c("nx", "ny", "nz", "m")]], fun = dipdir_calc_fun)
  dip <- lapp(component_rasters[[c("nz", "r")]], fun = dip_calc_fun)
  
  # Stack and return final results
  final_rasters <- c(dip_direction, dip)
  names(final_rasters) <- c("dip_direction", "dip")
  return(final_rasters)
}


# 3. MAIN EXECUTION BLOCK
# ------------------------------------------------------------------------------

# --- User-Defined Parameters ---

# Set to TRUE to calculate from points, FALSE to use existing rasters.
calculate_from_points <- TRUE

# Point-based workflow parameters
if (calculate_from_points) {
  input_point_file <- "data/structural_points.gpkg"
  dip_col_name <- "dip"
  dip_dir_col_name <- "dip_direction"
  template_raster_file <- "data/template_raster.tif"
  interpolation_method <- "kriging" # "idw" or "kriging"
  save_kriging_variance <- TRUE
  intermediate_dir <- "intermediate_rasters"
}

# Raster-based workflow parameters
if (!calculate_from_points) {
  raster_dir <- "rasters"
  path_nx <- file.path(raster_dir, "rx.tif")
  path_ny <- file.path(raster_dir, "ry.tif")
  path_nz <- file.path(raster_dir, "rz.tif")
}

# General output parameters
output_dir <- "output"

# --- Workflow Logic ---

main <- function() {
  if (calculate_from_points) {
    # --- Point-based Workflow ---
    if (!file.exists(input_point_file)) stop("Input point file not found.")
    if (!file.exists(template_raster_file)) stop("Template raster file not found.")
    
    # 1. Load data
    points <- vect(input_point_file)
    template <- rast(template_raster_file)
    
    # 2. Calculate normal vectors
    points_with_vectors <- calculate_normal_vectors(points, dip_col_name, dip_dir_col_name)
    
    # 3. Interpolate components
    component_rasters <- interpolate_components(
      points = points_with_vectors,
      template_raster = template,
      method = interpolation_method,
      save_variance = save_kriging_variance,
      output_dir = intermediate_dir
    )
    
  } else {
    # --- Raster-based Workflow ---
    message("Workflow selected: Use pre-existing rasters.")
    input_files <- c(path_nx, path_ny, path_nz)
    if (!all(file.exists(input_files))) {
      stop("One or more input component raster files are missing.")
    }
    component_rasters <- rast(input_files)
    names(component_rasters) <- c("nx", "ny", "nz")
  }
  
  # 4. Calculate final bedding geometry
  final_rasters <- calculate_bedding_geometry(component_rasters)
  
  # 5. Save final output
  if (!dir.exists(output_dir)) dir.create(output_dir)
  message("Saving final rasters to '", output_dir, "' directory...")
  
  writeRaster(final_rasters$dip_direction, file.path(output_dir, "dip_direction.tif"), overwrite = TRUE)
  writeRaster(final_rasters$dip, file.path(output_dir, "dip.tif"), overwrite = TRUE)
  
  message("Script finished successfully.")
}

# Run the main function
main()

# --- END OF SCRIPT ---
