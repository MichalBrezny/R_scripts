# Function loops through raster list and prints their extent, # of rows and cols and resolution

# Load necessary library
library(terra)

# Function to check and show details of all rasters in the list
show_rasters_details <- function(raster_list) {
  if (length(raster_list) == 0) {
    stop("The raster list is empty.")
  }
  
  details <- list()  # Initialize a list to store details
  
  # Loop through rasters and print details
  for (i in 1:length(raster_list)) {
    current_raster <- rast(raster_list[[i]])
    current_extent <- ext(current_raster)
    current_nrows <- nrow(current_raster)
    current_ncols <- ncol(current_raster)
    current_res <- res(current_raster)
    
    # Store the details in a list
    details[[raster_list[[i]]]] <- list(
      extent = current_extent,
      nrows = current_nrows,
      ncols = current_ncols,
      resolution = current_res
    )
    
    # Print raster details
    cat(sprintf("Details of raster '%s':\n", raster_list[[i]]))
    cat("Extent:\n")
    print(current_extent)
    cat(sprintf("Number of Rows: %d\n", current_nrows))
    cat(sprintf("Number of Columns: %d\n", current_ncols))
    cat(sprintf("Resolution (X, Y): %f, %f\n", current_res[1], current_res[2]))
    cat("\n")
  }
  
  return(details)  # Return a list of details for further use if needed
}

# Example usage
raster_list <- c("C:/path/to/raster1.tif",
                 "C:/path/to/raster2.tif",
                 "C:/path/to/raster3.tif")

# Show details of all rasters
raster_details <- show_rasters_details(raster_list)
