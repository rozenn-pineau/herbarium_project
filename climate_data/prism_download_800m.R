# =============================================================================
# PRISM 800m Climate Data: Download & Clip to Sample Regions
# =============================================================================
# PURPOSE (this stage):
#   Download PRISM 800m grids for all variables and time periods needed by
#   your samples, and clip them to a ~100 km buffer around each sample
#   location. No point extraction is done here — rasters are saved for later
#   re-use once sample coordinates have been revised.
#
# NEXT STAGE (separate script):
#   Once lat/long coordinates are finalised, run extract_climate_values.R
#   (see bottom of this file for a stub) to pull point values from the
#   already-clipped rasters. No re-downloading needed.
#
# OUTPUT
# ------
# data/prism_clipped/
#   One GeoTiff per variable × time step, clipped to the union of all
#   per-sample 100 km buffers. File naming (from PRISM):
#     PRISM_<var>_stable_800mD2_<YYYYMM>_bil_clipped.tif   (monthly)
#     PRISM_<var>_stable_800mD2_<YYYYMMDD>_bil_clipped.tif (daily)
#
# data/sample_aoi.gpkg
#   The AOI polygon (union of buffers) saved as a GeoPackage so you can
#   inspect it in QGIS / ArcGIS and confirm coverage before downloading.
#
# =============================================================================


# --- 0. Packages -------------------------------------------------------------

pkgs <- c("prism", "terra", "sf", "tidyverse", "googlesheets4")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install)

library(prism)
library(terra)
library(sf)
library(tidyverse)
library(googlesheets4)

rm(list = ls())
# --- 1. Configuration --------------------------------------------------------
setwd("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/4.Herbarium/7.Metadata/0.PRISM")

PRISM_RAW_DIR  <- "data/prism_raw"    # temporary; full US rasters land here
PRISM_CLIP_DIR <- "data/prism_clipped" # permanent; clipped rasters saved here
OUTPUT_DIR     <- "data"

DELETE_AFTER_CLIP <- TRUE   # strongly recommended — full US 800m files are large
KEEP_ZIP          <- FALSE
BUFFER_KM         <- 100    # radius around each sample point

VARIABLES <- c("ppt", "tmax", "tmin", "tmean", "tdmean", "vpdmin", "vpdmax")

dir.create(PRISM_RAW_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(PRISM_CLIP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR,     recursive = TRUE, showWarnings = FALSE)

prism_set_dl_dir(PRISM_RAW_DIR)


# --- 2. Load sample data -----------------------------------------------------

samples_raw <- read.table("/Users/rozenn/Library/CloudStorage/GoogleDrive-rozennpineau@uchicago.edu/My Drive/Work/9.Science/4.Herbarium/7.Metadata/herbarium_coordinates.csv", 
                          sep = ",", header = T)
samples <- samples_raw |>
  mutate(
    # Use original coords for AOI construction — revised coords come later
    use_lat = as.numeric(lat),
    use_lon = as.numeric(long),
    year    = as.integer(unlist(year)),
    month   = as.integer(month),
    sample_id = coalesce(
      na_if(as.character(barcode), ""),
      na_if(as.character(accession), ""),
      paste0("row_", row_number())
    )
  ) |>
  filter(
    !is.na(use_lat), !is.na(use_lon),
    !is.na(year), !is.na(month),
    between(use_lat, -90, 90),
    between(use_lon, -180, 180)
  )

message(sprintf("%d samples loaded with valid coordinates and dates.", nrow(samples)))


# --- 3. Build Area Of Interest: union of per-sample 100 km buffers -----------------------
# Buffering in a projected CRS (metres), then back to WGS84 for PRISM

message(sprintf("Building AOI: %d km buffer around each sample...", BUFFER_KM))

samples_sf <- st_as_sf(samples,
                        coords = c("use_lon", "use_lat"),
                        crs    = 4326)

# Project to equal-area CRS for accurate km buffering (Albers CONUS)
samples_proj <- st_transform(samples_sf, crs = 5070)

# Buffer each point individually, take the union, reproject to WGS84
aoi <- samples_proj |>
  st_buffer(dist = BUFFER_KM * 1000) |>   # metres
  st_union() |>
  st_transform(4326)

# Save AOI for inspection in GIS software
st_write(aoi, file.path(OUTPUT_DIR, "sample_aoi.gpkg"), delete_dsn = TRUE)
message("AOI saved to data/sample_aoi.gpkg — inspect this before downloading.")

# Report bounding box
bb <- st_bbox(aoi)
message(sprintf(
  "AOI bounding box: lon [%.2f, %.2f], lat [%.2f, %.2f]",
  bb["xmin"], bb["xmax"], bb["ymin"], bb["ymax"]
))


# --- 4. Helper: clip a PRISM raster to the AOI and save --------------------

clip_and_save <- function(prism_folder, aoi, clip_dir,
                          delete_original = TRUE) {
  rast_file <- list.files(prism_folder,
                          pattern = "\\.(tif|bil)$",
                          full.names = TRUE)[1]
  if (is.na(rast_file)) {
    warning("No raster found in: ", prism_folder)
    return(invisible(NULL))
  }

  r        <- rast(rast_file)
  aoi_proj <- st_transform(aoi, crs(r))
  aoi_vect <- vect(aoi_proj)
  r_clip   <- crop(r, aoi_vect)    # bounding box crop; use mask=TRUE for
                                    # precise polygon mask (slower, smaller)

  out_name <- paste0(basename(prism_folder), "_clipped.tif")
  out_path <- file.path(clip_dir, out_name)
  writeRaster(r_clip, out_path, overwrite = TRUE,
              gdal = c("COMPRESS=LZW", "TILED=YES"))

  message("  Saved: ", out_name)

  if (delete_original) {
    unlink(prism_folder, recursive = TRUE)
    message("  Deleted full US raster: ", basename(prism_folder))
  }

  return(invisible(out_path))
}


# --- 5. Identify required time periods from sample metadata -----------------

monthly_years  <- sort(unique(samples$year[samples$year>1895]))  #monthly PRISM starts 1895
daily_samples  <- samples |> filter(year >= 1981)   # daily PRISM starts 1981
daily_ym       <- daily_samples |> distinct(year, month) |> arrange(year, month)

message(sprintf(
  "\nCoverage needed:\n  Monthly: %d unique years (%d–%d)\n  Daily:   %d year×month combos (%d samples pre-1981 skipped)",
  length(monthly_years),
  min(monthly_years), max(monthly_years),
  nrow(daily_ym),
  sum(samples$year < 1981)
))


# --- 6. Download and clip MONTHLY rasters -----------------------------------

message("\n=== Downloading MONTHLY rasters ===")

for (var in VARIABLES) {
  message(sprintf("\nVariable: %s", var))

  get_prism_monthlys(
    type       = var,
    years      = monthly_years,
    resolution = "800m",
    keepZip    = KEEP_ZIP
  )

  folders <- prism_archive_subset(var, "monthly", resolution = "800m") |>
    pd_to_file() |> dirname() |> unique()

  # Only clip files not already clipped (allows safe re-runs)
  folders_to_clip <- folders[!file.exists(
    file.path(PRISM_CLIP_DIR,
              paste0(basename(folders), "_clipped.tif"))
  )]

  for (folder in folders_to_clip) {
    clip_and_save(folder, aoi, PRISM_CLIP_DIR,
                  delete_original = DELETE_AFTER_CLIP)
  }

  Sys.sleep(2)   # polite pause between variables
}


# --- 7. Download and clip DAILY rasters -------------------------------------

message("\n=== Downloading DAILY rasters ===")

for (var in VARIABLES) {
  message(sprintf("\nVariable: %s (daily)", var))

  for (i in seq_len(nrow(daily_ym))) {
    yr <- daily_ym$year[i]
    mo <- daily_ym$month[i]

    min_date <- sprintf("%04d-%02d-01", yr, mo)
    max_date <- format(
      seq(as.Date(min_date), by = "month", length.out = 2)[2] - 1,
      "%Y-%m-%d"
    )

    tryCatch({
      get_prism_dailys(
        type       = var,
        minDate    = min_date,
        maxDate    = max_date,
        resolution = "800m",
        keepZip    = KEEP_ZIP
      )

      folders <- prism_archive_subset(var, "daily", resolution = "800m") |>
        pd_to_file() |> dirname() |> unique()

      ym_str        <- sprintf("%04d%02d", yr, mo)
      folders_month <- folders[grepl(ym_str, basename(folders))]
      folders_new   <- folders_month[!file.exists(
        file.path(PRISM_CLIP_DIR,
                  paste0(basename(folders_month), "_clipped.tif"))
      )]

      for (folder in folders_new) {
        clip_and_save(folder, aoi, PRISM_CLIP_DIR,
                      delete_original = DELETE_AFTER_CLIP)
      }

      Sys.sleep(0.5)
    }, error = function(e) {
      warning(sprintf("Failed: daily %s %04d-%02d: %s",
                      var, yr, mo, e$message))
    })
  }

  Sys.sleep(2)
}

# Summary
clipped_files <- list.files(PRISM_CLIP_DIR, pattern = "_clipped\\.tif$")
message(sprintf(
  "\n=== Download complete ===\n%d clipped GeoTiffs saved in %s",
  length(clipped_files), PRISM_CLIP_DIR
))


# =============================================================================
# NEXT STAGE STUB — run this separately once coordinates are revised
# Save as: extract_climate_values.R
# =============================================================================
#
# library(terra); library(tidyverse)
#
# PRISM_CLIP_DIR <- "data/prism_clipped"
#
# # Load revised sample coordinates
# samples_revised <- read_csv("data/samples_revised.csv")
# # Expected columns: sample_id, lat_final, lon_final, year, month
#
# variables <- c("ppt","tmax","tmin","tmean","tdmean","vpdmin","vpdmax")
#
# extract_at_point <- function(tif_path, lon, lat) {
#   r  <- rast(tif_path)
#   pt <- vect(matrix(c(lon, lat), ncol=2), type="points", crs="EPSG:4326")
#   pt <- project(pt, crs(r))
#   as.numeric(extract(r, pt)[1, 2])
# }
#
# results <- samples_revised |>
#   rowwise() |>
#   mutate(across(
#     all_of(variables),
#     ~ {
#       pattern <- sprintf(".*_%s_.*800m.*%04d%02d.*_clipped\\.tif$",
#                          cur_column(), year, month)
#       tif <- list.files(PRISM_CLIP_DIR, pattern=pattern, full.names=TRUE)[1]
#       if (is.na(tif)) NA_real_ else extract_at_point(tif, lon_final, lat_final)
#     }
#   )) |>
#   ungroup()
#
# write_csv(results, "data/prism_extracted.csv")
