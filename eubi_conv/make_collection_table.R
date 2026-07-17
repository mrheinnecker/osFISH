suppressPackageStartupMessages({
  library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
arg_value <- function(name, default="") {
  hit <- which(args %in% paste0("--", name))
  if (length(hit) && hit[1] < length(args)) return(args[[hit[1] + 1]])
  hit <- grep(paste0("^--", name, "="), args)
  if (length(hit)) return(sub(paste0("^--", name, "="), "", args[[hit[1]]]))
  default
}

all_s3 <- arg_value("all_s3")
all_images <- arg_value("all_images")
s3_bucket <- arg_value("s3_bucket", "s3embl/temscreen/osFISH")
local_collection_table <- arg_value("local_collection_table", "osfish_collection_table.tsv")
sheet_mode <- arg_value("sheet_mode", "local")
google_key <- arg_value("google_key", "")
collection_table_url <- arg_value("collection_table_url", "")
collection_table_sheet <- arg_value("collection_table_sheet", "ct")
metadata_root <- arg_value("metadata_root", "")

if (all_s3 == "" || !file.exists(all_s3)) stop("--all_s3 is required")
if (all_images == "" || !file.exists(all_images)) stop("--all_images is required")

parse_mc_ls_path <- function(line) {
  parsed <- str_match(line, "^\\[.*?\\]\\s+\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2]
  ifelse(
    is.na(parsed),
    str_match(line, "^\\S+\\s+(?:STANDARD\\s+)?(.+)$")[, 2],
    parsed
  )
}

s3_public_prefix <- function(bucket) {
  bucket_path <- bucket %>%
    str_remove("^s3embl/") %>%
    str_remove("/$")
  file.path("https://s3.embl.de", bucket_path)
}

near_square_columns <- function(n) {
  if (n <= 1) return(1L)
  as.integer(ceiling(sqrt(n)))
}

channel_colors <- c("cyan", "magenta", "yellow", "white", "green", "blue", "red")

read_image_metadata <- function(path) {
  read_tsv(path, col_types=cols(.default=col_character())) %>%
    mutate(
      metadata_file=as.character(path),
      metadata_mtime=as.numeric(file.info(path)$mtime)
    )
}

metadata_files <- all_images
if (metadata_root != "" && dir.exists(metadata_root)) {
  discovered <- list.files(
    metadata_root,
    pattern="^all_images[.]tsv$",
    recursive=TRUE,
    full.names=TRUE
  )
  metadata_files <- c(metadata_files, discovered)
}
metadata_files <- metadata_files[file.exists(metadata_files)]
metadata_files <- unique(normalizePath(metadata_files, mustWork=TRUE))

images <- map_dfr(metadata_files, read_image_metadata)
if (!"run_name" %in% names(images)) images$run_name <- "unknown_run"
if (!"run_date" %in% names(images)) images$run_date <- ""

images <- images %>%
  mutate(
    channels=as.integer(channels),
    channels=if_else(is.na(channels) | channels < 1L, 1L, channels),
    source_name=output_name,
    run_name=coalesce(run_name, "unknown_run"),
    run_date=coalesce(run_date, ""),
    site=condition,
    view=paste(run_name, condition, sep="_"),
    grid=paste(run_name, condition, sep="_"),
    format="OmeZarr",
    blend="sum",
    exclusive=TRUE
  ) %>%
  arrange(desc(metadata_mtime)) %>%
  distinct(source_name, .keep_all=TRUE)

s3_entries <- read_lines(all_s3) %>%
  as_tibble() %>%
  mutate(
    s3_raw=parse_mc_ls_path(value),
    s3_raw=str_remove(s3_raw, "/$"),
    s3_prefix=str_match(s3_raw, "^([^/]+\\.ome\\.zarr)$")[, 2],
    source_name=str_remove(s3_prefix, "\\.ome\\.zarr$")
  ) %>%
  filter(!is.na(source_name), source_name != "") %>%
  distinct(source_name, .keep_all=TRUE) %>%
  transmute(
    source_name,
    uri=paste0(s3_public_prefix(s3_bucket), "/", s3_prefix, "/")
  )

base_table <- images %>%
  inner_join(s3_entries, by="source_name") %>%
  group_by(grid) %>%
  arrange(source_name, .by_group=TRUE) %>%
  mutate(
    grid_index=row_number() - 1L,
    grid_columns=near_square_columns(n()),
    grid_position=paste0("(", grid_index %% grid_columns, ",", grid_index %/% grid_columns, ")")
  ) %>%
  ungroup()

if (nrow(base_table) == 0) {
  stop("No uploaded osFISH OME-Zarr datasets matched the conversion metadata.")
}

collection_table <- base_table %>%
  mutate(channel=map(channels, ~seq_len(.x) - 1L)) %>%
  unnest(channel) %>%
  mutate(
    display=paste0("channel_", channel),
    color=channel_colors[(channel %% length(channel_colors)) + 1L],
    name=paste0(source_name, "_c", channel),
    contrast_limits="",
    x_scale_nm=x_scale_nm,
    y_scale_nm=y_scale_nm,
    z_scale_nm=z_scale_nm
  ) %>%
  select(
    uri,
    name,
    view,
    grid,
    grid_position,
    channel,
    display,
    color,
    contrast_limits,
    blend,
    format,
    exclusive,
    source_name,
    run_name,
    run_date,
    condition,
    image_name,
    x_scale_nm,
    y_scale_nm,
    z_scale_nm,
    width_px,
    height_px
  )

write_tsv(collection_table, local_collection_table)

if (tolower(sheet_mode) == "google") {
  if (collection_table_url == "") {
    stop("--collection_table_url is required when --sheet_mode google")
  }
  if (google_key == "" || !file.exists(google_key)) {
    stop("--google_key must point to an existing service-account JSON when --sheet_mode google")
  }

  suppressPackageStartupMessages({
    library(googlesheets4)
    library(googledrive)
  })
  gs4_auth(path=google_key)
  drive_auth(path=google_key)
  write_sheet(collection_table, ss=collection_table_url, sheet=collection_table_sheet)
}

write_tsv(tibble(done="done"), "done.tsv")
