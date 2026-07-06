library(tidyverse)
library(getopt)
library(googlesheets4)
library(googledrive)

spec <- matrix(c(
  "collection_table", "t", 1, "character",
  "google_key", "k", 1, "character",
  "collection_table_url", "u", 1, "character",
  "collection_table_sheet", "s", 1, "character"
), ncol=4, byrow=TRUE)
opt <- getopt(spec)

value_or_default <- function(value, default) {
  if (is.null(value) || is.na(value) || value == "") default else value
}

collection_table <- value_or_default(opt$collection_table, "")
google_key <- value_or_default(opt$google_key, "")
collection_table_url <- value_or_default(opt$collection_table_url, "")
collection_table_sheet <- value_or_default(opt$collection_table_sheet, "osfish_collection_table")

if (collection_table == "" || !file.exists(collection_table)) {
  stop("--collection_table must point to an existing TSV file")
}
if (google_key == "" || !file.exists(google_key)) {
  stop("--google_key must point to a Google service-account JSON key")
}
if (collection_table_url == "") {
  stop("--collection_table_url is required")
}

col_table <- read_tsv(collection_table, col_types=cols(.default=col_character()))

gs4_auth(path=google_key)
drive_auth(path=google_key)
write_sheet(col_table, ss=collection_table_url, sheet=collection_table_sheet)

write_tsv(tibble(done="done"), file="collection_table_upload_done.txt")
