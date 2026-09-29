rd <- list.files("man", pattern = "\\.Rd$")
topic <- sub("\\.Rd$", "", rd)
is_internal <- function(f) {
  x <- readLines(f, warn = FALSE)
  any(x == "\\keyword{internal}") || any(grepl("keyword\\{internal\\}", x, fixed = TRUE))
}
internal <- vapply(rd, function(f) is_internal(file.path("man", f)), logical(1))
pub <- sort(topic[!internal])

yml <- readLines("_pkgdown.yml", warn = FALSE)
in_contents <- FALSE
listed <- character()
for (l in yml) {
  if (grepl("^ *contents: *$", l)) { in_contents <- TRUE; next }
  if (!in_contents) next
  if (grepl("^ *- ", l)) {
    listed <- c(listed, sub("^ *- ?", "", l))
  } else if (grepl("^ *[^ ]", l) && !grepl("^ *- ", l)) {
    in_contents <- grepl("^ *contents: *$", l)
  }
}
listed <- unique(gsub("[`'\"]", "", trimws(listed)))
cat("public Rd topics:", length(pub), " | yml index entries:", length(listed), "\n")
miss <- setdiff(pub, listed)
cat("MISSING FROM INDEX (", length(miss), "):\n", sep = "")
cat(miss, sep = "\n")
cat("\nINDEX ENTRIES WITH NO MATCHING Rd (", length(setdiff(listed, topic)), "):\n", sep = "")
cat(setdiff(listed, topic), sep = "\n")
