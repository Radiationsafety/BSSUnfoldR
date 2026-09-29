files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
hits <- list()
for (f in files) {
  L <- readLines(f, warn = FALSE)
  idx <- which(grepl("^#'", L))
  if (!length(idx)) next
  # split into contiguous blocks
  brk <- which(diff(idx) > 1)
  starts <- c(idx[1], if (length(brk)) idx[brk + 1])
  ends <- c(if (length(brk)) idx[brk], tail(idx, 1))
  for (k in seq_along(starts)) {
    blk <- L[starts[k]:ends[k]]
    if (any(grepl("^#' *@export\\s*$", blk)) && any(grepl("^#' *@keywords internal", blk))) {
      fn <- if (ends[k] < length(L)) L[ends[k] + 1] else "?"
      hits[[length(hits) + 1]] <- data.frame(
        file = f, line = starts[k],
        fn = sub(" *<-.*", "", fn), stringsAsFactors = FALSE
      )
    }
  }
}
d <- do.call(rbind, hits)
cat("blocks with @export + @keywords internal:", nrow(d), "\n")
print(d, row.names = FALSE)
