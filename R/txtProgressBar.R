txtProgressBar=function (min = 0, max = 1, initial = 0, char = "=", width = NA,
          title, label, style = 1, file = "") {
  if (!identical(file, "") && !(inherits(file, "connection") &&
                                isOpen(file)))
    stop("'file' must be \"\" or an open connection object")
  if (!style %in% 1L:3L)
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L)
      width <- width - 10L
    width <- trunc(width/nw)
  }
  if (max <= min)
    stop("must have 'max' > 'min'")
  up <- function(value) {
    if (!is.finite(value) || value < min || value > max)
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc)
      return()
    cat(paste0(" Program is running..be patient...\n Bootstrap is running..be patient...\n\r  |",
               strrep(" ", nw * width +
                                 6)), file = file)
    cat(paste(c(" Program is running..be patient...\n Bootstrap is running..be patient...\n\r  |",
                rep.int(char, nb), rep.int(" ",
                                                    nw * (width - nb)), sprintf("| %3d%%", pc)),
              collapse = ""), file = file)
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n", file = file)
    flush.console()
    .killed <<- TRUE
  }
  up(initial)
  structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}
