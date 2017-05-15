#' @export
bed_chromsweep <- function(x, y, genome, ...){
  # arrange x and y chroms according to genome file
  # to ensure the same sort order
  x <- x[order(match(x$chrom, genome$chrom)), ]
  y <- y[order(match(y$chrom, genome$chrom)), ]
  res <- chromsweep_impl(x, y, genome)
  #filter for basic by strandedness, etc.
}
