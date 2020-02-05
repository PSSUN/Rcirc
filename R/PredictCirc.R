##' @export
##' @author Peisne Sun

PredictCirc <- function(sam, fa, out, gtf) {
  sam <- sam
  fa <- fa
  out <- out
  gtf <- gtf
  ciri <- system.file("extdata", "CIRI2.pl", package = "Rcirc")
  ciri_command <- paste("perl", ciri)
  input_sam <- paste("-I", sam)
  input_out <- paste("-O", out)
  input_fa <- paste("-F", fa)
  input_gtf <- paste("-A", gtf)
  full_command <- paste(ciri_command, input_sam, input_out, input_fa, input_gtf)
  system(full_command)
}
