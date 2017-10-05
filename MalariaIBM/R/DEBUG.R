
#' debug
#'
#' get the bugs out!
#'
#' @export
DEBUG_CLASS <- function(){

  MalariaIBM:::.R6_test$set(which = "public",name = "finalize",
    value = function(){cat("test being garbage collected\n",sep="")},overwrite = overwrite
  )

}
