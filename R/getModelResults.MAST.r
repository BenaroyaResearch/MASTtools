#' Calculate model statistics for MAST model fitting results
#'
#' Calculate log-fold-changes and significance values for a fitted MAST model, and format the results
#' for more straightforward viewing.
#' @param zlm_result The zlm fit results, generally the output of \code{MAST::zlm.SingleCellAssay}
#' @param lrTest_result The likelihood-ratio test results, generally the output of \code{MAST::lrTest}. Needed for determinations of significance.
#' @param contrast character, the name of the contrast for which to determine significance (and possibly logFC).
#' @param contrast.type character, the type of contrast being extracted. Must be specified, as logFC can only be extracted for binary and/or continuous contrasts. Fro categorical contrasts, only significance values are extracted.
#' @param method.p_adjust character, the method to be used for calculating adjusted p-values. Passed to \code{p.adjust}. Defaults to "BH".
#' @export
#' @return A data frame with one row per gene, and columns for the gene name, (optionally logFC,) p-value, and adjusted p-values.
#' @usage \code{
#' getModelResults.MAST(
#'   zlm_result, lrTest, contrast,
#'   contrast.type,
#'   method.p_adjust="BH")}
getModelResults.MAST <-
  function(zlm_result, lrTest, contrast,
           contrast.type,
           method.p_adjust="BH") {
    contrast.type <-
      match.arg(contrast.type,
                choices=c("binary", "continuous", "categorical"))
    if (contrast.type %in% c("binary", "continuous")) {
      results <- getLogFC(zlm_result)
      results <- results[results[,"contrast"]==contrast,]
      colnames(results)[colnames(results)=="primerid"] <- "gene"
    } else if (contrast.type %in% "categorical") {
      results <- data.frame(gene=rownames(zlm_result@coefC))
    }
      
    results[,"P.Val"] <- lrTest[,"hurdle","Pr(>Chisq)"]
    results[,"adj.P.Val"] <-
      p.adjust(results[,"P.Val"], method=method.p_adjust)
    
    rownames(results) <- results[,"gene"]
    
    results
}
