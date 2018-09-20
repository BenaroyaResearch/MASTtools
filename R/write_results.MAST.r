#' Write MAST model fitting results to file(s)
#'
#' Output file(s) with results from MAST model fitting
#' @param lfc_result The fit results, generally the output of \code{MASTtools::getModelResults.MAST}
#' @param write_full_results logical, whether to output the full results data frame to a file. Defaults to TRUE.
#' @param write_sigGenes logical, whether to output a list of the genes with adjusted p-values below a threshold. Defaults to TRUE.
#' @param thresholds.adj.P.Val numeric, a vector of adjusted p-value thresholds to use for gene significance cutoffs. One file is generated for each threshold value. Defaults to 0.05. Ignored if \code{write_sigGenes == FALSE}
#' @param write_topGenes logical, whether to output a sorted list of the genes with strongest significance. Defaults to TRUE.
#' @param n_topGenes numeric, a vector of numbers of most significant genes to output. One file is generated for each threshold value. Defaults to 100. Ignored if \code{write_topGenes == FALSE}
#' @param file_prefix character, the path to the file location. "MAST_..." will be appended to it.
#' @param test_id character, the name of the model test, to be included in the filename.
#' @export
#' @usage \code{
#' write_results.MAST(lfc_result,
#'                    write_full_results=TRUE,
#'                    write_sigGenes=TRUE, thresholds.adj.P.Val=0.05,
#'                    write_topGenes=TRUE, n_topGenes=100,
#'                    file_prefix, test_id)}
write_results.MAST <-
  function(lfc_result,
           write_full_results=TRUE,
           write_sigGenes=TRUE, thresholds.adj.P.Val=0.05,
           write_topGenes=TRUE, n_topGenes=100,
           file_prefix, test_id) {
    lfc_result <- lfc_result[order(lfc_result$P.Val),]
    if (write_full_results) {
      write.table(
        lfc_result,
        paste0(file_prefix, "MAST_Results.", test_id, ".txt"),
        quote = FALSE, col.names=TRUE, row.names=FALSE)
      
    }
    if (write_sigGenes) {
      for (threshold.tmp in thresholds.adj.P.Val) {
        genes.tmp <- lfc_result$gene[lfc_result$adj.P.Val < threshold.tmp]
        write.table(
          genes.tmp,
          paste0(file_prefix, "MAST_", threshold.tmp, "_SigGenes.", test_id, ".txt"),
          quote = FALSE, col.names=FALSE, row.names=FALSE)
      }
    }
    if (write_topGenes) {
      for (n.tmp in n_topGenes) {
        genes.tmp <- lfc_result$gene[1:n.tmp]
        write.table(
          genes.tmp,
          paste0(file_prefix, "MAST_top", n.tmp, "_genes.", test_id, ".txt"),
          quote = FALSE, col.names=FALSE, row.names=FALSE)
      }
    }
}
