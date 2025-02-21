#' Update Gene Symbols and Create a Seurat Object
#'
#' @description
#' \code{updateGenes} updates gene symbols (e.g., row names in a count matrix) to
#' their most current HGNC-approved symbols, aggregates counts for duplicated
#' updated symbols, and creates a Seurat object.
#'
#' @param tmp A matrix (or sparse matrix) of counts where each row is a gene and
#'   each column is a cell/sample.
#' @param sample A character string specifying the project name or sample name
#'   used for creating the Seurat object.
#' @param gene.symbols A data frame containing HGNC symbol information. By default,
#'   uses \code{hgnc.table} from the \pkg{HGNChelper} package.
#'
#' @details
#' This function attempts to map old or alias gene symbols to their updated,
#' HGNC-approved symbols. It then sums counts across genes that collapse to the
#' same updated symbol, and constructs a Seurat object from the resulting
#' aggregated matrix.
#'
#' Internally, the function:
#' \enumerate{
#'   \item Checks each gene symbol in \code{rownames(tmp)} for approval status.
#'   \item Updates non-approved symbols to their approved equivalent if a match
#'         is found in \code{gene.symbols}.
#'   \item Aggregates counts across duplicated updated symbols.
#'   \item Creates a \code{\link[Seurat]{Seurat}} object using the updated count matrix.
#' }
#'
#' @return A Seurat object with updated gene symbols (in the row names).
#'
#' @importFrom HGNChelper hgnc.table
#' @importFrom SeuratObject CreateAssayObject CreateSeuratObject
#' 
#' @seealso \code{\link[HGNChelper]{hgnc.table}} for the data on HGNC symbols,
#'   \code{\link[SeuratObject]{CreateSeuratObject}} for creating a Seurat object.
#'
#' @examples
#' \dontrun{
#' # Toy count matrix
#' mat <- matrix(
#'   c(10, 5, 0, 3, 2, 8),
#'   nrow = 3,
#'   ncol = 2,
#'   dimnames = list(c("GENEA", "GENEB", "OLDNAME"), c("Cell1", "Cell2"))
#' )
#' 
#' # Run the function
#' seurat_obj <- updateGenes(tmp = mat, sample = "ToySample")
#'
#' # You can inspect the updated Seurat object:
#' seurat_obj
#' }
#'
#' @export
updateGenes <- function(tmp,
                        sample,
                        gene.symbols = HGNChelper::hgnc.table) {
  

  # Helper function: Check and update gene names using a HGNC reference table
  CheckAndUpdateGenes <- function(genes, gene.symbols) {
    approved <- genes %in% gene.symbols$Approved.Symbol
    idx <- match(genes, gene.symbols$Symbol)
    to_update <- !approved & !is.na(idx)
    updated_genes <- genes
    updated_genes[to_update] <- gene.symbols$Approved.Symbol[idx[to_update]]
    not_found_genes <- genes[!approved & is.na(idx)]
    list(updated_genes = updated_genes, not_found_genes = not_found_genes)
  }
  

  # 1. Update the gene symbols
  updated.features <- CheckAndUpdateGenes(rownames(tmp), gene.symbols)
  new.features <- updated.features$updated_genes
  

  # 2. Sum counts for genes with the same updated symbol
  tmp_matrix <- as.matrix(tmp)
  summed_matrix <- rowsum(tmp_matrix, group = new.features, reorder = FALSE)
  
  # Create a Seurat Assay
  tmp_assay <- Seurat::CreateAssayObject(counts = summed_matrix)
  
  SeuratObj <- Seurat::CreateSeuratObject(
    counts  = tmp_assay,
    assay   = "RNA",
    project = sample
  )
  
  return(SeuratObj)
}

