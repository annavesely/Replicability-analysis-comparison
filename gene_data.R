
require(curatedOvarianData)
require(KEGGREST)


# -------------------------------------------------------


# P-values from single study
# eset: gene expression set with m genes
# output: vector of length m containing p-values for each gene

compute_pvals <- function(eset){
  
  group <- eset$sample_type
  group <- ifelse(group %in% c("tumor", "metastatic"), "tumor", "healthy")
  
  expr <- exprs(eset)
  
  out <- apply(expr, 1, function(x) t.test(x ~ group)$p.value)
  return(out)
}


# -------------------------------------------------------


# KEGG pathways
# pathway_id: identification number of the pathway
# output: vector of genes inside the pathway

get_kegg_genes <- function(pathway_id) {
  
  pathway <- KEGGREST::keggGet(pathway_id)[[1]]
  out <- sapply(pathway$GENE, function(x) x[1])
  out <- sub(";.*", "", out)
  names(out) <- NULL
  
  return(out)
}


# -------------------------------------------------------

data("GSE18520_eset")
data("GSE6008_eset")
data("GSE26712_eset")
data("GSE9891_eset")
data("GSE51088_eset")

pvals1 <- compute_pvals(GSE18520_eset)
pvals2 <- compute_pvals(GSE6008_eset)
pvals3 <- compute_pvals(GSE26712_eset)
pvals4 <- compute_pvals(GSE9891_eset)
pvals5 <- compute_pvals(GSE51088_eset)


common_genes <- Reduce(intersect, list(
  names(pvals1), names(pvals2), names(pvals3), names(pvals4), names(pvals5)
))


pmat0 <- cbind(
  pvals1[common_genes],
  pvals2[common_genes],
  pvals3[common_genes],
  pvals4[common_genes],
  pvals5[common_genes]
)

rownames(pmat0) <- common_genes


# -------------------------------------------------------


# PI3K-Akt signaling
pi3k_genes <- get_kegg_genes("hsa04151")

# p53 signaling
p53_genes <- get_kegg_genes("hsa04115")


# -------------------------------------------------------


save(pmat0, pi3k_genes, p53_genes, file="ovarian_cancer_data.RData")
