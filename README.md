# Replicability analysis comparison

## Description

This repository provides R code to compare three methods for replicability analysis:

1. **BHY** – Benjamini, Heller, and Yekutieli (2009): a partial conjunction FDR-based method that sequentially tests replicability levels across studies.  
2. **adaFilter** – Wang et al. (2022): an adaptive filtering approach controlling FDR for replicability at a pre-specified level.  
3. **ARI** – Goeman and Solari (2011): a closed testing procedure providing simultaneous lower confidence bounds on the number of studies in which a feature is active (true discovery guarantee).  

The code is implemented in **`code.R`** and can be applied to any matrix of p-values, where rows correspond to features and columns correspond to studies.  

An example dataset, as illustrated in the paper, is pre-processed in **`gene_data.R`** and saved as **`ovarian_cancer_data.RData`**.  

For detailed discussion of method properties, assumptions, and inferential guarantees, see Vesely (2026).

---

## Example analysis

After loading the example dataset and setting the significance level, apply the three methods.

```r
load("ovarian_cancer_data.RData")
alpha <- 0.05

res_bhy <- BHY(pmat0, alpha=alpha)
round(prop.table(table(res_bhy) * 100), 2)

res_ada <- adafilter(pmat0, alpha=alpha)
round(prop.table(table(res_ada) * 100), 2)

res_ari <- ari(pmat0, alpha=alpha)
round(prop.table(table(res_ari) * 100), 2)
```

ARI's results can also be reported for subsets. For instance, we consider the pathways PI3K-Akt and p53.

```r
genes <- rownames(pmat0)

# PI3K-Akt
sel1 <- which(genes %in% pi3k_genes)
round(prop.table(table(res_ari[sel1])) * 100, 2)

# p53
sel2 <- which(genes %in% p53_genes)
round(prop.table(table(res_ari[sel2])) * 100, 2)
```
