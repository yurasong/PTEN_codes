# PTEN_codes

**Publication**

This repository supports the analyses and visualizations for the manuscript:

> Jiang C., Song Y., *et al.* (2025). **Innate immunity and the NF-κB pathway control prostate stem cell
plasticity, reprogramming and tumor initiation.** _Nature Cancer_. **Accepted**, currently in final proof stage.

A DOI and link to the final publication will be provided here once available.

---

## Contents

- **02_Marker_Genes** – Scripts for identifying marker genes in scRNA-seq data.
- **04_Slingshot_scripts** – Pseudotime trajectory analyses using Slingshot.
- **06_SCENIC** – Gene regulatory network inference with pySCENIC.
- **07_Human_public_Data** – Processing pipelines for public human single-cell datasets.
- **08_bulKRNA-seq** – Bulk RNA-seq alignment and analysis pipelines.
- **09_bulkATAC-seq** – Bulk ATAC-seq processing: alignment, peak calling, and motif prediction.
- **R scripts in root: Related to scRNA-seq data analysis**:
  - `01_UMAP.R` – UMAP embedding and visualization.
  - `02_ScatterPlot.R` – Generation of scatter plots.
  - `03_marker_genes.R` – Visualization of marker gene expression.
  - `04-01_Slingshot_analysis_HY_BC_Prox.R` – Slingshot analysis for HY_BC_Prox lineage.
  - `04-02_Slingshot_analysis_HY_Nkx3-1.R` – Slingshot analysis for HY_Nkx3-1 lineage.
  - `05_DEG_analysis.R` – Differential expression analysis.

---

## Getting Started

Clone the repository:
```bash
git clone https://github.com/yurasong/PTEN_codes.git
cd PTEN_codes
```

Each code file contains its own instructions for running the specific analyses.

---

## Data Availability

All raw sequencing datasets that support the findings of this study have been deposited in the NCBI Gene Expression Omnibus (GEO) under the following accession numbers:

- **[GSE270187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270187)**: scRNA-seq (WT, PTEN 6W, PTEN 10M)
- **[GSE270189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270189)**: bulk RNA-seq, Rapamycin treated  
- **[GSE270190](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270190)**: bulk RNA-seq, cell fate upon PTEN deletion in BC  
- **[GSE286018](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286018)**: bulk RNA-seq of P65/Pten knock BCs  
- **[GSE286019](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286019)**: bulk RNA-seq of basal-derived and luminal-derived tumors  
- **[GSE270191](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270191)**: ATAC-seq of pooled BC/LC of WT and PTEN-deleted cells  
- **[GSE288787](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288787)**: ATAC-seq of BCs from different lobes upon PTEN deletion

Previously published human prostate cancer datasets re‐analyzed in this study were available via:

- Song, H., *et al*. Nat Commun (2022): GitHub repository of the author ([scRNA-seq-Analysis-of-Prostate-Cancer-Samples](https://github.com/franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples))  
- Tuong ZK., *et al*. Cell Rep (2021): the Prostate Cell Atlas ([prostatecellatlas.org](https://www.prostatecellatlas.org))  
- Hirz T., *et al*. Nat Commun (2023): [GSE181294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181294)  
- Wong, H.Y., *et al*. Nat Commun 13 (2022): [GSE185344](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE185344)

---

## License

All rights are reserved by the authors.
