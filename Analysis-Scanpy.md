# 🧬 Analysis and Visualization of Spatial Transcriptomics Data using Scanpy


---

# Introduction

Spatial transcriptomics is an advanced sequencing technology that allows researchers to study gene expression while preserving the exact spatial location of cells within tissues. Unlike traditional RNA sequencing methods, spatial transcriptomics helps identify not only which genes are expressed but also where they are expressed inside the tissue architecture. This provides important biological insights into tissue organization, cell communication, tumor microenvironments, developmental biology, and disease progression.

This tutorial demonstrates the analysis and visualization of spatial transcriptomics data using the Python library Scanpy. The workflow mainly focuses on 10x Genomics Visium data and additionally provides an example of MERFISH spatial data analysis. The complete workflow includes data loading, quality control, preprocessing, normalization, dimensionality reduction, clustering, visualization in spatial coordinates, and marker gene analysis.

---

#  Importing Required Libraries

```python
from __future__ import annotations
```

This enables postponed evaluation of type annotations in Python and improves compatibility and performance when using type hints.

---

```python
import matplotlib.pyplot as plt
```

Matplotlib is imported for generating visualizations such as histograms, UMAP plots, spatial plots, and heatmaps.

---

```python
import pandas as pd
```

Pandas is used for handling tabular datasets such as coordinate files, metadata tables, and count matrices.

---

```python
import scanpy as sc
```

Scanpy is the primary library used for single-cell and spatial transcriptomics analysis. It provides tools for preprocessing, clustering, dimensionality reduction, and visualization.

---

```python
import seaborn as sns
```

Seaborn is imported for creating statistical visualizations such as histograms and distribution plots.

---

#  Scanpy Configuration

```python
sc.logging.print_versions()

sc.set_figure_params(
    facecolor="white",
    figsize=(8, 8)
)

sc.settings.verbosity = 3
```

The first command prints the versions of Scanpy and all installed dependencies to ensure reproducibility of the analysis. The second command configures default plotting parameters by setting a white figure background and a default figure size of 8 × 8 inches. The final command increases the verbosity level so that Scanpy displays detailed information during each processing step, helping the user track the workflow more clearly.

---

#  Reading the Visium Spatial Transcriptomics Dataset

```python
adata = sc.datasets.visium_sge(
    sample_id="V1_Human_Lymph_Node"
)

adata.var_names_make_unique()

adata.var["mt"] = adata.var_names.str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    inplace=True
)
```

The dataset used in this tutorial is a publicly available 10x Genomics Visium human lymph node dataset. The `visium_sge()` function automatically downloads and loads the dataset into an AnnData object named `adata`. This object stores gene expression counts, metadata, spatial coordinates, and tissue images in a structured format.

The command `var_names_make_unique()` ensures that all gene names are unique because duplicate gene names can create problems during downstream analysis. Next, mitochondrial genes are identified by checking whether gene names begin with `"MT-"`. These genes are stored in a new variable called `"mt"` because mitochondrial gene expression is commonly used during quality control to identify damaged or stressed cells.

The `calculate_qc_metrics()` function computes several quality control metrics including total counts per spot, number of detected genes, mitochondrial read percentages, and other sequencing statistics. These values are stored within the metadata of the AnnData object for later filtering and visualization.

---

#  Structure of the AnnData Object

```python
adata
```

```python
AnnData object with n_obs × n_vars = 4035 × 36601
```

The AnnData object contains 4035 spatial spots and 36601 genes. Spatial transcriptomics data in Scanpy is organized into different components. The `obs` section stores information related to spots or cells such as counts and QC metrics. The `var` section stores gene-related information. The `uns` section contains unstructured data including tissue images, while the `obsm` section stores multidimensional information such as spatial coordinates and dimensionality reduction embeddings.

---

#  Quality Control and Visualization

```python
fig, axs = plt.subplots(1, 4, figsize=(15, 4))

sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0]
)

sns.histplot(
    adata.obs["total_counts"][
        adata.obs["total_counts"] < 10000
    ],
    kde=False,
    bins=40,
    ax=axs[1],
)

sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    bins=60,
    ax=axs[2],
)

sns.histplot(
    adata.obs["n_genes_by_counts"][
        adata.obs["n_genes_by_counts"] < 4000
    ],
    kde=False,
    bins=60,
    ax=axs[3],
)
```

This section visualizes quality control metrics using histograms. The first histogram displays the distribution of total sequencing counts across all spatial spots. The second histogram zooms into spots having fewer than 10,000 counts to provide a clearer view of lower-count distributions. The third histogram visualizes the number of detected genes per spot, while the fourth histogram focuses specifically on spots with fewer than 4000 detected genes. These visualizations help identify low-quality spots, sequencing artifacts, and outliers before proceeding with downstream analysis.
<img width="2025" height="608" alt="imagy one" src="https://github.com/user-attachments/assets/3b374f5c-33f5-450c-a19e-b367a93e2d83" />


---

#  Filtering Low-Quality Spots and Genes

```python
sc.pp.filter_cells(
    adata,
    min_counts=5000
)

sc.pp.filter_cells(
    adata,
    max_counts=35000
)

adata = adata[
    adata.obs["pct_counts_mt"] < 20
].copy()

print(f"#cells after MT filter: {adata.n_obs}")

sc.pp.filter_genes(
    adata,
    min_cells=10
)
```

Quality filtering removes low-quality spots and genes that may introduce noise into the dataset. The first filtering step removes spots with fewer than 5000 total counts because such spots often contain insufficient biological information. The second filtering step removes spots with more than 35,000 counts because unusually high counts may indicate technical artifacts or doublets.

The dataset is further filtered by removing spots with more than 20% mitochondrial gene expression. High mitochondrial content is generally associated with damaged or dying cells. After spot filtering, genes expressed in fewer than 10 spots are removed because extremely rare genes contribute little useful information and increase computational noise.

---

# Data Normalization and Highly Variable Genes

```python
sc.pp.normalize_total(
    adata,
    inplace=True
)

sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    flavor="seurat",
    n_top_genes=2000
)
```

Normalization adjusts sequencing counts so that all spots become comparable despite differences in sequencing depth. The `normalize_total()` function scales counts across spots to a common total expression level. Following normalization, logarithmic transformation is applied using `log1p()`, which calculates log(x + 1). This transformation reduces the influence of highly expressed genes and stabilizes variance across the dataset.

The `highly_variable_genes()` function identifies the top 2000 genes showing the greatest variability across spots. Highly variable genes contain the most biologically informative patterns and are commonly used for dimensionality reduction and clustering analyses.

---

#  PCA, Neighborhood Graph Construction, UMAP, and Clustering

```python
sc.pp.pca(adata)

sc.pp.neighbors(adata)

sc.tl.umap(adata)

sc.tl.leiden(
    adata,
    key_added="clusters",
    flavor="igraph",
    directed=False,
    n_iterations=2
)
```

Principal Component Analysis (PCA) reduces the dimensionality of the dataset while preserving the major sources of biological variation. This simplifies downstream computations and visualization. After PCA, the `neighbors()` function constructs a nearest-neighbor graph based on similarities between spots in PCA space.

The UMAP algorithm then projects the high-dimensional data into two-dimensional space for visualization while preserving local neighborhood relationships. Finally, Leiden clustering identifies groups of spots with similar gene expression profiles. These cluster labels are stored within the AnnData object and represent biologically similar tissue regions.

---

# UMAP Visualization

```python
plt.rcParams["figure.figsize"] = (4, 4)

sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "clusters"
    ],
    wspace=0.4
)
```

<img width="2239" height="599" alt="imagy2" src="https://github.com/user-attachments/assets/fe71e975-c6e1-4d5f-8f5f-01f43ddc5664" />
This visualization displays the UMAP embedding colored according to total sequencing counts, number of detected genes, and cluster assignments. Examining these patterns helps determine whether clustering results are driven by biological variation or influenced by technical artifacts such as sequencing depth or low-quality spots.

---

# Spatial Visualization

```python
plt.rcParams["figure.figsize"] = (8, 8)

sc.pl.spatial(
    adata,
    img_key="hires",
    color=[
        "total_counts",
        "n_genes_by_counts"
    ]
)
```
<img width="2200" height="1000" alt="imagy3" src="https://github.com/user-attachments/assets/ade08f80-631a-4dfe-b88a-b68b20ce991e" />

This section visualizes sequencing counts and detected gene numbers directly on the tissue image. The `spatial()` function overlays circular spots on top of the high-resolution histology image provided by 10x Genomics. This allows researchers to observe how sequencing quality and gene expression vary across different tissue regions while preserving tissue morphology.

---

# Spatial Cluster Visualization

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    size=1.5
)
```
<img width="800" height="500" alt="imagy4" src="https://github.com/user-attachments/assets/10fd1bd0-18c6-42ef-ba86-801b8468e8be" />


This visualization displays the spatial arrangement of Leiden clusters across the tissue section. Spots belonging to the same transcriptional cluster often appear close together spatially, indicating organized tissue structures and biologically related regions. Spatial clustering provides important insights into tissue architecture and potential cell-cell interactions.

---

#  Zooming into Specific Tissue Regions

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color="clusters",
    groups=["5", "9"],
    crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5,
    size=1.3,
)
```
<img width="659" height="500" alt="imagy5" src="https://github.com/user-attachments/assets/c9c3ceb1-c8ca-48f5-8f8a-10a1e8c2df89" />

This command focuses on a specific region of the tissue image by cropping selected coordinates. Only clusters 5 and 9 are visualized, making it easier to examine their spatial relationships. Transparency is adjusted using the alpha parameter so the underlying histological tissue image remains visible. This type of zoomed visualization helps researchers investigate localized tissue organization in greater detail.

---

#  Marker Gene Identification

```python
sc.tl.rank_genes_groups(
    adata,
    "clusters",
    method="t-test"
)

sc.pl.rank_genes_groups_heatmap(
    adata,
    groups="9",
    n_genes=10,
    groupby="clusters"
)
```
<img width="657" height="935" alt="imagy6" src="https://github.com/user-attachments/assets/fcfffbd2-ed6c-4877-a38d-7a6909afde86" />

Marker gene analysis identifies genes that are significantly enriched in specific clusters compared to others. The `rank_genes_groups()` function performs differential expression analysis using a t-test. The heatmap visualization displays the top 10 marker genes associated with cluster 9 across all clusters. These marker genes help characterize the biological identity and functional properties of tissue regions.

---

#  Spatial Expression of Marker Genes

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color=["clusters", "CR2"]
)
```
<img width="2287" height="1091" alt="imagy7" src="https://github.com/user-attachments/assets/2349344f-1144-4b85-ac32-3e6d3d912004" />

This visualization compares cluster assignments with the spatial expression pattern of the marker gene `CR2`. The spatial localization of CR2 closely matches the regions occupied by specific clusters, demonstrating how marker genes can validate cluster identities and reveal tissue-specific expression patterns.

---

# Visualization of Additional Genes

```python
sc.pl.spatial(
    adata,
    img_key="hires",
    color=["COL1A2", "SYPL1"],
    alpha=0.7
)
```
<img width="2232" height="1091" alt="imagy8" src="https://github.com/user-attachments/assets/4aa264b4-2fb9-4713-bd90-16c9e326a9e7" />

This section visualizes the spatial expression patterns of the genes `COL1A2` and `SYPL1`. Different genes exhibit distinct localization patterns across the tissue section, reflecting differences in cellular composition and biological function across tissue regions.

---

# MERFISH Spatial Transcriptomics Example

MERFISH is a fluorescence imaging-based spatial transcriptomics technology that provides single-cell resolution. Unlike Visium, which captures gene expression at spot level, MERFISH directly measures transcript abundance within individual cells while preserving precise spatial locations.

---

#  Reading MERFISH Data

```python
coordinates = pd.read_excel(
    "./data/pnas.1912459116.sd15.xlsx",
    index_col=0
)

counts = sc.read_csv(
    "./data/pnas.1912459116.sd12.csv"
).transpose()

adata_merfish = counts[
    coordinates.index,
    :
].copy()

adata_merfish.obsm["spatial"] = coordinates.to_numpy()
```

The coordinate file containing spatial positions is loaded using Pandas, while the gene expression count matrix is read using Scanpy. The transpose operation converts rows and columns into the correct orientation for analysis. A new AnnData object called `adata_merfish` is created, and the spatial coordinates are stored inside the `obsm["spatial"]` slot. This prepares the MERFISH dataset for downstream spatial analysis.

---

# MERFISH Preprocessing and Clustering

```python
sc.pp.normalize_per_cell(
    adata_merfish,
    counts_per_cell_after=1e6
)

sc.pp.log1p(adata_merfish)

sc.pp.pca(
    adata_merfish,
    n_comps=15
)

sc.pp.neighbors(adata_merfish)

sc.tl.umap(adata_merfish)

sc.tl.leiden(
    adata_merfish,
    key_added="clusters",
    resolution=0.5,
    n_iterations=2,
    flavor="igraph",
    directed=False,
)
```

The MERFISH dataset undergoes standard preprocessing steps including normalization, logarithmic transformation, PCA, neighbor graph construction, UMAP embedding, and Leiden clustering. These steps reduce dimensionality, identify transcriptionally similar cells, and group cells into clusters based on shared gene expression patterns.

---

# MERFISH Visualization

```python
sc.pl.umap(
    adata_merfish,
    color="clusters"
)

sc.pl.embedding(
    adata_merfish,
    basis="spatial",
    color="clusters"
)
```
<img width="700" height="600" alt="imagy9" src="https://github.com/user-attachments/assets/6b7fba79-031a-4a58-b51e-8fabf5c8d03d" />

The first plot visualizes clusters in UMAP space, where transcriptionally similar cells appear close together. The second plot visualizes the same clusters in actual spatial coordinates. Since the experiment used cultured U2-OS cells rather than organized tissue sections, strong spatial organization is not expected. Instead, the clusters primarily represent different cellular states such as stages of the cell cycle.

---

# Conclusion

This tutorial demonstrated a complete workflow for analyzing and visualizing spatial transcriptomics data using Scanpy. The workflow included dataset loading, quality control, filtering, normalization, dimensionality reduction, clustering, marker gene identification, and spatial visualization. Both 10x Genomics Visium and MERFISH spatial transcriptomics technologies were explored.

---
