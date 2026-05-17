# 🧬 Analyze Xenium Data using SpatialData, Scanpy, and Squidpy

---

# Introduction

Xenium is a high-resolution spatial transcriptomics platform developed by 10x Genomics that enables direct imaging of RNA molecules inside tissue sections at single-cell resolution. Unlike traditional bulk RNA sequencing or even conventional single-cell RNA sequencing, Xenium preserves the spatial position of every transcript and every cell within the tissue architecture. This allows researchers to study not only gene expression but also the physical organization of cells and their interactions inside tissues.

Spatial transcriptomics is especially important in cancer biology because tumors contain many heterogeneous cell populations that interact with immune cells, stromal cells, and blood vessels. By preserving spatial information, Xenium makes it possible to investigate tumor microenvironments, immune infiltration patterns, tissue heterogeneity, and localized gene expression patterns within complex tissues.

This tutorial demonstrates a complete spatial transcriptomics workflow using a Human Lung Cancer Xenium dataset together with SpatialData, Scanpy, Squidpy, SpatialData-Plot, and Napari-SpatialData.

---

# Import Required Libraries

```python
import spatialdata as sd
from spatialdata_io import xenium

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import squidpy as sq
```

- `spatialdata` manages multimodal spatial datasets.
- `spatialdata_io` reads Xenium datasets and converts them into SpatialData objects.
- `matplotlib` creates scientific visualizations.
- `seaborn` provides statistical plotting tools.
- `scanpy` performs preprocessing and single-cell analysis.
- `squidpy` performs spatial transcriptomics analysis.

---

# Define Dataset Paths

```python
xenium_path = "./Xenium"

zarr_path = "./Xenium.zarr"
```

- `xenium_path` stores the location of the Xenium dataset.
- `zarr_path` defines where the processed Zarr dataset will be saved.

Zarr format improves loading speed and memory efficiency for large spatial datasets.

---

# Load Xenium Dataset

```python
sdata = xenium(xenium_path)
```

The Xenium reader imports the complete dataset and converts it into a SpatialData object.

The object contains:

- Tissue morphology images
- Transcript coordinates
- Segmentation masks
- Cell boundaries
- Gene expression matrices

This unified structure keeps imaging and molecular data connected throughout the analysis.

---

# Convert Dataset into Zarr Format

```python
sdata.write(zarr_path)
```

The dataset is saved into Zarr format for:

- Faster loading
- Efficient storage
- Better memory handling
- Scalable processing of large datasets

---

# Read Dataset from Zarr Store

```python
sdata = sd.read_zarr(zarr_path)
```

The processed dataset is loaded directly from the Zarr store instead of raw Xenium files.This reduces preprocessing time and improves reproducibility.

---

#Inspect the SpatialData Object

```python
sdata
```

The SpatialData object contains multiple biological layers:

- Images → tissue morphology
- Labels → segmentation masks
- Points → transcript coordinates
- Shapes → cell boundaries
- Tables → gene expression matrices

---

# Extract the AnnData Object

```python
adata = sdata.tables["table"]
```

The AnnData object stores transcriptomics information required for downstream analysis.

It contains:

- Gene expression matrix
- Cell metadata
- Spatial coordinates
- Clustering information
- Dimensionality reduction embeddings

---

# Inspect the AnnData Object

```python
adata
```

Important AnnData components:

- `n_obs` → number of cells
- `n_vars` → number of genes
- `.obs` → cell metadata
- `.var` → gene metadata
- `.obsm` → embeddings and spatial coordinates

---

# Access Spatial Coordinates

```python
adata.obsm["spatial"]
```

This matrix stores the x and y coordinates of every cell in the tissue.
These coordinates are used for:

- Spatial visualization
- Neighborhood analysis
- Spatial graph construction
- Tissue organization analysis

---

# Calculate Quality-Control Metrics

```python
sc.pp.calculate_qc_metrics(
    adata,
    percent_top=(10, 20, 50, 150),
    inplace=True
)
```

This step calculates quality-control (QC) metrics for every cell in the dataset.

Quality control is important because spatial transcriptomics datasets may contain:

- Damaged cells
- Empty droplets
- Segmentation errors
- Doublets (two cells counted as one)

The function computes several statistics:

- Total transcript counts per cell
- Number of detected genes
- Percentage contribution of highly expressed genes
- Sequencing complexity

### Parameters

- `adata`
  → AnnData object containing expression data.

- `percent_top=(10, 20, 50, 150)`
  → Calculates the percentage contribution of the top expressed genes.

Example:
- Top 10 genes
- Top 20 genes
- Top 50 genes
- Top 150 genes

This helps determine whether a few genes dominate total expression.

- `inplace=True`
  → Saves results directly inside the AnnData object.

---

# Calculate Control Probe Percentages

```python
cprobes = (
    adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100
)

cwords = (
    adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100
)

print(f"Negative DNA probe count % : {cprobes}")
print(f"Negative decoding count % : {cwords}")
```

This section measures technical background noise in the dataset.

### Control Probes

Control probes are artificial probes that should not bind to real biological transcripts.

If control probe counts are high:
- Experimental noise may be high
- Signal quality may be poor

### Control Codewords

Control codewords are artificial decoding signals.

They help estimate:
- False positive transcript detection
- Decoding errors

### Formula

```python
(control counts / total transcript counts) × 100
```

### Output

```python
Negative DNA probe count % : 0.0050
Negative decoding count % : 0.0025
```

Very low percentages indicate:
- Good sequencing quality
- Low technical noise
- Reliable transcript detection

---

# Visualize Quality-Control Distributions

```python
fig, axs = plt.subplots(1, 4, figsize=(15, 4))

axs[0].set_title("Total transcripts per cell")
sns.histplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)

axs[1].set_title("Unique transcripts per cell")
sns.histplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    ax=axs[1],
)

axs[2].set_title("Area of segmented cells")
sns.histplot(
    adata.obs["cell_area"],
    kde=False,
    ax=axs[2],
)

axs[3].set_title("Nucleus ratio")
sns.histplot(
    adata.obs["nucleus_area"] / adata.obs["cell_area"],
    kde=False,
    ax=axs[3],
)
```
<img width="1250" height="393" alt="ha1" src="https://github.com/user-attachments/assets/270fa197-9548-49b8-9d82-69aba460b6c3" />

This block visualizes important quality-control metrics using histograms.

---

## 1️⃣Total Transcripts per Cell

```python
adata.obs["total_counts"]
```

Measures:
- Total RNA molecules detected in each cell

Low counts may indicate:
- Dead cells
- Poor capture efficiency

Extremely high counts may indicate:
- Doublets
- Segmentation artifacts

---

## 2️⃣Unique Transcripts per Cell

```python
adata.obs["n_genes_by_counts"]
```

Measures:
- Number of unique genes detected per cell

Higher diversity generally indicates:
- Better cell quality
- Rich biological information

Very low values may indicate:
- Damaged cells

---

## 3️⃣ Area of Segmented Cells

```python
adata.obs["cell_area"]
```

Measures:
- Physical size of segmented cells

Helps identify:
- Segmentation errors
- Extremely small or unusually large cells

---

## 4️⃣ Nucleus Ratio

```python
adata.obs["nucleus_area"] / adata.obs["cell_area"]
```

Measures:
- Fraction of cell occupied by nucleus

Abnormal ratios may indicate:
- Incorrect segmentation
- Damaged cells

---

# Filter Low-Quality Cells and Genes

```python
sc.pp.filter_cells(adata, min_counts=10)

sc.pp.filter_genes(adata, min_cells=5)
```

This step removes low-quality biological data.

---

## Filter Cells

```python
min_counts=10
```

Cells with fewer than 10 transcripts are removed.

Reason:
- Very low transcript counts often represent noise.

---

## Filter Genes

```python
min_cells=5
```

Genes expressed in fewer than 5 cells are removed.

Reason:
- Rare genes contribute little useful information
- Reduces computational complexity

---

# ⚙️ Normalize and Preprocess Data

```python
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, inplace=True)

sc.pp.log1p(adata)

sc.pp.pca(adata)

sc.pp.neighbors(adata)

sc.tl.umap(adata)

sc.tl.leiden(adata)
```

---

# Store Raw Counts

```python
adata.layers["counts"] = adata.X.copy()
```

## Explanation

Saves original raw counts before normalization.

Useful because:
- Some analyses require raw counts later.

---

# Normalize Total Counts

```python
sc.pp.normalize_total(adata, inplace=True)
```

Normalization makes transcript counts comparable across cells.

Without normalization:
- Cells with higher sequencing depth dominate analysis.

Normalization adjusts all cells to similar total counts.

---

# Log Transformation

```python
sc.pp.log1p(adata)
```

Applies logarithmic transformation:


::contentReference[oaicite:0]{index=0}


Purpose:
- Reduces effect of extremely highly expressed genes
- Stabilizes variance
- Improves clustering

---

# Principal Component Analysis (PCA)

```python
sc.pp.pca(adata)
```

## Explanation

PCA reduces high-dimensional gene expression data into smaller dimensions.

Benefits:
- Faster computation
- Noise reduction
- Better visualization

The most important biological variation is preserved.

---

# Construct Neighbor Graph

```python
sc.pp.neighbors(adata)
```

## Explanation

Builds a graph connecting biologically similar cells.

Cells with similar expression profiles become neighbors.

This graph is used for:
- UMAP
- Clustering
- Spatial analysis

---

# Generate UMAP Embedding

```python
sc.tl.umap(adata)
```

UMAP creates a 2D visualization of cellular relationships.

Similar cells cluster together in low-dimensional space.

This helps visualize:
- Cell populations
- Biological heterogeneity
- Cluster structure

---

# Leiden Clustering

```python
sc.tl.leiden(adata)
```

Leiden algorithm identifies cell clusters.

Each cluster may represent:
- Cell type
- Cell state
- Biological niche

Cluster labels are stored in:

```python
adata.obs["leiden"]
```

---

# Visualize UMAP Clusters

```python
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)
```
<img width="2280" height="431" alt="ha2" src="https://github.com/user-attachments/assets/a08b233b-237f-4d69-9f20-695ad1439f39" />

This creates UMAP plots colored by:

- Total transcript counts
- Number of genes
- Leiden clusters

Purpose:
- Inspect clustering quality
- Identify batch effects
- Observe biological structure

---

# Visualize Spatial Clusters

```python
sq.pl.spatial_scatter(
    adata,
    library_id="spatial",
    shape=None,
    color=[
        "leiden",
    ],
    wspace=0.4,
)
```
<img width="651" height="252" alt="ha3" src="https://github.com/user-attachments/assets/121cc341-75dc-4948-ac81-8257f145be8e" />

Plots clusters directly on tissue coordinates.

Purpose:
- Observe tissue organization
- Detect spatially localized cell populations
- Study tumor microenvironment structure

---

# Build Spatial Neighbor Graph

```python
sq.gr.spatial_neighbors(
    adata,
    coord_type="generic",
    delaunay=True
)
```
Constructs a spatial connectivity graph based on physical cell positions.

### Delaunay Triangulation

```python
delaunay=True
```

Creates connections between nearby cells.

Used for:
- Spatial statistics
- Neighborhood analysis
- Interaction studies

---

# Compute Centrality Scores

```python
sq.gr.centrality_scores(
    adata,
    cluster_key="leiden"
)
```

Calculates network centrality metrics for each Leiden cluster.

---

## Closeness Centrality

Measures:
- How close a cluster is to all other clusters

Higher values indicate:
- Strong tissue integration

---

## Degree Centrality

Measures:
- Number of neighboring interactions

Higher values indicate:
- Highly connected clusters

---

## Clustering Coefficient

Measures:
- How tightly neighboring cells cluster together

Higher values indicate:
- Dense local organization

---

# Visualize Centrality Scores

```python
sq.pl.centrality_scores(
    adata,
    cluster_key="leiden",
    figsize=(16, 5)
)
```
<img width="1611" height="511" alt="ha4" src="https://github.com/user-attachments/assets/842e7c7d-f9da-4300-bc27-95f91c54d899" />


Displays spatial interaction statistics for each cluster.

Helps identify:
- Dominant tissue regions
- Highly interactive cell populations
- Spatially organized communities

---

# Create Subsample Dataset

```python
sdata.tables["subsample"] = sc.pp.subsample(
    adata,
    fraction=0.5,
    copy=True
)

adata_subsample = sdata.tables["subsample"]
```
Creates a smaller dataset using 50% of cells.

Purpose:
- Faster computation
- Reduced memory usage
- Efficient spatial analysis

---

# Compute Co-Occurrence Probability

```python
sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)
```

Measures how frequently different cell clusters occur near each other.

High co-occurrence suggests:
- Biological interaction
- Shared tissue niches
- Spatial association

---

# Visualize Co-Occurrence Patterns

```python
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(10, 10),
)
```

Visualizes interaction probability around cluster 12.

Helps identify:
- Neighboring cell populations
- Tumor-immune interactions
- Tissue organization patterns

---

# Spatial Scatter Plot

```python
sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
)
```
<img width="900" height="1011" alt="ha6" src="https://github.com/user-attachments/assets/c5ecabae-ec50-4ae1-8cd4-fb5ce4a9a266" />
<img width="651" height="252" alt="ha7" src="https://github.com/user-attachments/assets/9eea9696-db78-40f4-b3f6-2c158136696e" />

Displays Leiden clusters across physical tissue coordinates.

Purpose:
- Visualize tissue architecture
- Observe spatial cell arrangements

---

# Neighborhood Enrichment Analysis

```python
sq.gr.nhood_enrichment(
    adata,
    cluster_key="leiden"
)
```

Measures whether specific clusters appear near each other more often than expected by chance.

Uses:
- Permutation testing
- Z-score calculation

High enrichment:
- Strong biological association

Negative enrichment:
- Spatial separation

---

# Visualize Neighborhood Enrichment

```python
fig, ax = plt.subplots(1, 2, figsize=(13, 7))

sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)

sq.pl.spatial_scatter(
    adata_subsample,
    color="leiden",
    shape=None,
    size=2,
    ax=ax[1]
)
```
Creates:
- Heatmap of neighborhood enrichment
- Spatial cluster visualization
<img width="1224" height="494" alt="haha8" src="https://github.com/user-attachments/assets/ae1a3d59-53ef-46c1-ac3e-ff1d03756e84" />

This helps compare:
- Statistical enrichment
- Actual tissue arrangement

---

# Compute Moran’s I Spatial Autocorrelation

```python
sq.gr.spatial_neighbors(
    adata_subsample,
    coord_type="generic",
    delaunay=True
)

sq.gr.spatial_autocorr(
    adata_subsample,
    mode="moran",
    n_perms=100,
    n_jobs=1,
)
```


Moran’s I measures spatial autocorrelation.

Purpose:
- Detect spatially variable genes

Interpretation:

- Positive Moran’s I
  → clustered expression

- Negative Moran’s I
  → dispersed expression

- Near zero
  → random distribution

---

# View Top Spatially Variable Genes

```python
adata_subsample.uns["moranI"].head(10)
```

## Explanation

Displays genes with strongest spatial organization.

Examples:
- AREG
- MET
- EPCAM

These genes show non-random spatial expression patterns.

---

#  Visualize Spatial Gene Expression

```python
sq.pl.spatial_scatter(
    adata_subsample,
    library_id="spatial",
    color=[
        "AREG",
        "MET",
    ],
    shape=None,
    size=2,
    img=False,
)
```
Displays spatial expression of genes across tissue.
<img width="1197" height="398" alt="haha1" src="https://github.com/user-attachments/assets/7f5af0bb-4d3c-4f4a-a806-ad98b0539195" />

Purpose:
- Identify tumor regions
- Detect localized expression
- Study tissue heterogeneity

---

# Visualize Expression on Morphology Images

```python
import spatialdata_plot

gene_name = ["AREG", "MET"]

for name in gene_name:

    sdata.pl.render_images("morphology_focus").pl.render_shapes(
        "cell_circles",
        color=name,
        table_name="table",
        use_raw=False,
    ).pl.show(
        title=f"{name} expression over Morphology image",
        coordinate_systems="global",
        figsize=(10, 5),
    )
```

Overlays gene expression directly on tissue morphology images.
<img width="1011" height="511" alt="haha2" src="https://github.com/user-attachments/assets/b83402fb-ffc0-49bd-8c0d-76e5f17e729f" />
<img width="1011" height="511" alt="haha3" src="https://github.com/user-attachments/assets/5aa4184d-8c80-48ab-91af-5d30376fd700" />

Purpose:
- Compare molecular expression with tissue structure
- Identify tumor regions
- Observe spatial localization

---

# Interactive Visualization with Napari

```python
from napari_spatialdata import Interactive

Interactive(sdata)
```
<img width="2560" height="1377" alt="haha4" src="https://github.com/user-attachments/assets/98c26247-c047-4641-80a5-91ceae87b650" />

Launches an interactive graphical interface for spatial exploration.

Users can:
- Zoom into tissue regions
- Explore gene expression
- Visualize segmentation masks
- Interactively inspect spatial organization

Napari provides high-resolution visualization for large spatial transcriptomics datasets.

---
# Conclusion

This project demonstrated a complete spatial transcriptomics analysis workflow using Xenium data together with SpatialData, Scanpy, and Squidpy. The workflow integrated molecular transcriptomics information with spatial tissue organization to study cellular heterogeneity and spatial interactions within a human lung cancer tissue sample.The analysis began with importing the Xenium dataset and converting it into Zarr format for efficient storage and faster computation. The SpatialData framework enabled unified management of tissue images, transcript coordinates, segmentation masks, cell boundaries, and gene expression matrices inside a single structured object.
Spatial analysis using Squidpy revealed important tissue organization patterns. Spatial neighbor graphs and centrality scores identified highly connected cellular communities and spatially influential clusters. Co-occurrence analysis and neighborhood enrichment analysis further highlighted biologically meaningful spatial interactions between cell populations.


---
