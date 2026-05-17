# 🧬 Analyze Visium H&E Data using Squidpy
---

# Introduction

Spatial transcriptomics is a modern technology that combines gene expression analysis with spatial tissue organization. Unlike traditional RNA sequencing, spatial transcriptomics preserves the exact location of gene expression inside tissue sections. This enables researchers to study cellular organization, tissue architecture, and molecular interactions directly within biological tissues.

This tutorial demonstrates the analysis of Visium H&E spatial transcriptomics data using the Squidpy library. The dataset consists of a coronal section of the mouse brain generated using 10x Genomics Visium technology. A preprocessed dataset is provided in AnnData format together with the histology image stored in a Squidpy ImageContainer object.

The workflow combines transcriptomics analysis, image feature extraction, graph-based spatial analysis, ligand-receptor interaction analysis, and spatially variable gene detection.

The tutorial includes:

- Loading Visium H&E datasets
- Spatial cluster visualization
- Extraction of image features
- Feature-based clustering
- Spatial graph analysis
- Neighborhood enrichment analysis
- Cluster co-occurrence analysis
- Ligand-receptor interaction analysis
- Spatial autocorrelation analysis using Moran’s I

---

# Importing Required Libraries

```python
%matplotlib inline
```

This command ensures that all figures and plots are displayed directly inside the Jupyter Notebook environment.

---

```python
import numpy as np
```

NumPy is imported for numerical operations, array processing, and mathematical computations used during spatial analysis.

---

```python
import pandas as pd
```

Pandas is used for handling feature matrices, metadata tables, and tabular outputs generated during analysis.

---

```python
import anndata as ad
```

AnnData provides the main data structure used for storing transcriptomics datasets including gene expression matrices, metadata, embeddings, and image-derived features.

---

```python
import scanpy as sc
```

Scanpy is used for transcriptomics preprocessing, PCA, neighborhood graph construction, clustering, and statistical analysis.

---

```python
import squidpy as sq
```

Squidpy is the primary library used for spatial transcriptomics analysis, image processing, graph analysis, and spatial visualization.

---

# Displaying Package Versions

```python
sc.logging.print_header()

print(f"squidpy=={sq.__version__}")
```

The `print_header()` function displays the versions of Scanpy and all related dependencies installed in the environment. The second command prints the installed Squidpy version.

The output confirms that the analysis uses:

- Scanpy 1.9.2
- Squidpy 1.2.3
- AnnData 0.8.0
- NumPy 1.23.5
- Pandas 1.5.1
- SciPy 1.9.3

Recording software versions is important for reproducibility because different library versions may produce different results.

---

# Loading the Visium H&E Dataset

```python
img = sq.datasets.visium_hne_image()

adata = sq.datasets.visium_hne_adata()
```

The first command loads the high-resolution H&E tissue image into an ImageContainer object called `img`.

The second command loads the corresponding transcriptomics dataset into an AnnData object called `adata`.

The dataset contains:

- Gene expression counts
- Spatial coordinates
- Cluster annotations
- Tissue metadata
- Histology image

The tissue section represents a coronal mouse brain slice with multiple anatomically distinct regions.

---

# Visualizing Spatial Clusters

```python
sq.pl.spatial_scatter(
    adata,
    color="cluster"
)
```
<img width="651" height="218" alt="gege1" src="https://github.com/user-attachments/assets/6971d81a-4805-4441-bae0-c846e658600d" />

This command visualizes transcriptional clusters directly in spatial coordinates.Each Visium spot is colored according to its cluster identity. The resulting figure shows how different molecularly distinct tissue regions are spatially organized across the mouse brain section.

Several brain regions become clearly visible, including:

- Hippocampus
- Cortex
- Fiber tracts
- Pyramidal layers

Neighboring spots often belong to similar clusters, demonstrating strong spatial organization within the tissue.The spatial distribution confirms that gene expression patterns closely reflect anatomical brain structures.

---

#  Importance of Image Features

Visium datasets contain high-resolution tissue images captured during sample preparation. These images contain valuable biological information that complements transcriptomics data.

Image features can describe:

- Tissue morphology
- Texture patterns
- Intensity distributions
- Structural organization
- Cellular density

Some morphological features correlate with gene expression patterns, while others provide additional biological insights that are not directly observable from RNA counts alone.

By integrating image features with transcriptomics data, tissue characterization becomes more detailed and biologically informative.

---

#  Extracting Multi-Scale Summary Features

```python
for scale in [1.0, 2.0]:

    feature_name = f"features_summary_scale{scale}"

    sq.im.calculate_image_features(
        adata,
        img.compute(),
        features="summary",
        key_added=feature_name,
        n_jobs=4,
        scale=scale,
    )
```

This loop extracts image summary features at two different spatial scales.

The scale parameter controls how much surrounding tissue context is included during feature extraction.

- `scale=1.0` extracts features mainly from the Visium spot region
- `scale=2.0` includes more neighboring tissue context

The extracted summary features include measurements such as:

- Mean intensity
- Standard deviation
- Statistical image properties

The features are stored separately inside:

```python
adata.obsm[feature_name]
```

Using multiple scales allows the analysis to capture both local tissue properties and broader tissue architecture.

---

# Combining Feature Tables

```python
adata.obsm["features"] = pd.concat(
    [
        adata.obsm[f]
        for f in adata.obsm.keys()
        if "features_summary" in f
    ],
    axis="columns",
)
```

This command combines all extracted feature tables into one large feature matrix.

The combined matrix contains image features from multiple spatial scales, allowing richer representation of tissue morphology.

---

```python
adata.obsm["features"].columns = ad.utils.make_index_unique(
    adata.obsm["features"].columns
)
```

This command ensures that all feature names are unique after merging multiple feature tables.

Unique feature names prevent conflicts during downstream clustering and analysis.

---

# Feature-Based Clustering Function

```python
def cluster_features(
    features: pd.DataFrame,
    like=None
) -> pd.Series:
```

This function clusters Visium spots based on image-derived features rather than gene expression profiles.

---

```python
if like is not None:
    features = features.filter(like=like)
```

This filters the feature matrix according to specific feature categories.

---

```python
adata = ad.AnnData(features)
```

A temporary AnnData object is created from the image feature matrix.

---

```python
sc.pp.scale(adata)
```

Feature values are standardized so that all features contribute equally during dimensionality reduction.

Scaling is important because image features may exist on very different numerical ranges.

---

```python
sc.pp.pca(
    adata,
    n_comps=min(10, features.shape[1] - 1)
)
```

Principal Component Analysis reduces dimensionality while preserving major variation in the image features.

---

```python
sc.pp.neighbors(adata)
```

A nearest-neighbor graph is constructed based on feature similarity between Visium spots.

---

```python
sc.tl.leiden(adata)
```

Leiden clustering groups Visium spots according to similarities in image morphology.

---

# Generating Feature-Based Clusters

```python
adata.obs["features_cluster"] = cluster_features(
    adata.obsm["features"],
    like="summary"
)
```

This command generates image-feature-based clusters and stores them inside:

```python
adata.obs["features_cluster"]
```

The clustering is based entirely on histology image morphology rather than transcriptomics data.

---

# Comparing Image Clusters and Gene Clusters

```python
sq.pl.spatial_scatter(
    adata,
    color=[
        "features_cluster",
        "cluster"
    ]
)
```
<img width="1474" height="434" alt="gege2" src="https://github.com/user-attachments/assets/782c4b97-ebe9-4f60-a620-dd8274a17d50" />

This visualization compares:
- Image-feature clusters
- Gene-expression clusters

Several regions show strong agreement between both clustering approaches.
For example:
- Fiber tract regions appear similarly clustered
- Hippocampal structures are partially recapitulated

However, some differences are also observed.
The gene-expression clusters capture layered cortical organization more strongly, while image-feature clusters identify broader morphological regions within the cortex.This demonstrates that image morphology and gene expression provide complementary biological information.Image-derived features capture structural tissue properties that may not be fully represented by transcriptomics alone.

---

# Spatial Statistics and Graph Analysis

Spatial transcriptomics data naturally contains neighborhood relationships between spots. Graph-based spatial analysis helps identify how clusters interact spatially within tissues.

Squidpy provides several graph-based methods to study tissue organization and cellular interactions.

---

# Building the Spatial Neighbor Graph

```python
sq.gr.spatial_neighbors(adata)
```

This command constructs a spatial connectivity graph.Each Visium spot becomes connected to nearby neighboring spots according to spatial coordinates.The resulting graph is stored inside the AnnData object and forms the basis for downstream spatial analyses.

---

#  Neighborhood Enrichment Analysis

```python
sq.gr.nhood_enrichment(
    adata,
    cluster_key="cluster"
)
```

This function calculates neighborhood enrichment scores between clusters.The method determines whether certain cluster pairs tend to occur near each other more often than expected by chance.The analysis uses permutation testing to estimate statistical significance.Clusters that frequently neighbor each other receive high enrichment scores, while clusters that rarely interact become depleted.

---

# Visualizing Neighborhood Enrichment

```python
sq.pl.nhood_enrichment(
    adata,
    cluster_key="cluster"
)
```
<img width="700" height="600" alt="gege3" src="https://github.com/user-attachments/assets/acbdd819-941c-4c59-83f4-a81d15c42225" />

The resulting heatmap visualizes enrichment relationships between tissue regions.The hippocampal clusters show strong enrichment relationships with pyramidal layers.
In particular:

- Hippocampus
- Pyramidal_layer
- Pyramidal_layer_dentate_gyrus

frequently occur near each other spatially.This result matches known anatomical organization of the mouse brain.Strong enrichment indicates biologically meaningful spatial interactions between tissue regions.

---

#  Co-Occurrence Analysis

Co-occurrence analysis evaluates how often clusters appear together across increasing spatial distances.Unlike neighborhood enrichment, this analysis operates directly on original spatial coordinates instead of graph connections.The co-occurrence score is defined as:

.. math:: \frac{p(exp|cond)}{p(exp)}

where 𝑝⁡(𝑒⁢𝑥⁢𝑝|𝑐⁢𝑜⁢𝑛⁢𝑑) is the conditional probability of observing a cluster 𝑒𝑥𝑝 conditioned on the presence of a cluster 𝑐⁢𝑜⁢𝑛⁢𝑑, whereas 𝑝⁡(𝑒⁢𝑥⁢𝑝) is the probability of observing 𝑒𝑥𝑝 in the radius size of interest. The score is computed across increasing radii size around each observation (i.e. spots here) in the tissue.Higher values indicate stronger spatial co-occurrence.

---

# Calculating Co-Occurrence Scores

```python
sq.gr.co_occurrence(
    adata,
    cluster_key="cluster"
)
```

This command computes co-occurrence scores for all cluster combinations across multiple spatial distances.

---

# Visualizing Co-Occurrence Patterns

```python
sq.pl.co_occurrence(
    adata,
    cluster_key="cluster",
    clusters="Hippocampus",
    figsize=(8, 4),
)
```
<img width="811" height="411" alt="gege4" src="https://github.com/user-attachments/assets/34412b4d-e98d-408a-9fe4-5172c314020c" />

This plot visualizes co-occurrence relationships involving the Hippocampus cluster.The curves show how likely other clusters are to occur near Hippocampus spots across increasing spatial radii.The Pyramidal_layer cluster strongly co-occurs with the Hippocampus at short distances.This confirms the anatomical proximity of these tissue structures within the mouse brain.The analysis demonstrates how tissue regions interact spatially beyond simple clustering.

---

# Ligand-Receptor Interaction Analysis

Spatial transcriptomics not only reveals tissue organization but also enables investigation of potential cellular communication.Ligand-receptor analysis identifies genes that may mediate signaling interactions between neighboring tissue regions.Squidpy implements a fast version of CellPhoneDB and integrates interaction databases from Omnipath.

---

# Calculating Ligand-Receptor Interactions

```python
sq.gr.ligrec(
    adata,
    n_perms=100,
    cluster_key="cluster",
)
```

This function computes statistically significant ligand-receptor interactions between all cluster pairs.Permutation testing is used to estimate interaction significance.

The analysis evaluates:

- Ligand expression
- Receptor expression
- Cluster-specific interactions
- Statistical enrichment

---

# Visualizing Ligand-Receptor Interactions

```python
sq.pl.ligrec(
    adata,
    cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=[
        "Pyramidal_layer",
        "Pyramidal_layer_dentate_gyrus"
    ],
    means_range=(3, np.inf),
    alpha=1e-4,
    swap_axes=True,
)
```

The resulting dot plot visualizes candidate ligand-receptor interactions between:

- Hippocampus
- Pyramidal_layer
- Pyramidal_layer_dentate_gyrus

Large dots indicate stronger interaction strength, while color intensity reflects expression levels.The analysis identifies signaling molecules that may contribute to communication between neighboring brain regions.These interactions represent candidate molecular pathways involved in tissue organization and neuronal signaling.The results provide biologically meaningful hypotheses for future experimental validation.

---

# Spatially Variable Gene Analysis using Moran’s I

Some genes exhibit strong spatial expression patterns across tissues.Spatial autocorrelation analysis identifies genes whose expression is spatially structured rather than randomly distributed.Squidpy implements Moran’s I statistic for detecting spatially variable genes.Higher Moran’s I values indicate stronger spatial organization.

---

# Selecting Highly Variable Genes

```python
genes = adata[
    :,
    adata.var.highly_variable
].var_names.values[:1000]
```

This selects the first 1000 highly variable genes for spatial autocorrelation analysis.Highly variable genes are more likely to exhibit biologically meaningful spatial patterns.

---

# Calculating Moran’s I

```python
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
    n_perms=100,
    n_jobs=1,
)
```

This function computes Moran’s I statistics for each selected gene.
The analysis estimates:

- Spatial autocorrelation strength
- Statistical significance
- Permutation-based p-values
- Adjusted false discovery rates

The results are stored inside:

```python
adata.uns["moranI"]
```

---

# Top Spatially Variable Genes

```python
adata.uns["moranI"].head(10)
```

The top-ranked genes include:

- Olfm1
- Plp1
- Itpka
- Snap25
- Nnat
- Ppp3ca
- Chn1
- Mal
- Tmsb4x
- Cldn11

These genes exhibit strong spatial organization across the tissue.

High Moran’s I values indicate that neighboring spots tend to have similar expression levels for these genes.

---

# Visualizing Spatially Variable Genes

```python
sq.pl.spatial_scatter(
    adata,
    color=[
        "Olfm1",
        "Plp1",
        "Itpka",
        "cluster"
    ]
)
```
<img width="2591" height="600" alt="gege5" src="https://github.com/user-attachments/assets/6b1b32d1-be3f-4dc0-a0b6-fc72d3b69aae" />

This visualization displays spatial expression patterns for selected highly variable genes.The genes show distinct tissue-specific localization patterns.
Some genes are highly expressed in:
- Pyramidal layers
- Fiber tract regions
- Hippocampal structures
The expression patterns strongly correspond with anatomical tissue organization.

For example:

- `Plp1` shows enrichment in fiber tract regions
- `Olfm1` and `Itpka` display strong hippocampal localization

These spatial expression patterns demonstrate biologically meaningful tissue specialization.

---

# Conclusion

This tutorial demonstrated a complete workflow for analyzing Visium H&E spatial transcriptomics data using Squidpy.
The workflow included:

- Loading Visium H&E datasets
- Spatial cluster visualization
- Multi-scale image feature extraction
- Feature-based clustering
- Spatial neighbor graph construction
- Neighborhood enrichment analysis
- Co-occurrence analysis
- Ligand-receptor interaction analysis
- Spatial autocorrelation analysis using Moran’s I

The results demonstrated that combining transcriptomics data with spatial image analysis provides deeper biological insight into tissue organization, anatomical structure, and molecular interactions.

Graph-based spatial analyses revealed strong neighborhood relationships between hippocampal and pyramidal-layer regions. Ligand-receptor analysis identified candidate molecular communication pathways between tissue structures. Moran’s I analysis detected genes with highly structured spatial expression patterns corresponding to distinct anatomical regions.

---
