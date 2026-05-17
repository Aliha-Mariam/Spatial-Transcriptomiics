# 🧬 Analyze Visium Fluorescence Data using Squidpy



---

# Introduction

Spatial transcriptomics combines gene expression analysis with spatial information from tissue sections. In this tutorial, Visium fluorescence data is analyzed using the Squidpy library. The workflow demonstrates how image analysis features can be integrated with transcriptomics data to better understand tissue organization and cellular morphology.

The dataset used in this tutorial is a Visium slide containing a coronal section of the mouse brain. The original dataset is publicly available from the 10x Genomics dataset portal. A preprocessed version of the dataset is used here, where clusters have already been annotated. The tutorial combines fluorescence image analysis, segmentation, feature extraction, clustering, and spatial visualization.

The fluorescence image contains three biological channels:

- DAPI → stains DNA and nuclei
- anti-NEUN → labels neurons
- anti-GFAP → labels glial cells

The tutorial demonstrates how to:

- Load Visium fluorescence datasets
- Visualize clusters spatially
- Analyze fluorescence channels
- Segment tissue images
- Extract image-based features
- Cluster spots based on image morphology
- Compare image-derived clusters with gene-expression clusters

---

# Importing Required Libraries

```python
%matplotlib inline
```

This command ensures that plots and visualizations appear directly inside the Jupyter Notebook instead of opening in external windows.

---

```python
import pandas as pd
```

Pandas is imported for handling tabular data structures such as feature matrices and metadata tables.

---

```python
import matplotlib.pyplot as plt
```

Matplotlib is used for plotting figures, images, and segmentation results.

---

```python
import anndata as ad
```

AnnData is the main data structure used for storing transcriptomics datasets including gene expression matrices, metadata, embeddings, and image features.

---

```python
import scanpy as sc
```

Scanpy is used for preprocessing, clustering, PCA, neighborhood graph construction, and transcriptomics analysis.

---

```python
import squidpy as sq
```

Squidpy is the primary library used for spatial transcriptomics image analysis and spatial visualization.

---

# Display Package Versions

```python
sc.logging.print_header()

print(f"squidpy=={sq.__version__}")
```

The `print_header()` function displays versions of Scanpy and all related dependencies. This helps ensure reproducibility of the workflow. The second command specifically prints the installed version of Squidpy.

The output confirms that the environment contains:

- Scanpy 1.9.2
- Squidpy 1.2.3
- AnnData 0.8.0
- NumPy 1.23.5
- Pandas 1.5.1
- SciPy 1.9.3

These package versions are important because spatial transcriptomics workflows can behave differently across software versions.

---

# Loading the Dataset

```python
img = sq.datasets.visium_fluo_image_crop()

adata = sq.datasets.visium_fluo_adata_crop()
```

The first command loads the fluorescence tissue image into an `ImageContainer` object named `img`. This object stores the microscopy image and all image-processing layers.

The second command loads the transcriptomics dataset into an AnnData object named `adata`. This dataset contains:

- Gene expression counts
- Spatial coordinates
- Cluster annotations
- Metadata for Visium spots

The provided dataset is a cropped region of a mouse brain section. A smaller crop is used to reduce computational time during image processing and feature extraction.

---

# Visualizing Spatial Clusters

```python
sq.pl.spatial_scatter(
    adata,
    color="cluster"
)
```

This command visualizes the annotated clusters directly in spatial coordinates. Each Visium spot is colored according to its assigned cluster identity.


<img width="644" height="491" alt="mage1" src="https://github.com/user-attachments/assets/1ebfc391-0b71-4cce-bac4-cdc5fd2ff45b" />

The resulting figure shows how different transcriptional clusters are organized spatially within the mouse brain tissue. Neighboring spots often belong to similar clusters, indicating biologically organized tissue regions.The hippocampus, cortex, and fiber tract regions can already be distinguished based on the cluster distribution. Spatial visualization allows researchers to connect molecular expression patterns with anatomical tissue structures.

---

#  Visualizing Fluorescence Channels

```python
img.show(channelwise=True)
```
<img width="790" height="287" alt="mage2" src="https://github.com/user-attachments/assets/477ba125-77e2-4b3e-8591-4540b86a8a5f" />

This command visualizes each fluorescence channel separately.

The fluorescence image contains three channels:

| Channel | Biological Marker | Function |
|---|---|---|
| DAPI | DNA stain | Labels nuclei |
| anti-NEUN | Neuronal marker | Labels neurons |
| anti-GFAP | Glial marker | Labels glial cells |

The resulting plots show the spatial distribution of each biological marker across the tissue.The DAPI channel highlights cell nuclei throughout the tissue. The anti-NEUN channel shows neuron-rich regions, while the anti-GFAP channel identifies areas enriched with glial cells.These fluorescence patterns provide complementary biological information that cannot always be observed directly from gene expression data alone.

---

# Importance of Image Features

Visium datasets contain high-resolution tissue images alongside gene expression measurements. By extracting image features from these images, additional biological information can be integrated with transcriptomics data.

Image features can capture:

- Tissue morphology
- Cell density
- Texture patterns
- Fluorescence intensity
- Structural organization

Some image features may overlap with gene expression patterns, while others provide entirely complementary information.

For example:

- Regions with different cell types may have distinct morphologies
- Cell segmentation can estimate cell numbers
- Fluorescence intensity can identify enriched cell populations

This integration improves tissue characterization beyond transcriptomics alone.

---

# Image Preprocessing and Segmentation

```python
sq.im.process(
    img=img,
    layer="image",
    method="smooth",
)
```

This command preprocesses the fluorescence image by applying smoothing. Smoothing reduces image noise and improves segmentation quality.

The processed image is stored in a new image layer called `"image_smooth"`.

---

```python
sq.im.segment(
    img=img,
    layer="image_smooth",
    method="watershed",
    channel=0,
    chunks=1000
)
```

This command performs image segmentation using the watershed algorithm.

The segmentation uses:

- The smoothed image layer
- Channel 0 (DAPI channel)
- Chunk processing for memory efficiency

The watershed algorithm identifies individual nuclei by separating neighboring objects based on intensity gradients.

The segmentation result is automatically stored in:

```python
img["segmented_watershed"]
```

Each segmented object receives a unique integer label.

---

#  Visualizing Segmentation Results

```python
fig, ax = plt.subplots(1, 2)

img_crop = img.crop_corner(
    2000,
    2000,
    size=500
)

img_crop.show(
    layer="image",
    channel=0,
    ax=ax[0]
)

img_crop.show(
    layer="segmented_watershed",
    channel=0,
    ax=ax[1],
)
```
<img width="515" height="265" alt="mage3" src="https://github.com/user-attachments/assets/172c713b-a847-477c-86f6-dd82b4ab5259" />

The first plot displays the original DAPI fluorescence image for a cropped tissue region.The second plot visualizes the segmentation output produced by the watershed algorithm.In the segmentation image, each nucleus is represented as a separate labeled object. Proper segmentation indicates that nuclei boundaries are accurately detected and neighboring nuclei are successfully separated.This segmentation forms the basis for extracting morphology-based image features.

---

#  Extracting Segmentation Features

```python
features_kwargs = {
    "segmentation": {
        "label_layer": "segmented_watershed"
    }
}
```

This dictionary specifies which segmentation layer should be used during feature extraction.

---

```python
sq.im.calculate_image_features(
    adata,
    img,
    features="segmentation",
    layer="image",
    key_added="features_segmentation",
    n_jobs=1,
    features_kwargs=features_kwargs,
)
```

This command extracts segmentation-based image features for every Visium spot.

The extracted features include:

- Number of segmented cells
- Morphological measurements
- Mean fluorescence intensities
- Channel-specific statistics

The results are stored inside:

```python
adata.obsm["features_segmentation"]
```

These image-derived features can later be integrated with transcriptomics data.

---

# Visualizing Segmentation Features

```python
sq.pl.spatial_scatter(
    sq.pl.extract(
        adata,
        "features_segmentation"
    ),
    color=[
        "segmentation_label",
        "cluster",
        "segmentation_ch-0_mean_intensity_mean",
        "segmentation_ch-1_mean_intensity_mean",
    ],
    frameon=False,
    ncols=2,
)
```
<img width="1089" height="827" alt="mage4" src="https://github.com/user-attachments/assets/4a4df78f-6db8-46d6-afed-eb4460679bf1" />

This visualization compares segmentation-derived image features with transcriptomics clusters.The first panel shows the number of segmented cells within each Visium spot. Cell-rich regions appear with higher values.The hippocampal pyramidal layer contains many densely packed cells, resulting in higher segmentation counts compared to surrounding tissue areas.The second panel displays the original gene-expression cluster annotation.The third panel visualizes mean intensity for fluorescence channel 0. High intensities indicate stronger DAPI staining and greater nuclear density.The fourth panel visualizes mean intensity for channel 1, corresponding to anti-NEUN staining. Regions with stronger signal contain more neurons.The Cortex_1 and Cortex_3 regions exhibit higher anti-NEUN intensity, indicating neuron-rich tissue areas.Fiber tracts and lateral ventricles show stronger anti-GFAP signal, indicating enrichment of glial cells.These fluorescence-derived measurements provide biologically meaningful information that complements gene-expression clustering.

---

#  Extracting Summary, Histogram, and Texture Features

```python
params = {
    "features_orig": {
        "features": [
            "summary",
            "texture",
            "histogram"
        ],
        "scale": 1.0,
        "mask_circle": True,
    },

    "features_context": {
        "features": [
            "summary",
            "histogram"
        ],
        "scale": 1.0,
    },

    "features_lowres": {
        "features": [
            "summary",
            "histogram"
        ],
        "scale": 0.25,
    },
}
```

This dictionary defines multiple feature extraction strategies.

The extracted features include:

| Feature Type | Biological Meaning |
|---|---|
| Summary | Mean intensity and statistical properties |
| Histogram | Distribution of pixel intensities |
| Texture | Spatial texture and structural organization |

Different scales and crop sizes provide both local and contextual tissue information.

---

#  Calculating Image Features

```python
for feature_name, cur_params in params.items():

    sq.im.calculate_image_features(
        adata,
        img,
        layer="image",
        key_added=feature_name,
        n_jobs=1,
        **cur_params
    )
```

This loop calculates image features using all parameter combinations defined earlier.

The extracted features are stored separately inside:

```python
adata.obsm[feature_name]
```

Each Visium spot receives a large set of image-derived numerical features describing morphology and fluorescence properties.

---

#  Combining Feature Tables

```python
adata.obsm["features"] = pd.concat(
    [
        adata.obsm[f]
        for f in params.keys()
    ],
    axis="columns"
)
```

This command combines all extracted feature tables into one large feature matrix.

The combined matrix contains summary, histogram, and texture features extracted at multiple scales.

---

```python
adata.obsm["features"].columns = ad.utils.make_index_unique(
    adata.obsm["features"].columns
)
```

This command ensures that all feature names are unique after combining multiple feature tables.

Unique column names prevent conflicts during downstream clustering and analysis.

---

#  Clustering Image Features

```python
def cluster_features(
    features: pd.DataFrame,
    like=None
):
```

This function clusters Visium spots based on image-derived features rather than gene expression data.

---

```python
if like is not None:
    features = features.filter(like=like)
```

This filters the feature table based on specific feature types such as summary, histogram, or texture.

---

```python
adata = ad.AnnData(features)
```

A temporary AnnData object is created using the image features.

---

```python
sc.pp.scale(adata)
```

Feature values are scaled so that all features contribute equally during PCA and clustering.

---

```python
sc.pp.pca(
    adata,
    n_comps=min(10, features.shape[1] - 1)
)
```

PCA reduces the dimensionality of the image feature matrix.

---

```python
sc.pp.neighbors(adata)
```

A nearest-neighbor graph is constructed based on feature similarity.

---

```python
sc.tl.leiden(adata)
```

Leiden clustering groups Visium spots according to image morphology.

---

# Creating Feature-Based Clusters

```python
adata.obs["features_summary_cluster"] = cluster_features(
    adata.obsm["features"],
    like="summary"
)

adata.obs["features_histogram_cluster"] = cluster_features(
    adata.obsm["features"],
    like="histogram"
)

adata.obs["features_texture_cluster"] = cluster_features(
    adata.obsm["features"],
    like="texture"
)
```

This section creates separate cluster annotations based on:

- Summary features
- Histogram features
- Texture features

Each clustering method captures different image characteristics.

---

# Visualizing Feature-Based Clusters

```python
sc.set_figure_params(
    facecolor="white",
    figsize=(8, 8)
)

sq.pl.spatial_scatter(
    adata,
    color=[
        "features_summary_cluster",
        "features_histogram_cluster",
        "features_texture_cluster",
        "cluster",
    ],
    ncols=3,
)
```
<img width="3440" height="2213" alt="mage5" src="https://github.com/user-attachments/assets/6766f445-1b58-406d-8512-511b61932fed" />

This visualization compares image-derived clusters with transcriptomics-based clusters.The feature-space clusters are highly spatially coherent, meaning neighboring spots tend to belong to similar image-derived groups.The hippocampus is subdivided into multiple fine-grained structural regions by the image features. This level of structural detail is greater than what was observed using gene-expression clustering alone.Texture-based clusters capture tissue organization patterns and layered structures within the cortex.Histogram-based clusters separate regions according to fluorescence intensity distributions.Summary-feature clusters identify broader morphological differences between tissue regions.Compared with gene-space clustering, image features provide additional spatial and structural information that improves tissue characterization.

---

#  Conclusion

This tutorial demonstrated how Squidpy can integrate fluorescence image analysis with spatial transcriptomics data from Visium experiments.

The workflow included:

- Loading fluorescence Visium datasets
- Visualizing transcriptomics clusters
- Analyzing fluorescence channels
- Image preprocessing
- Nuclei segmentation
- Extracting segmentation features
- Calculating summary, histogram, and texture features
- Clustering tissue regions using image morphology
- Comparing image-based clusters with gene-expression clusters

The results show that fluorescence image features provide highly informative biological and structural information that complements transcriptomics data. Image-based clustering captures fine tissue structures, neuron-rich regions, glial-cell enrichment, and cortical organization that may not be fully resolved using gene-expression data alone.

---
