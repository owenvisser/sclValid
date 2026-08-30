# sclValid

`sclValid` is an R package for clustering and validating single-cell RNA sequencing data. It provides a unified workflow for applying multiple clustering methods, calculating internal, external, and stability-based validation measures, scaling validation results, and aggregating them into an overall ranking of clustering solutions.

The methodology implemented in this package is based on:

> O. Visser and S. Datta, “Integrating Multiple Clustering Techniques and Performance Measures via Ranking for scRNA-Seq Data,” *Statistics in Medicine* 44, no. 28-30 (2025): e70331. https://doi.org/10.1002/sim.70331

The package is built around `SingleCellExperiment` objects. Expression data are stored in the `"data"` assay, clustering assignments are stored in `colData`, and validation and ranking results are stored in `metadata`.

## Installation

The development version can be installed from GitHub with:

```r
install.packages("remotes")
remotes::install_github("owenvisser/sclValid")
```

Then load the package:

```r
library(sclValid)
```

## Basic workflow

`sclValid` expects a numeric matrix with features in rows and cells in columns.

A `SingleCellExperiment` object can be created with:

```r
sce <- make_sce(
  example_data,
  example_labels
)
```

The package includes a small example dataset derived from the single-cell RNA-seq data of Biase et al. The example contains 1,000 sampled features and 49 cells.

## Clustering

Several clustering procedures are available.

### pcaReduce

```r
sce <- run_pcareduce(
  sce,
  cluster_sizes = 2:4,
  method = "M"
)

sce <- run_pcareduce(
  sce,
  cluster_sizes = 2:4,
  method = "S"
)
```

### RaceID

```r
sce <- run_raceid(
  sce,
  cluster_sizes = 2:4
)
```

### PCA- and t-SNE-based clustering

```r
sce <- run_reduction_clustering(
  sce,
  reductions = c("PCA", "TSNE"),
  methods = c(
    "kmeans",
    "hierarchical",
    "louvain",
    "leiden"
  ),
  cluster_sizes = 2:4,
  knn = c(5, 10),
  resolution = c(0.5, 1)
)
```

Clustering assignments are added as columns in:

```r
SummarizedExperiment::colData(sce)
```

## Validation

Validation measures can be calculated with:

```r
sce <- run_validation(
  sce,
  measures = c(
    "AD",
    "ADM",
    "APN",
    "ARI",
    "BHI",
    "BSI",
    "CN",
    "DI",
    "IGP",
    "SW"
  )
)
```

The available measures are:

* **AD**: Average Distance
* **ADM**: Average Distance Between Means
* **APN**: Average Proportion of Non-Overlap
* **ARI**: Adjusted Rand Index
* **BHI**: Biological Homogeneity Index
* **BSI**: Biological Stability Index
* **CN**: Connectivity
* **DI**: Dunn Index
* **IGP**: In-Group Proportion
* **SW**: Silhouette Width

ARI, BHI, and BSI require known labels. AD, ADM, APN, and BSI use perturbed datasets and repeated reclustering to evaluate stability.

By default, stability-based validation displays a progress bar while suppressing detailed clustering output.

Raw validation results are stored in:

```r
S4Vectors::metadata(sce)$validation
```

## Scaling validation measures

Validation measures have different ranges and directions. They can be transformed to a common scale with:

```r
sce <- scale_validation(sce)
```

After scaling, larger values consistently indicate better clustering performance.

Scaled results are stored separately from the raw validation measures:

```r
S4Vectors::metadata(sce)$validation_scaled
```

## Rank aggregation

The scaled validation measures can be combined into an overall ranking of clustering solutions:

```r
sce <- run_rank_aggregation(sce)
```

Rank aggregation is performed in two stages. Clustering solutions are first compared within each cluster size, and the highest-ranked solution for each size is retained. These selected solutions are then aggregated to produce the final ranking.

The final ranking is stored in:

```r
S4Vectors::metadata(sce)$ranking
```

## Complete example

```r
library(sclValid)

sce <- make_sce(
  example_data,
  example_labels
)

sce <- run_pcareduce(
  sce,
  cluster_sizes = 2:4,
  method = "M"
)

sce <- run_pcareduce(
  sce,
  cluster_sizes = 2:4,
  method = "S"
)

sce <- run_raceid(
  sce,
  cluster_sizes = 2:4
)

sce <- run_reduction_clustering(
  sce,
  reductions = c("PCA", "TSNE"),
  methods = c(
    "kmeans",
    "hierarchical",
    "louvain",
    "leiden"
  ),
  cluster_sizes = 2:4,
  knn = c(5, 10),
  resolution = c(0.5, 1)
)

sce <- run_validation(sce)

sce <- scale_validation(sce)

sce <- run_rank_aggregation(sce)

S4Vectors::metadata(sce)$ranking
```

## Methodology

The overall framework implemented in `sclValid` follows the methodology introduced by Visser and Datta (2025). Multiple clustering procedures are evaluated using a diverse set of validation measures representing different aspects of clustering performance. These measures are standardized to a common direction and then combined using weighted rank aggregation to identify clustering solutions that perform consistently across validation criteria.

Several stability-based validation measures in the package are adaptations of earlier clustering validation methods. Full methodological references are provided in the documentation for the corresponding functions.

## Citation

If you use `sclValid` in your research, please cite:

Visser O, Datta S. Integrating Multiple Clustering Techniques and Performance Measures via Ranking for scRNA-Seq Data. *Statistics in Medicine*. 2025;44(28-30):e70331. https://doi.org/10.1002/sim.70331

A package-specific citation entry may also be added in a future release.

## References

The package includes implementations or adaptations of established clustering validation measures, including methods based on work by Dunn, Rousseeuw, Hubert and Arabie, Datta and Datta, Kapp and Tibshirani, and Pihur and colleagues.

The clustering procedures implemented or supported in the package include pcaReduce, RaceID, K-means, hierarchical clustering, Louvain clustering, and Leiden clustering.

See the individual function documentation for complete references.

## License

License information is provided in the package `DESCRIPTION` file.

## Author

Owen Visser
