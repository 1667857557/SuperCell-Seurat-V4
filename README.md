# SuperCell

SuperCell builds metacells from single-cell RNA, CITE-seq, and multiome Seurat objects. The current code defaults to Seurat v4-style assay outputs while remaining compatible with Seurat v5 objects.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("1667857557/SuperCell_Seurat_V4")
library(SuperCell)
```

Key runtime dependencies include Seurat, SeuratObject, Signac, data.table, fs, foreach, doParallel, future, future.apply, and pbapply.

## Minimal Seurat workflow

### 1. Preprocess the single-cell object

Run the usual Seurat preprocessing first. `SCimplify_for_Seurat()` expects reductions such as `pca`, `lsi`, `apca`, or other user-provided embeddings to already exist.

```r
DefaultAssay(obj) <- "RNA"
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
```

### 2. Build metacells

Unimodal RNA metacells:

```r
mc <- SCimplify_for_Seurat(
  seurat = obj,
  assay = "RNA",
  reduction = list("pca"),
  dims = list(1:30),
  gamma = 20
)
validate_metacell_output(mc)
```

Multimodal RNA + ADT metacells:

```r
mc <- SCimplify_for_Seurat(
  seurat = obj,
  assay = c("RNA", "ADT"),
  reduction = list("pca", "apca"),
  dims = list(1:30, 1:18),
  gamma = 20
)
validate_metacell_output(mc)
```

Optional label-aware construction uses any categorical metadata column. For example, a cell-type label prevents metacells from mixing known cell types; `NA` labels are allowed for partial annotation:

```r
mc <- SCimplify_for_Seurat(
  seurat = obj,
  assay = c("RNA", "ADT"),
  reduction = list("pca", "apca"),
  dims = list(1:30, 1:18),
  label = "celltype",
  gamma = 20
)
```

### Independent cell-type WNN graphs with conditions pooled within cell type

Use `SCimplify_by_graph_group()` for paired RNA+ATAC data when every broad cell type must have an independent multimodal graph while all conditions within that cell type share one graph geometry.

```r
sc <- SCimplify_by_graph_group(
  seurat = obj,
  cell.graph.group = obj$cell_type,
  cell.split.condition = obj$condition,
  assay = c("RNA", "ATAC"),
  reduction = list("pca", "lsi"),
  dims = list(1:30, 2:30),
  gamma = 20,
  k.knn = 30
)
```

The contract is:

- one independent multimodal WNN graph is built for each `cell.graph.group` value;
- RNA and ATAC modality weights are learned adaptively by the native SuperCell/Seurat WNN implementation within that graph group;
- all conditions in a graph group jointly participate in neighbour search and Walktrap clustering;
- `cell.split.condition` is applied only after clustering, so final metacells are condition-pure without constructing condition-specific graphs;
- returned metacell IDs are globally unique and the cell-level table records both parent WNN membership and final condition-pure membership.

Do not split the input by condition before calling this function. RNA PCA and ATAC LSI reductions must be defined on a shared coordinate system, not recomputed independently for each condition.

`sample_col` is not a SuperCell builder argument. The relevant interfaces are:

- `SCimplify_for_Seurat()` builds one graph on the supplied Seurat object; `label` creates label-restricted graph components.
- `SCimplify_by_graph_group()` builds one native multimodal WNN graph per graph group and applies condition only after clustering.
- `SCimplify()` and `SCimplify_from_embedding()` remain matrix/embedding builders for workflows that do not require the grouped Seurat WNN contract.

The returned Seurat object from `SCimplify_for_Seurat()` contains aggregated assays, metacell size in `mc$size`, categorical metadata assignments, purity columns, and run metadata in `mc@misc`. Run `validate_metacell_output(mc)` after construction to check assay colnames, metacell size metadata, optional membership tables, and optional fragment manifests.

## Metacell expression aggregation

Use `MetacellExpression()` when metacell IDs already exist in metadata. Output column names preserve the grouping IDs by default.

```r
obj$metacell_id <- c("MC_1", "MC_1", "MC_2")
mat <- MetacellExpression(
  object = obj,
  assays = "RNA",
  group.by = "metacell_id",
  slot = "counts",
  return.seurat = FALSE
)[["RNA"]]
```

Custom output names can be supplied with `metacell.names`.

## Fragment aggregation for ATAC/multiome

`AggregateFragmentFile()` expects a named vector mapping single-cell barcodes to final metacell barcodes. Unmatched barcodes are removed by default.

```r
membership <- c(
  "AAAC-1" = "Metacell_1",
  "AAAG-1" = "Metacell_2"
)

mc_fragments <- AggregateFragmentFile(
  input_file = "fragments.tsv.gz",
  membership = membership,
  output_path = "metacell_fragments",
  nb_cl = 4
)
```

For `SCimplify_for_Seurat()`, pass fragment files by chromatin assay. Multiple files per assay are supported.

```r
mc <- SCimplify_for_Seurat(
  seurat = obj,
  assay = c("RNA", "ATAC"),
  reduction = list("pca", "lsi"),
  dims = list(1:30, 2:30),
  fragmentFiles = list(ATAC = c("sample1/fragments.tsv.gz", "sample2/fragments.tsv.gz")),
  gamma = 20
)
```

`bgzip` and `tabix` must be available in `PATH` or passed with `bgzip_path` and `tabix_path`. For fragment-producing workflows, validate with `validate_metacell_output(mc, require_fragments = TRUE)`.

## Useful plotting helpers

```r
DimPlotSC(obj, mc, reduction = "umap", metacell.col = "celltype")
DimPlot.SuperCell(mc, reduction = "umap", group.by = "celltype")
VlnPlot.SuperCell(mc, features = c("CD3D", "MS4A1"), group.by = "celltype")
FeatureScatter.SuperCell(mc, feature1 = "rna_CD14", feature2 = "adt_CD14")
```

## Conversion from legacy SuperCell objects

```r
seurat_mc <- supercell_2_Seurat(
  SC.GE = SC.GE,
  SC = SC,
  fields = c("ident"),
  output.assay.version = "v4"
)
```

## Tutorials

Long-form rendered tutorials are available under `docs/tutorials/`. The package vignette `vignettes/a_SuperCell.Rmd` contains a short runnable example.

## Citation

If you use SuperCell, please cite the SuperCell publications listed in the manuscript and package documentation.
