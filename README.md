# SuperCell

SuperCell builds metacells from single-cell Seurat objects, including paired RNA+ATAC data.

## Installation

```r
remotes::install_github("1667857557/SuperCell_Seurat_V4")
library(SuperCell)
```

## Standard Seurat builder

Precompute the requested reductions, then run:

```r
mc <- SCimplify_for_Seurat(
  seurat = obj,
  assay = "RNA",
  reduction = list("pca"),
  dims = list(1:30),
  gamma = 30
)
```

For multimodal input:

```r
mc <- SCimplify_for_Seurat(
  seurat = obj,
  assay = c("RNA", "ATAC"),
  reduction = list("pca", "lsi"),
  dims = list(1:30, 2:30),
  gamma = 30
)
```

`label` can restrict graph construction by an existing categorical metadata column. `sample_col` and `condition_col` are not formals of `SCimplify_for_Seurat()`.

`return.graph = TRUE` is an internal/advanced compact-return option: it requires `return.seurat = FALSE` and returns the exact graph used for Walktrap. The grouped builder uses it transiently for repair and discards the graph afterward.

## Shared-WNN grouped builder

Use `SCimplify_by_graph_group()` when each broad cell type requires its own RNA+ATAC WNN while conditions must share that graph geometry:

```r
sc <- SCimplify_by_graph_group(
  seurat = obj,
  cell.graph.group = obj$cell_type,
  cell.split.condition = obj$condition,
  assay = c("RNA", "ATAC"),
  reduction = list("pca", "lsi"),
  dims = list(1:30, 2:30),
  gamma = 30,
  k.knn = 30
)
```

The grouped contract is:

```text
one WNN per cell.graph.group
  -> all conditions jointly in WNN + Walktrap
  -> split parent membership by condition
  -> optional small-metacell repair on the same original WNN
  -> final condition-pure membership
```

To enforce a minimum final size, provide an explicit affinity threshold:

```r
sc <- SCimplify_by_graph_group(
  seurat = obj,
  cell.graph.group = obj$cell_type,
  cell.split.condition = obj$condition,
  assay = c("RNA", "ATAC"),
  reduction = list("pca", "lsi"),
  dims = list(1:30, 2:30),
  min_metacell_size = 10L,
  min_metacells_per_stratum = 3L,
  min_merge_affinity = 0.05,
  unresolved_small_policy = "error"
)
```

Repair candidates are restricted to the same condition and graph group. `membership_table` contains final `metacell_id` and original `parent_metacell_id` provenance.

## Existing memberships

Use `MetacellExpression()` to aggregate assays from an existing grouping column:

```r
mat <- MetacellExpression(
  object = obj,
  assays = "RNA",
  group.by = "metacell_id",
  slot = "counts",
  return.seurat = FALSE
)[["RNA"]]
```

## ATAC fragments

`AggregateFragmentFile()` accepts a named single-cell-to-metacell membership vector:

```r
mc_fragments <- AggregateFragmentFile(
  input_file = "fragments.tsv.gz",
  membership = membership,
  output_path = "metacell_fragments",
  nb_cl = 4
)
```

`SCimplify_for_Seurat()` also accepts fragment files through `fragmentFiles` for chromatin assays. Use `validate_metacell_output()` to validate generated metacell objects and fragment manifests.

## Documentation

Use the package Rd pages for the complete current function signatures and arguments.
