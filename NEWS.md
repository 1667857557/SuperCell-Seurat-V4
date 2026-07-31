# SuperCell 2.1.0

- Adds `SCimplify_by_graph_group()` as the single canonical grouped multimodal metacell builder.
- Builds one native RNA+ATAC WNN graph per graph group while pooling all conditions during neighbour construction and Walktrap clustering.
- Applies condition only after clustering to produce condition-pure final memberships on a shared graph geometry.
- Removes the obsolete `SCimplify_by_graph_group_from_embedding()` API, export, help page and embedding-specific tests without retaining a compatibility wrapper.
- Unifies the default graining level at `gamma = 30` for `SCimplify_for_Seurat()` and `SCimplify_by_graph_group()`.
- Adds installed-package CI for grouped-WNN API, membership purity, default parameters and obsolete-API removal.
