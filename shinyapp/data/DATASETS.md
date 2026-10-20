# Demonstration datasets

The app ships compact, fully runnable demonstration inputs:

- `beatAML2_demo.csv`: bulk/pseudobulk workflow demonstration.
- `bulk_expression_demo.csv`: gene-by-sample bulk matrix sent through the bundled MuSiC basis.
- `Van_scRNA_demo.csv`: annotated scRNA-seq aggregation workflow demonstration.
- `Visium_spot_demo.csv`: annotated Visium spot workflow demonstration (144 spots on a small tissue-style grid).

They are intentionally compact so the app can be tested locally and deployed without a large download. Replace them with the corresponding public matrices or Seurat/Visium `meta.data` exports after mapping annotations to the nine frozen RNHI lineages.
