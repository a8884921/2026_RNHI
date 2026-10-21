# Data files

| File | Content |
|---|---|
| `bulk_compendium_composition.csv` | Nine-lineage compositions (wide format) of the 45 healthy reference donors of the BoneMarrowMap atlas (`healthy_reference = 1`; lineage cell counts, at most 150 cells per lineage and donor) and of the AML transcriptomes of the bulk compendium (3,377 transcriptomes, 11 public datasets) in which at least one lineage is detected (n = 3,310; `healthy_reference = 0`; deconvolved proportions). The 67 transcriptomes without any detected lineage cannot be scored and do not affect the detection limits. The `cohort` column gives the source dataset; some datasets are different releases of the same cohort (OHSU2022 and BeatAML2; GSE106291 and GSE146173; GSE12417 and GSE37642). Provides the default healthy center and the fixed detection limits for bulk data. |
| `joint_reference_basis.rds` | Reference for MuSiC deconvolution: signature matrix `D`, cross-subject variance `Sigma` and average library sizes `S` for 12 normal hematopoietic and 7 AML cell states (7,000 genes). |
| `simulated_composition_demo.csv` | Simulated nine-lineage compositions (8 samples). |
| `simulated_expression_demo.csv` | Simulated gene-by-sample expression matrix (1,200 genes, 8 samples). |
| `simulated_scrnaseq_demo.csv` | Simulated scRNA-seq lineage counts in long format (6 samples). |
| `simulated_visium_demo.csv` | Simulated Visium spot compositions with grid coordinates (144 spots). |

Simulated files demonstrate the input formats and do not represent patient data.
