# Coleoptera BUSCO Phylogenomics Pipeline

Author and maintainer: **Hehao Ye** (hehao.ye16@ic.ac.uk) — Imperial College London / Natural History Museum, London.

Code accompanying the manuscript:

> **An updated phylogenomic tree of Coleoptera (beetles).**
> Hehao Ye, Zichen Zhou, Fiona L. Carpenter, Huaxi Liu, Alfried P. Vogler. *Molecular Phylogenetics and Evolution* (in review, 2026).

Pipeline to assemble Coleoptera genomes/transcriptomes, extract BUSCO orthologs, apply locus filtering, and infer species trees under multiple substitution models. The published tree uses **263 Coleoptera + 13 outgroup taxa** with **2,150 loci** at 60% taxon occupancy under the PMSF model (LG+C20+FO+R).

---

## Pipeline overview

```
01_assembly  →  02_busco  →  03_extraction_filtering  →  04_tree_building  →  05_visualization
```

## Repository structure

```
.
├── README.md
├── LICENSE                              MIT
├── environment.yml                      conda environments
│
├── 01_assembly/                         de novo assembly per data type
│   ├── Genome_Assembly_Pair-end_short.sh    Illumina paired-end
│   ├── Genome_Assembly_Single-end_short.sh  Illumina single-end
│   ├── Genome_Assembly_Single_long_HiFi.sh  PacBio HiFi
│   ├── Genome_Assembly_Single_long_CLR.sh   PacBio CLR
│   ├── Genome_Assembly_Single_long_ONT.sh   Oxford Nanopore
│   └── Transcriptome_Assembly.sh            Illumina RNA-seq
│
├── 02_busco/run_busco.sh                BUSCO v6 against coleoptera_odb12
│
├── 03_extraction_filtering/
│   ├── BUSCO_extraction.sh              align, trim, TreeShrink, gappy filter,
│   │                                    8-criteria locus filter, FASconCAT-G
│   └── filter_by_coverage.py            per-sequence coverage filter
│
├── 04_tree_building/tree_building.sh    IQ-TREE (LG, GHOST, PMSF, Dayhoff6),
│                                        ASTRAL-IV, wASTRAL; gCF/sCF, AU test,
│                                        RF/GRF distances, tree statistics
│
└── 05_visualization/
    ├── tree_main_with_support.R         Figure 5 (PMSF 60% retained tree)
    ├── grf_dis_plot_clean.R             compute GRF distances (TreeDist)
    ├── plot_rf_distance.py              Figure 3 (standard RF MDS)
    └── plot_grf_distance.py             Suppl Fig 2 (GRF MDS panels)
```

## Pipeline-to-paper map

| Step | Script(s) | Paper |
|------|-----------|-------|
| Assembly + BUSCO | `01_assembly/*.sh`, `02_busco/run_busco.sh` | Methods 'Specimen selection and sources of genome data'; Fig 1; Suppl Fig 1; Suppl Table 1 |
| Locus filtering & matrix concatenation | `03_extraction_filtering/*` | Methods 'Ortholog Extraction and Filtering' (filtering) + 'Phylogenomic Inference' (concatenation); Fig 2; Table 1 |
| Tree inference | `04_tree_building/tree_building.sh` | Methods 'Phylogenomic Inference'; Fig 5; Suppl Figs 4–6 |
| gCF / sCF | `tree_building.sh` | Fig 5 (support values); Table 3 |
| AU / SH / KH | `tree_building.sh` | Table 2; Suppl Table 4 |
| RF / GRF + MDS | `tree_building.sh`, `plot_rf_distance.py`, `grf_dis_plot_clean.R` + `plot_grf_distance.py` | Fig 3; Suppl Fig 2 |
| Tree statistics (ABS, DVMC, treeness) | `tree_building.sh` | Table 3; Suppl Table 5 |

Not produced by code in this repository (manual / external):
- Fig 4 (per-tree runtime / memory plot)
- Suppl Fig 3 (collapsed-gCF<15 tree derived from `PMSF_Retained_60`)
- Fig 6 (topological comparison with previously published Coleoptera studies)

## Data

- **Raw reads & reference assemblies**: NCBI BioProject `PRJNA1462652`. Provenance summary in Suppl Table 2.
- **Filtered alignments, gene trees, supermatrices, species trees, taxonomy CSV**: Zenodo `doi.org/10.5281/zenodo.20111661`.

The pipeline expects `${PROJECT_ROOT}/Coleoptera.csv` (master taxonomy) — provided in the Zenodo deposit.

## Software

See `environment.yml` for a ready-to-use conda specification. Auxiliary environments (BESST, redundans, medaka) are documented inline because of conflicting Python pins.

Some tools require manual install: FASconCAT-G, ASTRAL v5.7.8 (jar), wASTRAL/ASTER, Phylogears (`pgrecodeseq.pl`), HiFiAdapterFilt.

## Usage

Scripts are documented walkthroughs intended to be **run section by section** (each `# n.` block) after setting `PROJECT_ROOT` and editing input lists / resource flags. Per-step memory and runtime are commented inline (e.g. PMSF 60% required ~228 GB and ~6 days on a single 32-core node).

Quickstart for the visualisation step (the four scripts in `05_visualization/` can be run end-to-end once trees and distance matrices have been produced):

```bash
export PROJECT_ROOT=/path/to/your/project

# Standard RF distance plot (Figure 3)
python 05_visualization/plot_rf_distance.py \
    --rfdist rf_all.rfdist --map all_trees.map \
    --out rf_mds.png --coords rf_mds_coords.tsv

# Generalised RF: compute in R, plot in Python (Suppl Fig 2)
Rscript 05_visualization/grf_dis_plot_clean.R \
    --tree all_trees.tree --coords grf_mds_coords.tsv
python 05_visualization/plot_grf_distance.py \
    --coords_in grf_mds_coords.tsv --map all_trees.map --out grf_mds.png

# Main PMSF 60% retained tree (Figure 5)
Rscript 05_visualization/tree_main_with_support.R
```

## Citation

```bibtex
@article{Ye2026Coleoptera,
  title   = {An updated phylogenomic tree of Coleoptera (beetles)},
  author  = {Ye, Hehao and Zhou, Zichen and Carpenter, Fiona L. and Liu, Huaxi and Vogler, Alfried P.},
  journal = {Molecular Phylogenetics and Evolution},
  year    = {2026},
  doi     = {TO BE ADDED}
}
```

## License

MIT — see [LICENSE](LICENSE).
