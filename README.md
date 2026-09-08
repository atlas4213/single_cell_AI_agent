# Containerized Single-Cell RNA-seq Pipeline

An R-based single-cell RNA-seq analysis pipeline built with Seurat.

## Pipeline stages

1. Load 10X count matrices
2. Remove ambient RNA with SoupX
3. Calculate quality-control metrics
4. Filter low-quality cells
5. Detect doublets with DoubletFinder
6. Normalize expression
7. Select highly variable genes
8. Scale expression
9. Run PCA
10. Perform clustering and visualization

## Main tools

- R
- Seurat
- SoupX
- DoubletFinder
- DropletUtils
- Docker
- renv

## Build the Docker image

```bash
docker build -t single-cell-pipeline:1.0 .
