# Deprecated Pipeline Notice

**Important:** The DGE pipeline is now deprecated and has been replaced by a DSL2 version, available at [lconde-ucl/DGE2](https://github.com/lconde-ucl/DGE2)

<br>
<br>


# DGE
**Run DESeq2/GSEA on the featureCounts matrix or kallisto counts obtained from the nextflow_rnaseq pipeline**

[![Build Status](https://travis-ci.org/nf-core/deseq2.svg?branch=master)](https://travis-ci.org/nf-core/deseq2)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)


### Introduction
The pipeline is built using [Nextflow](https://www.nextflow.io). It is a basic pipeline for differential gene expression analysis that is meant to be run after the data has been preprocessed with the 
[rnaseq](https://github.com/UCL-BLIC/rnaseq) pipeline


### Documentation
The nf-core/SGE pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)
4. [Troubleshooting](docs/troubleshooting.md)
