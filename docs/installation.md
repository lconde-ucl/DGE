# DGE: Installation

## Installation at UCL

The DGE pipeline is already installed in myriad. You just need to load it and start using it:

```bash
module load blic-modules
module load nextflow_dge

nextflow_dge --inputdir results_rnaseq --metadata metadata.txt --outdir results_DGE 
```

By default, the pipeline runs with the `legion` configuration profile [`conf/legion.config`](../conf/legion.config) if you submit it from legion, and with the `myriad` config 
[`conf/myriad.config`](../conf/myriad.config) if you send the job from myriad.

The 'standard' configuration (using the `local` executor) is not enabled.


## Installation at other sites

You'll need to install 
* nextflow, eg via bioconda or https://www.nextflow.io/index.html
* the R packages mentioned in bin/run_deseq2.R

