# DGE: Usage

## Table of contents

* [Running the pipeline](#running-the-pipeline)
* [Main arguments](#main-arguments)
    * [`-profile`](#-profile-single-dash)
        * [`legion`](#legion)
        * [`myriad`](#myriad)
    * [`--inputdir`](#--inputdir)
    * [`--metadata`](#--metadata)
* [Arguments - kallisto mode](#kallisto-arguments)
    * [`--kallisto`](#--kallisto)
    * [`--assembly`](#--assembly)
* [Arguments - deseq2 model](#deseq2-arguments)
    * [`--design`](#--design)
    * [`--condition`](#--condition)
    * [`--treatment`](#--treatment)
    * [`--control`](#--control)
* [Arguments - gsea (human only)](#gsea-arguments)
    * [`--skip_gsea`](#--skip_gsea)
    * [`--gmx`](#--gmx)
    * [`--gmx_ensembl`](#--gmx_ensembl)
    * [`--min_set`](#--min_set)
    * [`--max_set`](#--max_set)
    * [`--perm`](#--perm)
* [Job Resources](#job-resources)
* [Automatic resubmission](#automatic-resubmission)
* [Other command line parameters](#other-command-line-parameters)
    * [`--outdir`](#--outdir)
    * [`--pval`](#--pval)
    * [`--fc`](#--fc)
    * [`-resume`](#-resume-single-dash)
    * [`--max_memory`](#--max_memory)
    * [`--max_time`](#--max_time)
    * [`--max_cpus`](#--max_cpus)


## Running the pipeline
This is a basic pipeline for differential gene expression analysis that is meant to be run after the data has been processed with the 
[nextflow_rnaseq](https://github.com/UCL-BLIC/rnaseq) pipeline (or the [nfcore rnaseq](https://github.com/nf-core/rnaseq) pipeline, tested on v1.4.2) 
and therefore a featureCounts gene counts file (or kallisto abundance files) have been generated.
 
The typical command for running the pipeline  with the "nextflow_mergefastq' alias is as follows:
```bash
module load blic-modules
module load nextflow_dge

nextflow_dge --inputdir results_rnaseq --metadata metadata.txt --outdir results_DGE
```

This will launch the pipeline with the `legion` or `myriad` configuration profile, depending on where you submit the job from.

Note that the pipeline will create the following files in your working directory:

```bash
work                     # Directory containing the nextflow working files
merged_fastq_files       # Finished results (configurable, see below)
.nextflow_log            # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Main Arguments

To see all the available arguments, use the `--help` flag
```bash
nextflow_dge --help
```

The main arguments are:

### `-profile`
This parameter is NOT necessary as the shortcut `nextflow_merge_fastq` takes care of selecting the appropiate configuration profile. But just for your information, profiles are used to give 
configuration presets for different compute environments.

* `legion`
    * A generic configuration profile to be used with the UCL cluster legion
* `myriad`
    * A generic configuration profile to be used with the UCL cluster myriad

### `--inputdir`
This is used to specify the location of the results folder obtained after running the [nextflow_rnaseq](https://github.com/UCL-BLIC/rnaseq) pipeline (or the [nfcore rnaseq](https://github.com/nf-core/rnaseq) pipeline, 
tested on v1.4.2). For example:

```bash
--inputdir 'path/to/data/results_from_nextflow_rnaseq_pipeline/'
```

If left unspecified, it will look for the default dir: './results'

The DGE pipeline assumes that the `--inputdir` folder contains either a `[inputdir]/featureCounts/merged_gene_counts.txt` file (if run on normal mode), or one or more `[inputdir]/kallisto/SAMPLENAME/abundance.h5` abundance files (if 
run on kallisto mode). 

Please note that running the nextflow_rnaseq pipeline is not mandatory, as long as you have a featureCounts file or kallisto abundance.h5 files, you can run the DGE pipeline, just organize the files in a folder structure like the 
above. The featureCounts file shuld be formatted as outputted by the [nextflow_rnaseq](https://github.com/UCL-BLIC/rnaseq) pipeline ("ENSEMBL_ID" column with gene names, followed by sample counts columns), or the
[nfcore rnaseq](https://github.com/nf-core/rnaseq) pipeline ("Geneid" and "gene_name" columns with gene names, followed by sample counts columns)

### `--metadata`
This should be a txt file where the first column are the sample IDs, and the other (1 or more) columns displays the conditions for each sample. The samples must match those in the featureCounts matrix data located in inputdir.

Format:
```
SampleID	Status	Levels
sample_1	ctr	high
sample_2	ctr	high
sample_3	ctr	med
sample_4	case	low
sample_5	case	low
sample_6	case	low
```

If left unspecified, it will look for the default dir: './metadata.txt'

## Arguments - kallisto mode

### `--kallisto`
Run DESEq2 on kallisto abundance files instead of on featureCounts matrix. Requires specifying the assembly.
Not used by default.

### `--assembly`
Required when in kallisto mode, should be the same assembly used when running kallisto. Possible values are hg19, hg38, or mm10.
Not used by default.


## Arguments - deseq2 model

By default, the DGE pipeline will run differential gene expression analysis on each possible combination of conditions using a design with all the conditions. For example, for the `metadata.txt` file above, the pipeline will run the 
following analysis:

```
Design: ~ Status + Levels
Comparisons:
	cases vs. controls (status)
	high vs. medium (levels)
	high vs. low (levels)
	medium vs. low (levels)
```

This default behaviour (all possible comparisons) can be overrided and the user can choose the design and comparison of interest by specifying the following arguments:
 
### `--design`
Specifies DESeq2 design. If defined, --condition, --treatment and --control must also be defined.
Not used by default.

### `--condition`
Specifies 'condition' for the DESeq2 contrast. Requires --design to be specified.
Not used by default.

### `--treatment`
Specifies 'treatment' for the DESeq2 contrast. Requires --design to be specified.
Not used by default.

### `--control`
Specifies 'control' for the DESeq2 contrast. Requires --design to be specified.
Not used by default.


## Arguments - GSEA model

For each comparison above, a GSEA analysis using the hallmark gene sets from MSigDB will be performed. Please note that the hallmark dataset contains HUGO IDs. If your gene counts contain
Ensembl IDs (it will depend on what GFT file you used in the featurecounts step), you need to add the --gmx_ensembl flag. Also, if your data is from a species other than human, the default hallmark 
gene set will not work for your data, and you will have to either skip GSEA with the --skip_gsea flag, or add an appropiate gene set with the --gmx argument.

 
### `--skip_gsea`
Skip GSEA step, otherwise it will run GSEA on each result file.
Not used by default.

### `--gmx`
File with gene sets in GMX format. If not specified, it will use the hallmark gene sets from MSigDB (human HUGO IDs).

### `--gmx_ensembl`
Use a version of the MSigDB hallmark gene set with Ensembl IDs, obtained using the msigdbr R package. This flag overriddes the --gmx argument. 

### `--min_set NUM`
Ignore gene sets that contain less than NUM genes.
Default = 15

### `--max_set NUM`
Ignore gene sets that contain more than NUM genes.
Default = 500

### `--perm NUM`
Number of permutations for the NES calculation.
Default = 1000


## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

## Other command line parameters

### `--pval`
Pval threshold to display gene labels in the output volcano plot.
Default:  1e-50

### `--fc`
FC threshold to display gene labels in the output volcano plot.
Default: 3

### `--outdir`
The output directory where the results will be saved.
Default: results_DGE

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
Please note that since this pipeline only runs one process, the -resume option is not useful here. This might change if more processes are added to the pipeline in the future.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'``

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each process.
Should be a string in the format integer-unit. eg. `--max_cpus 1`

