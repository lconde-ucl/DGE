#!/usr/bin/env nextflow

/*
========================================================================================
                         nf-core/DGE
========================================================================================
 nf-core/DGE Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/DGE
----------------------------------------------------------------------------------------
*/



def helpMessage() {
    log.info """
    =======================================================
                                             ,--./,-.
              ___     __   __   __   ___    /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__        }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     nf-core/DGE v${workflow.manifest.version}
    =======================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow_dge --inputdir results_rnaseq --metadata metadata.txt --outdir results_DGE 

    Mandatory arguments:
      --inputdir                    Path to nextflow_rnaseq results folder [results]
      --metadata                    Path to metadata. This should be a txt file where the first column are the sample IDs, 
                                        and the other (1 or more) columns displays the conditions for each sample. The samples 
                                        must match those in the featureCounts matrix data located in inputdir [metadata.txt] 
    Options - kallisto mode:
      --kallisto                    Run DESEq2 on kallisto abundance files instead of on featureCounts matrix. Requires specifying the assembly [null]
      --assembly                    Required when in kallisto mode, should be the same assembly used when running kallisto. Possible values are hg19, hg38, or mm10 [null]

    Options - deseq2 model:
      --design                      Specifies DESeq2 design. If defined, --condition, --treatment and --control must also be defined [null]
      --condition                   Specifies 'condition' for the DESeq2 contrast. Requires --design to be specified [null]
      --treatment                   Specifies 'treatment' for the DESeq2 contrast. Requires --design to be specified [null]
      --control                     Specifies 'control' for the DESeq2 contrast. Requires --design to be specified [null]

    Options - gsea (human only):
      --skip_gsea                   Skip GSEA step, otherwise it will run GSEA on each result file [false]
      --gmt                         File with gene sets in GMX format. If not specified, it will use the hallmark gene sets from MSigDB [null]
      --min_set NUM                 Ignore gene sets that contain less than NUM genes [15]";
      --max_set NUM	            Ignore gene sets that contain more than NUM genes [500]";
      --perm NUM                    Number of permutations [1000]";
      

    Options - other:
      --outdir                      The output directory where the results will be saved [results_DGE]
      --pval                        Pval threshold to display gene labels in the volcano plot [1e-50]
      --fc                          FC threshold to display gene labels in the volcano plot [3]
    """.stripIndent()
}


// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// Configurable variables and validate inputs
params.inputdir="results"
params.metadata="metadata.txt"
params.outdir = "results_DGE"
params.kallisto = false
params.assembly = false
params.design = false
params.condition = false
params.treatment = false
params.control = false
params.skip_gsea = false
params.perm = 1000
params.min_set = 15
params.max_set = 500
params.pval = 1e-50
params.fc = 3

if( !params.inputdir ){
    exit 1, "No inputdir specified!"
}
inputdir = file(params.inputdir)
if( !inputdir.exists() ) exit 1, "Inputdir folder not found: ${params.inputdir}. Specify path with --inputdir."

if( !params.metadata ){
    exit 1, "No metadata specified!"
}
metadata = file(params.metadata)
if( !metadata.exists() ) exit 1, "Metadata file not found: ${params.metadata}. Specify path with --metadata."

if ((params.design) && (!params.condition || !params.treatment || !params.control)){
    exit 1, "Invalid arguments: --design \'${params.design}\' requieres --condition, --treatment and --control. Pelase specify all of them or run the pipeline without specifying any design"
}
if ((params.condition) && (!params.design || !params.treatment || !params.control)){
    exit 1, "Invalid arguments: --condition \'${params.condition}\' requieres --design, --treatment and --control. Pelase specify all of them or run the pipeline without specifying any design"
}
if ((params.treatment) && (!params.design || !params.condition || !params.control)){
    exit 1, "Invalid arguments: --treatment \'${params.treatment}\' requieres --design, --condition and --control. Pelase specify all of them or run the pipeline without specifying any design"
}
if ((params.control) && (!params.design || !params.treatment || !params.condition)){
    exit 1, "Invalid arguments: --control \'${params.control}\' requieres --design, --condition and --treatment. Pelase specify all of them or run the pipeline without specifying any design"
}
if (params.kallisto  && !params.assembly){
    exit 1, "Running the pipeline in kallisto mode requires specifying an assembly. Valid options: 'hg19', 'hg38', 'mm10'"
}
if (params.kallisto && (params.assembly!= 'hg19' && params.assembly != 'hg38' && params.assembly != 'mm10')){
    exit 1, "Invalid assembly option: ${params.assembly}. Valid options: 'hg19', 'hg38', 'mm10'"
}
if(params.kallisto && params.assembly){
	params.tx2gene = params.assembly ? params.genomes[ params.assembly ].tx2gene ?: false : false
	tx2gene = file(params.tx2gene)
	if( !tx2gene.exists() ) exit 1, "tx2gene file not found: ${params.tx2gene}"
}else{
	tx2gene='-'
}
if(!params.skip_gsea){
	if (params.assembly && params.assembly != 'hg19' && params.assembly != 'hg38'){
	    exit 1, "GSEA can only be run on human (assembly hg19 or hg38). You are indicating that your assembly is ${params.assembly}. Please correct the assembly if this was a mistake, or use --skip_gsea to skip the GASE step"
	}else{
		gmx = file(params.gmx)
		if( !gmx.exists() ) exit 1, "GMX file not found: ${params.gmx}"
	}
}else{
	gmx='-'
}


// Header 
println "========================================================"
println "                                       ,--./,-.         "
println "          ___     __  __   __  ___    /,-._.--~\'       "
println "    |\\ | |__  __ /  `/  \\ |__)|__        }  {         "
println "    | \\| |       \\__ \\__/ |  \\|___   \\`-._,-`-,    "
println "                                      `._,._,\'         "
println "                                                        "
println "            D G E    P I P E L I N E                    "
println "========================================================"
println "['Pipeline Name']     = nf-core/DGE"
println "['Pipeline Version']  = $workflow.manifest.version"
println "['Inputdir']          = $params.inputdir"
if(params.design != "-"){
	println "['DESeq2 design']     = $params.design"
	println "['DESeq2 condition']  = $params.condition"
	println "['DESeq2 treatment']  = $params.treatment"
	println "['DESeq2 control']    = $params.control"
}else{
	println "['DESeq2 design']     = No design specified"
}
if(params.kallisto){
	println "['Read counts mode']  = kallisto"
	println "['Assembly']          = $params.assembly"
}else{
	println "['Read counts mode']  = featureCounts"
}
if(!params.skip_gsea){
	println "['GSEA step']         = True"
	println "['GMX set']           = $params.gmx"
	if(!params.assembly){
		println "['Species']           = *Assuming Human*"
	}
}else{
	println "['GSEA step']         = False"
}
if(params.pval != "-"){
	println "['Volcano plot Pval threshold'] = $params.pval"
}
if(params.fc != "-"){
	println "['Volcano plot FC threshold']   = $params.fc"
}
println "['Output dir']        = $params.outdir"
println "['Working dir']       = $workflow.workDir"
println "['Container Engine']  = $workflow.containerEngine"
println "['Current home']      = $HOME"
println "['Current user']      = $USER"
println "['Current path']      = $PWD"
println "['Working dir']       = $workflow.workDir"
println "['Script dir']        = $workflow.projectDir"
println "['Config Profile']    = $workflow.profile"
println "========================================================"



/*
 * STEP 1 - DESeq2
 */

process deseq2 {
    publishDir "${params.outdir}", mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("_MAplot.png") > 0) "deseq2_PlotsAndFiles/MAplots/$filename"
            else if (filename.indexOf("_VolcanoPlot.png") > 0) "deseq2_PlotsAndFiles/volcanoPlots/$filename"
            else if (filename.indexOf("_heatmap.png") > 0) "deseq2_PlotsAndFiles/heatmapPlots/$filename"
            else if (filename.indexOf(".txt") > 0) "deseq2_PlotsAndFiles/files/$filename"
            else if (filename == "PCAplot.png") "deseq2_PlotsAndFiles/$filename"
	    else "$filename"
    }

    output:
    file "*.txt" into stats_for_gsea
    file "*{_MAplot.png,_VolcanoPlot.png,_heatmap.png}"
    file "PCAplot.png"
    file "deseq2_report"

    script:
    """

cat > deseq2.conf << EOF
[deseq2]

INPUTDIR = $inputdir
METADATA = $metadata
DESIGN = $params.design
CONDITION = $params.condition
TREATMENT = $params.treatment
CONTROL = $params.control
PVAL = $params.pval
FC = $params.fc
KALLISTO = $params.kallisto
TX2GENE = $tx2gene

EOF

    run_deseq2.R $inputdir $metadata deseq2.conf
    mv *_AllPlot.png deseq2_report/figuresDESeq2_nextflow_pipeline_results/.
    cp PCAplot.png deseq2_report/figuresDESeq2_nextflow_pipeline_results/.
    """
}


/*
 * STEP 2 - GSEA
 */

process gsea {
    tag "$contrast"
    publishDir "${params.outdir}/gsea_results/$contrast", mode: 'copy'

    when:
    !params.skip_gsea

    input:
    file stats from stats_for_gsea.flatten()

    output:
    file "gsea.*"

    script:        
    contrast = stats.toString() - ~/_results.txt$/
    """
    echo $stats
    get_metric.pl $stats $contrast
    run_gsea.pl --rnk ${contrast}.rnk --gmx $gmx --perm $params.perm --min_set $params.min_set --max_set $params.max_set
    mv gsea_results*/* .
    plot_gsea.R gsea_table.txt $params.perm
    """
}


workflow.onComplete { 
    println ( workflow.success ? "Done! Wrapping up..." : "Oops .. something went wrong" )
    log.info "[nf-core/test] Pipeline Complete"

}


