#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/deseq2
========================================================================================
 nf-core/deseq2 Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/nf-core/deseq2
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
     nf-core/deseq2 v${workflow.manifest.version}
    =======================================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/deseq2 --inputdir results_rnaseq --metadata metadata.txt --outdir results_deseq2 

    Mandatory arguments:
      --inputdir                    Path to nextflow_rnaseq results folder [results]
      --metadata                    Path to metadata. This should be a txt file where the first column are the sample IDs, 
                                        and the other (1 or more) columns displays the conditions for each sample. The samples 
                                        must match those in the featureCounts matrix data located in inputdir [metadata.txt] 
    Options:
      --design                      Specifies DESeq2 design. If defined, --condition, --treatment and --control must also be defined [null]
      --condition                   Specifies 'condition' for the DESeq2 contrast. Requires --design to be specified [null]
      --treatment                   Specifies 'treatment' for the DESeq2 contrast. Requires --design to be specified [null]
      --control                     Specifies 'control' for the DESeq2 contrast. Requires --design to be specified [null]
      --outdir                      The output directory where the results will be saved [results_deseq2]
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
params.outdir = "results_deseq2"
params.design = "-"
params.condition = "-"
params.treatment = "-"
params.control = "-"
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

if ((params.design != "-") && (params.condition == "-" || params.treatment == "-" || params.control == "-")){
    exit 1, "Invalid arguments: --design \'${params.design}\' requieres --condition, --treatment and --control. Pelase specify all of them or run the pipeline without specifying any design"
}
if ((params.condition != "-") && (params.design == "-" || params.treatment == "-" || params.control == "-")){
    exit 1, "Invalid arguments: --condition \'${params.condition}\' requieres --design, --treatment and --control. Pelase specify all of them or run the pipeline without specifying any design"
}
if ((params.treatment != "-") && (params.design == "-" || params.condition == "-" || params.control == "-")){
    exit 1, "Invalid arguments: --treatment \'${params.treatment}\' requieres --design, --condition and --control. Pelase specify all of them or run the pipeline without specifying any design"
}
if ((params.control) && (!params.design || !params.treatment || !params.condition)){
    exit 1, "Invalid arguments: --control \'${params.control}\' requieres --design, --condition and --treatment. Pelase specify all of them or run the pipeline without specifying any design"
}



// Header 
println "========================================================"
println "                                       ,--./,-.         "
println "          ___     __  __   __  ___    /,-._.--~\'       "
println "    |\\ | |__  __ /  `/  \\ |__)|__        }  {         "
println "    | \\| |       \\__ \\__/ |  \\|___   \\`-._,-`-,    "
println "                                      `._,._,\'         "
println "                                                        "
println "           D E S E Q 2    P I P E L I N E               "
println "========================================================"
println "['Pipeline Name']     = nf-core/deseq2"
println "['Pipeline Version']  = workflow.manifest.version"
println "['Inputdir']          = $params.inputdir"
if(params.design != "-"){
	println "['DESeq2 design']     = $params.design"
	println "['DESeq2 condition']  = $params.condition"
	println "['DESeq2 treatment']  = $params.treatment"
	println "['DESeq2 control']    = $params.control"
}else{
	println "['DESeq2 design']     = No design specified"
}
if(params.pval != "-"){
	println "['Volcano plot Pval threshold'] = $params.pval"
}
if(params.fc != "-"){
	println "['Volcano plot FC threshold']   = $params.fc"
}
println "['Output dir']        = $params.outdir"
println "['Working dir']       = workflow.workDir"
println "['Container Engine']  = workflow.containerEngine"
println "['Current home']      = $HOME"
println "['Current user']      = $USER"
println "['Current path']      = $PWD"
println "['Working dir']       = workflow.workDir"
println "['Script dir']        = workflow.projectDir"
println "['Config Profile']    = workflow.profile"
println "========================================================"



/*
 * STEP 1 - DESeq2
 */

process deseq2 {
    publishDir "${params.outdir}", mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf("_MAplot.png") > 0) "plots/MAplots/$filename"
            else if (filename.indexOf("_VolcanoPlot.png") > 0) "plots/volcanoPlots/$filename"
            else if (filename.indexOf(".txt") > 0) "files/$filename"
            else if (filename == "PCAplot.png") "plots/$filename"
	    else "$filename"
    }

    output:
    file "*.txt"
    file "*{_MAplot.png,_VolcanoPlot.png}"
    file "PCAplot.png"
    file "report"

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

EOF

    run_deseq2.R $inputdir $metadata deseq2.conf
    mv *_AllPlot.png report/figuresRNAseq_analysis_with_DESeq2/.
    cp PCAplot.png report/figuresRNAseq_analysis_with_DESeq2/.
    """
}


workflow.onComplete { 
    println ( workflow.success ? "Done! Wrapping up..." : "Oops .. something went wrong" )
    log.info "[nf-core/test] Pipeline Complete"

}


