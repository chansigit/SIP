#!/usr/bin/env nextflow
/*
======================================================================
======================================================================
                Single-cell Interchangeable Pipeline
----------------------------------------------------------------------
  ### SIP. Started February 26, 2018
      SIP is a pipeline dedicated to single-cell data processing,
      whose components are interchangeable.

      "For a machine to run smoothly and predictably, 
      its parts must be standard and hence replaceable."
                                        –Charles Eisenstein
  ### Authors:
        Sijie Chen
        Zheng Wei
        Hao Sun
======================================================================
======================================================================
*/

// Usage
// nextflow run /data/hca/chensijie/code/nextflow_test/sip.groovy --readqc=fastp --use_email=false

/*
****************************************
     Groovy colorful print utilities
        +~~~~~~+~~~~~~+~~~~~~~~~~~+
        |  fg  |  bg  |  color    |
        +~~~~~~+~~~~~~+~~~~~~~~~~~+
        |  30  |  40  |  black    |
        |  31  |  41  |  red      |
        |  32  |  42  |  green    |
        |  33  |  43  |  yellow   |
        |  34  |  44  |  blue     |
        |  35  |  45  |  magenta  |
        |  36  |  46  |  cyan     |
        |  37  |  47  |  white    |
        |  39  |  49  |  default  |
        +~~~~~~+~~~~~~+~~~~~~~~~~~+
****************************************
*/

def cprint(fg,bg,content){
    fgcolor=["black"  :30, "red"    :31, "green"  :32,
             "yellow" :33, "blue"   :34, "magenta":35,
             "cyan"   :36, "white"  :37, "default":39 ]

    bgcolor=["black"  :40, "red"    :41, "green"  :42,
             "yellow" :43, "blue"   :44, "magenta":45,
             "cyan"   :46, "white"  :47, "default":49 ]
    (fg, bg) = [ fgcolor[fg], bgcolor[bg] ]
     style = "${(char)27}[$fg;$bg" + "m"
     print style+content
     print "${(char)27}[39;49" + "m"
}

def cprintln(fg,bg,content){
    fgcolor=["black"  :30, "red"    :31, "green"  :32,
             "yellow" :33, "blue"   :34, "magenta":35,
             "cyan"   :36, "white"  :37, "default":39 ]

    bgcolor=["black"  :40, "red"    :41, "green"  :42,
             "yellow" :43, "blue"   :44, "magenta":45,
             "cyan"   :46, "white"  :47, "default":49 ]
    (fg, bg) = [ fgcolor[fg], bgcolor[bg] ]
     style = "${(char)27}[$fg;$bg" + "m"
     print style+content
     println "${(char)27}[39;49" + "m"
}






/*
****************************************
    Command-line parameters
****************************************
*/
// Inputs and outputs
params.cell_home   = "/data/hca/chensijie/SMART1/SmallFASTQfiles"
params.output      = "/data/hca/chensijie/SMART1_products"

// analysis step configs
params.skip_readqc = false
params.readqc      = "fastp"
params.quant       = "salmon"

// Reference
params.salmon_index="/data/hca/chensijie/Reference/transcriptome/SALMON21_Homo_sapiens.GRCh38.cdna.all"

// notification options
params.use_email   = false
params.email       = "chensj16@mails.tsinghua.edu.cn"




/*
****************************************
    Command-line parameters
****************************************
*/
if (params.skip_readqc){
    cprintln("magenta","yellow","Warning: Read QC skipped!")
}

Cells = Channel.fromPath("${params.cell_home}/**/*.fastq")

process ReadQC{
	executor "local"
    when:
        !params.skip_readqc

    input:
        file    ReadFile   from   Cells

    output:
        //stdout  CleanReads
        file    "${ReadFile}_clean.fastq"       into CleanReadFile
        //file    "${ReadFile}_readqc.{json,html}" into ReadQCReport    
        file    "${ReadFile}_report*" into ReadQCReport

    script:
        if (params.readqc=="fastp"){
            """
            # printf "%s" ${ReadFile.toString().trim()}
            fastp --thread 8 \
                  --json ${ReadFile}_report_fastp.json \
                  --html ${ReadFile}_report_fastp.html \
                  --in1  ${ReadFile} \
                  --out1 ${ReadFile}_clean.fastq
            """
        }
        else if (params.readqc=="fastqc"){
            """
            echo "Under construction!"
            """
        }
        else{
            cprintln("yellow","red","Opoos! Unsupported Read-QC tool.")
            error "Invalid Parameter: --readqc=$params.readqc"
        }

}




/*
CleanReadFile.subscribe {
    if (!params.skip_readqc){
        println it.toString()
    }
}
*/


process AlignmentFreeQuantification{
	executor "local"
    input:
        file ReadFile from CleanReadFile
    output:
        file "${ReadFile}_quant" into CellQuantReport

    script:
        """
        salmon quant --index ${params.salmon_index} \
                     --libType A\
                     --unmatedReads ${ReadFile} \
                     --output ${ReadFile}_quant
        """

}



process ReadQCReport{
    publishDir "${params.output}", mode: "copy"
    input:
        file "readqc_seperated/*" from ReadQCReport.collect()
        
    output:
        file "ReadQC.${params.readqc}/multiqc_report.html" into ReadQCReportOutput
        
    script:
        """
        multiqc --force --module fastp \
        		--dirs   readqc_seperated \
        		--outdir ReadQC.${params.readqc}
        """

}

process QuantReport{
	publishDir "${params.output}", mode:"copy"
	input:
		file "quant_seperated/*" from CellQuantReport.collect()

	output:
		file "Quant.${params.quant}/*" into CellQuantReportOutput

	script:
		"""
		multiqc --force --module salmon \
				--dirs quant_seperated \
				--outdir Quant.${params.quant}
		"""
}

if (params.use_email){
	workflow.onComplete {

	    def msg = """\
	        Pipeline execution summary
	        ---------------------------
	        Completed at: ${workflow.complete}
	        Duration    : ${workflow.duration}
	        Success     : ${workflow.success}
	        workDir     : ${workflow.workDir}
	        exit status : ${workflow.exitStatus}
			output      : ${params.output}
	        """
	        .stripIndent()

	    sendMail(to: "${params.email}", subject: 'My pipeline execution', body: msg)
	}	
}
