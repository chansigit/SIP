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
//params.cellhome    = "/data/hca/chensijie/SMART1/SmallFASTQfiles"
//params.output      = "/data/hca/chensijie/SMART1_products"

// analysis step configs
params.readqc      = "fastp"
params.quant       = "salmon"
params.quant_unit  = "count"
params.downstream  = "seurat"



// notification options
params.use_email   = false
params.email       = "chensj16@mails.tsinghua.edu.cn"

// modifiable built-in options
/// Reference
params.salmon_index="/data/hca/chensijie/Reference/transcriptome/SALMON21_Homo_sapiens.GRCh38.cdna.all"
params.tx2gene = "/data/hca/chensijie/code/SIP/scripts/tx2gene.SALMON21_Homo_sapiens.GRCh38.cdna.all.tsv"


// Seurat Analysis Parameters
params.seuratMinCells = "1"
params.seuratMinGenes = "200"
params.seuratNGeneLowerBound = "200"
params.seuratNGeneUpperBound = "Inf"
params.seuratMitoLowerBound  = "-Inf"
params.seuratMitoUpperBound  = "0.30"
params.seuratSubgroupResolution = "1.03"
params.seuratTSNEPerplexity  = "4"
params.seuratDEGeneKept      = "10"
params.seuratDEGeneViz       = "3"
params.seuratDEMinDetectRatio= "0.15"
params.seuratDEpvalCut       = "0.05"

/// read pattern
params.PEpattern = "/**/*_{1,2}.fastq"
params.SEpattern = "/**/*.fastq"
PE_load_fail     = false


Channel
    .fromFilePairs(params.cellhome+params.PEpattern, size:-1)
    .ifEmpty{
        cprintln("black","yellow", "Paired-end reads undetected!")
        PE_load_fail=true
    }
    .into {ReadChannel}



if (PE_load_fail){
    ReadChannel.close()
    Channel
        .fromFilePairs(params.cellhome+params.SEpattern, size:-1)
        .ifEmpty{ 
            cprintln("black","yellow", "Single-end reads undetected!")
        }
        .into {ReadChannel}
}



process ReadQC{
    tag   {ongoing}
    cache true
    executor "local"
    input:
    set val(cellID), file(readPath) from ReadChannel

    output:
    stdout ReadQCSTDOUT
    file "*clean*.fastq" into CleanReadChannel
    file "*fastp*"       into ReadQCReport
    val   cellID         into CellIdentifierChannel

    script:
    def single= readPath instanceof Path
    ongoing="Ongoing: "+cellID
    if (!single){
        """
        echo Paired-end
        echo cellName=${cellID}
        echo read0=${readPath[0]} read1=${readPath[1]}

        fastp --thread 8 \
              --json ${cellID}_fastp.json \
              --html ${cellID}_fastp.html \
              --in1  ${readPath[0]} --in2 ${readPath[0]}\
              --out1 ${cellID}_clean_1.fastq  --out2 ${cellID}_clean_2.fastq 
        """
    }else{
        """
        echo Single-end
        echo cellName=${cellID}
        echo ${readPath[0]}

        fastp --thread 8 \
              --json ${cellID}_fastp.json \
              --html ${cellID}_fastp.html \
              --in1  ${readPath[0]} \
              --out1 ${cellID}_clean.fastq 
        """
    }

}




process GenerateReadQCReport{
    cache true
    publishDir "${params.output}", mode: "copy", overwrite:false
    input:
        file "readqc_seperated/*" from ReadQCReport.collect()
        
    output:
        file "ReadQC.${params.readqc}/*.html" into ReadQCReportOutput
        
    script:
        """
        multiqc --force --module fastp \
            --dirs   readqc_seperated/ \
            --outdir ReadQC.${params.readqc}
        cp -r readqc_seperated/*.html    ReadQC.${params.readqc}/
        """

}




/*
CleanReadChannel.subscribe{
    println it
}
*/

process Quant{
    tag {layout}
    cache true
    maxForks 48
    input:
    file readPath from CleanReadChannel
    val  cellID   from CellIdentifierChannel

    output:
    stdout Quantstdout
    file("*_quant*") into QuantResultDir

    script:
    def paired= readPath instanceof nextflow.util.BlankSeparatedList

    if (!paired){
        //cprintln("red","yellow", cellID)
        layout=cellID+" single-end"
        //println ">>"+cellID+": quantifying single-end reads"
        """
        salmon quant --index ${params.salmon_index} \
                     --libType A \
                     --unmatedReads ${readPath} \
                     --output ${cellID}_quant
        """
    }else{
        layout=cellID+" paired-end"
        //println ">>"+cellID+": quantifying paired-end reads"
        """
        salmon quant --index ${params.salmon_index} \
                     --libType A \
                     --mates1 ${readPath[0]} \
                     --mates2 ${readPath[1]} \
                     --output ${cellID}_quant
        """
    }

}

/* for debug 
Quantstdout.subscribe{
    println it

}
*/

QuantResQueue1 = Channel.create()
QuantResQueue2 = Channel.create()
QuantResultDir.separate(QuantResQueue1,QuantResQueue2){ path->
    [path,   path.toString()+"/quant.sf" ]
}

//QuantResQueue1.subscribe { println "Channel 1: $it" }
//QuantResQueue2.subscribe { println "Channel 2: $it" }

process GenerateQuantQCReport{
    cache      true
    publishDir "${params.output}", mode:"copy", overwrite:false
    input:
    file("quant_separated/*")  from QuantResQueue1.collect()

    output:
    file "Quant.${params.quant}/*" into QuantQCReportChannel
    stdout GQQRstdout

    script:
    //println x
    """
    mkdir -p Quant.${params.quant}/individuals
    cp -r quant_separated/*  Quant.${params.quant}/individuals

    multiqc  --force --module salmon   \
        --dirs   quant_separated       \
        --outdir Quant.${params.quant}
    """
}
//GQQRstdout.println()





process QuantMerge{
    cache      true
    publishDir "${params.output}", mode:"copy", overwrite:false
    input: 
        //it is important to declare the quantPaths as a Value input
        // otherwise nextflow will treat each * /quant.sf as a file and
        // re-assign some names for them, which is unwanted.
        val(quantPaths) from QuantResQueue2.collect()

    output:
        stdout QMstdout
        file "Quant.${params.quant}/PathListFile.txt"  into PathListFileChannel
        file "Quant.${params.quant}/DigitalGeneExpression.tsv"  into DGEChannel

    script:
        //println quantPaths.getClass()
        //for (String path in quantPaths){
        //    println path
        //}
        
        if (params.quant_unit=="count"){
            """
            mkdir -p Quant.${params.quant}/
            echo "${quantPaths}" > Quant.${params.quant}/PathListFile.txt
            Rscript ${workflow.projectDir}/scripts/QuantMerge.R \
                    Quant.${params.quant}/PathListFile.txt  \
                    Quant.${params.quant}/DigitalGeneExpression.tsv \
                    ${params.tx2gene}  \
                    count
            """  
        }
        
}
//QMstdout.println()


process SeuratAnalysis{
    cache      true   
    publishDir "${params.output}", mode:"copy", overwrite:false
    input:
        file "DigitalGeneExpression.tsv" from DGEChannel
    output:
        stdout SAstdout
        file "Downstream.${params.downstream}/CellQC.pdf"                  into Seurat_CellQC_Channel
        file "Downstream.${params.downstream}/MeanDispersionPlot.pdf"      into Seurat_MeanDispersionPlot_Channel
        file "Downstream.${params.downstream}/VizPCA.pdf"                  into Seurat_VizPCA_Channel
        file "Downstream.${params.downstream}/GroupInfo.tsv"               into Seurat_GroupInfo_Channel
        file "Downstream.${params.downstream}/tSNE.pdf"                    into Seurat_tSNE_Channel
        file "Downstream.${params.downstream}/SeuratObjectSerialized.Robj" into Seurat_SeuratObjectSerialized_Channel
        file "Downstream.${params.downstream}/DE_Genes.tsv"                into Seurat_DE_Genes_Channel
        file "Downstream.${params.downstream}/DE_Genes_ViolinPlot.pdf"     into Seurat_DE_Genes_ViolinPlot_Channel
        file "Downstream.${params.downstream}/DE_Genes_FeaturePlot.pdf"    into Seurat_DE_Genes_FeaturePlot_Channel
        file "Downstream.${params.downstream}/DE_Genes_Heatmap.pdf"        into Seurat_DE_Genes_Heatmap_Channel
    when :
        params.downstream == "seurat"    
    script:
        """
        mkdir -p  Downstream.${params.downstream}/
        Rscript ${workflow.projectDir}/scripts/SeuratAnalysis.R \
                --filepath="DigitalGeneExpression.tsv" \
                --output=Downstream.${params.downstream} \
                --mincells=${params.seuratMinCells} --mingenes=${params.seuratMinGenes} \
                --nGeneLowerBound=${params.seuratNGeneLowerBound} --nGeneUpperBound=${params.seuratNGeneUpperBound} \
                --mitoLowerBound=${params.seuratMitoLowerBound}   --mitoUpperBound=${params.seuratMitoUpperBound} \
                --resolution=${params.seuratSubgroupResolution} --tsneperp=${params.seuratTSNEPerplexity} \
                --DEGeneKept=${params.seuratDEGeneKept}  --DEGeneViz=${params.seuratDEGeneViz} \
                --DEminDetectRatio=${params.seuratDEMinDetectRatio} --DEpvalcut=${params.seuratDEpvalCut}
        """
}

SAstdout.println()

