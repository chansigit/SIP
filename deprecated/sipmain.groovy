cellHome = "/data/hca/chensijie/SMART1/SmallFASTQfiles"
//productHome = "/data/hca/chensijie/SMART1_products"



params.readqc = 'fastqc'


cellList = {  cellHome ->
    cell_list = []
	file(cellHome).eachDir(){dir -> cell_list << dir.toAbsolutePath()}    
	cell_list
}


/*
****************************************
       Groovy color print utilities
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




def isPairedEnd(dir){
	Channel.fromFilePairs("${dir}/*_{1,2}.fastq").count().subscribe{ 
		if (it>0) {return true} else {return false} 
	}
}

cellFlow = Channel.from(cellList(cellHome))
//cellFlow.println()


process ReadQualityControl{
	input:
	file cellName from cellFlow
	
	output:
	val "ReadQualityControl.${params.readqc}/$cellName" into ReadQCReport

	script:
	if (params.readqc=="fastqc")
		"""
		mkdir -p ReadQualityControl.${params.readqc}/$cellName

		fastqc --threads 8     --format fastq \
	           --outdir  ReadQualityControl.${params.readqc}/$cellName \
	           $cellHome/$cellName/*.fastq
	    
		"""
	else
		"""
		printf fastp
		"""


}

ReadQCReport.subscribe{
	cprintln(33,42, it+" done!")
}



process ReadQCMerging{
	input:
	val "ReadQualityControl.${params.readqc}"

	script:
	"""
	multiqc --dirs   ReadQualityControl.${params.readqc} \
	        --outdir ReadQualityControl.${params.readqc} \
	"""
}