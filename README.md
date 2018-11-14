# Quick start

SIP is built based on Nextflow framework, thus requiring POSIX-compatible operating system (Linux, UNIX, MacOS X, or Windows subsystem for Linux). Apart form the OS requirements, SIP needs Java environment and Java version should be 1.8 or above.

SIP supports analyzing both under local environment and in docker containers. In the local mode, users should pre-install all needed bioinformatics tools and configure the environment variables properly so that these tools are available at command lines. In the container mode, users should download the docker image we prepared and set up a container so that SIP can execute analysis commands with tools prepared in the container image.

## 0. Basic concepts

- Cell Home Directory: directory holding all raw sequencing fastq files
- Output Directory: final analysis results are placed here
- Working directory: intermediate results are placed here

## 1. Installation
Enter the working directory, run the following command to load Nextflow framework from internet.

```bash
curl -s https://get.nextflow.io | bash
```

Download SIP's code in the Working directory:

```bash
git clone https://github.com/chansigit/SIP.git
```

## 2. Prepare data

Please put all data to be analyzed in the Cell Home Directory, and arrange them in the following way:

Suppose the path to Cell Home Directory is `/data/TestDataset1` , users need to create multiple subfolders under `/data/TestDataset1` and name them after each cell's identifier. 

Each subfolder contains the raw sequencing fastq file(s) of each cell. For single-end experiments, the fastq file's name should be `{CellName}.fastq`; while for paired-end experiments, the fastq file pairs should be named as `{CellName}_1.fastq` and `{CellName}_2.fastq`.

As a summary, the arrangement of reads shoule be `/data/TestDataset1/CellName/Cell1Name.fastq` for a single-end experiment; and `/data/TestDataset1/CellName/Cell1Name_1.fastq` , `/data/TestDataset1/CellName/Cell1Name_2.fastq` for a paired-end experiment.



## 3. Run SIP

We choose `fastp` as the read quality control tool, `salmon` as the gene-expression quantification tool, and use *read count* as the quantification unit. We used the *Seurat* package for downstream analysis after quantification. 

```bash
./nextflow run SIP.nf  -resume \
        --readqc=fastp --quant=salmon --quantunit=count \
        --downstream=seurat \
        --cellhome={Cell Home Directory Path} \
        --output={Output Directory Path}
```

