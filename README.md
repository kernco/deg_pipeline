# DEG Pipeline

## Introduction

This pipeline performs gene expression quantification and identifies
differentially expressed genes from paired-end RNA-seq data. Currently,
single-end data is not supported. This pipeline uses STAR for read mapping
and DESeq2 for differentially expressed gene analysis.

## Installation and Requirements

### Requirements

The following software must be installed on your machine:

- [Python v3.x](https://www.python.org/downloads/)
- [Snakemake](http://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

If you have [Conda](https://conda.io/docs/) available, the pipeline can setup
and run all the necessary programs in a Conda environment. This is the
recommended way to execute this pipeline. Without Conda, the following are
also required:

- [Trim Galore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
- [Samtools](http://samtools.sourceforge.net/)
- [STAR](https://github.com/alexdobin/STAR)
- [HTSeq](https://htseq.readthedocs.io/en/release_0.9.1/install.html)
- [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

### Installation from ZIP file

Download a zip file from the Github [project page](https://github.com/kernco/deg_pipeline/archive/master.zip) and unzip it in
the desired location.

### Installation using git

If [git](https://git-scm.com/downloads) is available, execute the following
command while in the directory where you want to install the pipeline.

    git clone https://github.com/kernco/deg-pipeline.git

## Running the Pipeline

### Setting up the working directory

Create a working directory where all the input files, output files, logs, etc.
will be stored:

    mkdir MyProject
    cd MyProject

Within this directory, put your raw read files in fastq format into a
directory named "Raw_Reads":

    mkdir Raw_Reads
    mv ~/Path/To/Data/* Raw_Reads

If your read files end in .fastq.gz, rename them to .fq.gz:
 
    for file in Raw_Reads/*.fastq.gz; do mv $file ${file/.fastq.gz/.fq.gz}; done;

### Creating the configuration file

Copy the configuration template into your directory:

    cp ~/Path/To/Pipeline/config.yaml .

Open the file in a text editor such as nano or vim and edit the parameters to
match your genome, data, etc.
```
    genome      The path to the genome in fasta format. The
                pipeline will create a STAR index for the genome in the same
                directory as the fasta file, so the user running the pipeline
                needs to have write permission for the directory.
    annotation  The path to the gene annotation in GTF format. Currently, GFF
                files (such as those provided by NCBI) require modification
                of the htseq-count command in the Snakefile file of the pipeline
                to work properly.
    mapq        The MAPQ quality score cutoff used for filtering alignments.
    threads     The number of threads to use for tasks that can utilize multiple
                threads. Note that Snakemake will execute multiple tasks at the
                same time, and each task may use this many threads.
    factors     This is a list of the various factors that determine the grouping
                of your samples. These factors should match specific parts of
                the fastq filenames. See below for examples of this parameter.
    comparisons The comparisons you want to perform between groups. See below
                for examples of this parameter.
```

#### Simple example

We have 8 RNA-seq libraries from chicken liver, 4 from chickens infected with
avian influenza virus, and 4 from uninfected control birds. Our fastq files
are named with the bird ID and the group:

    Control4.fq.gz, Control8.fq.gz, Control15.fq.gz, Control23.fq.gz, Infected6.fq.gz, Infected9.fq.gz, Infected13.fq.gz, Infected27.fq.gz

We will define 1 factor, condition, in our configuration file:
```
    factors:
      condition: [Control, Infected]
```
The pipeline will know which fastq file belongs to which group by matching the
strings 'Control' and 'Infected' defined in the configuration file to the filenames.
Since we want to know the genes that are differentially expressed between the
control and infected groups, we put these as a comparison in the configuration file:
```
    comparisons:
      - [Infected, Control]
```
The order of the groups here will determine the meaning of the fold change in the
DEG table produced. A positive log2 fold change means higher expression in the
first group, while a negative log2 fold change means lower expression in the
first group, or higher expression in the second group.

#### A more complex example

We have now expanded our AIV experiment to include 2 tissues: liver and muscle. We
also are now using two breeds of chicken: a layer and a broiler.
We are interested in seeing the difference between control and infected for each
breed and each tissue, as well as the difference between the breeds.
Each file must be named so this information is contained in the filename. For
example, we could name our files as BirdID_Breed_Tissue_Condition.fq.gz, so we
might have the file 32_Broiler_Muscle_Infected.fq.gz and our configuration would
look like
```
    factors:
      condition: [Control, Infected]
      tissue: [Liver, Muscle]
      breed: [Layer, Broiler]
    comparisons:
      - [Infected, Control]
      - [Layer, Broiler]
```
The Infected-Control comparison will create 4 sets of differentially expressed
genes:

    - Layer Liver Infected vs. Layer Liver Control
    - Layer Muscle Infected vs. Layer Muscle Control
    - Broiler Liver Infected vs. Broiler Liver Control
    - Broiler Muscle Infected vs. Broiler Muscle Control

The Layer-Broiler comparison will create another 4 sets:

    - Layer Liver Control vs. Broiler Liver Control
    - Layer Muscle Control vs. Broiler Muscle Control
    - Layer Liver Infected vs. Broiler Liver Infected
    - Layer Muscle Infected vs. Broiler Muscle Infected

Note that DESeq2 is not run with a multifactor
design, rather these are just used to group the samples.

### Executing the Pipeline

Once your raw data files are in place and your configuration file is set up,
you can now execute the pipeline. It is recommended to first run the pipeline
in "dry run" mode which will only list the steps that will be execute without
running anything. This will help to ensure the configuration file is correct
and all the input files can be found:

    snakemake -s ~/Path/To/Pipeline/Snakefile --configfile config.yaml -n

If there are no errors, a lot will be printed. At the bottom of the output
will be a count of how many of each step will be run. You should check that
the number of "star_align" steps matches the number of samples you have, and
that the number of "deseq" steps matches the number of DEG lists you expect.
If your terminal does not support color output, you may need to add the '--nocolor'
option to the command to see all the output.

If everything looks good, you can run the pipeline by removing the '-n' from
the above command. If you are running the pipeline using Conda, add '--use-conda'
to the command.

Snakemake executes as many tasks in parallel as possible. To limit the number
of tasks being run at the same time, add the '-j' option, e.g. '-j 8'. This limits
the number of CPUs being used by Snakemake tasks to 8. If running in cluster
mode (see below), this limits the number of simultaneous jobs to 8.

If you are using a computing cluster, there is a way to run
the pipeline so that every task is submitted as a separate job to the cluster's
job queue. The run.sh script included with the pipeline contains an example
command to execute the pipeline this way, but it will likely need to be modified along
with the cluster.yaml file for your cluster's environment.

### Output Files

After executing the pipeline, a number of files and directories will be created
in the working directory.
```
    Alignment_Summary.txt   A text file containing a table of the read count,
                            alignment rate, etc. for each sample.
    Alignment_Summary.csv   The above table, but in a format that can be opened
                            more easily in Excel.
    DESeq/                  This folder contains tables from the DESeq analysis.
                            They are sorted by adjust p-value (q-value) so
                            a cutoff can be easily applied.
    HTSeq/                  This directory contains the raw expression values
                            for each sample.
    Aligned_Reads/          Contains BAM files of the filtered alignments.
    Trimmed_Reads/          FastQ files after trimming
    STAR_Output/            STAR logs and BAM files before filtering alignments.
    Groups/                 Text files showing the samples included in each group.
                            These can be useful for troubleshooting problems and
                            confirming that the correct comparisons were done.
    Logs/                   Logs from running each step, if running on a cluster. 
                            Look here for specific errors if the pipeline didn't work.
```


