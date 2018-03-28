"""
Differential gene expression analysis with STAR and DESeq2.
See README.md for more information
"""

__author__ = "Colin Kern (https://github.com/kernco)"
__license__ = "MIT"

RAW_FILES, = glob_wildcards("Raw_Reads/{file}_1.fq.gz")

def comparisons(wildcards):
    import itertools
    groups = ['_'.join(x) for x in itertools.product(*config['factors'].values())]
    deseqfiles = []
    for comparison in config['comparisons']:
        fixed = set([x.replace(comparison[0], '').replace(comparison[1], '').replace('__', '_').strip('_') for x in groups])
        deseqfiles.extend(["DESeq/" + comparison[0] + "_vs_" + comparison[1] + "_" + x + ".tsv" for x in fixed])
    return deseqfiles    

def deseq_inputs(wildcards):
    import itertools
    f1 = wildcards.fixed.split('_') + [wildcards.condition]
    f2 = wildcards.fixed.split('_') + [wildcards.control]
    inputs = {}
    for group in itertools.product(*config['factors'].values()):
        if sorted(group) == sorted(f1):
            inputs['condition'] = 'Groups/' + '_'.join(group) + '.txt'
        elif sorted(group) == sorted(f2):
            inputs['control'] = 'Groups/' + '_'.join(group) + '.txt'
    return inputs


rule all:
    input:
        comparisons
        
rule star_index:
    input: 
        genome = config['genome'], 
        gtf = config['annotation']
    output: 
        'Genome_Index/chrName.txt'
    conda:
        'env.yaml'
    threads: config['threads']
    shell: 
        'mkdir -p Genome_Index STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir Genome_Index --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf}'

rule trim_paired_reads:
    input: 
        'Raw_Reads/{library}_1.fq.gz', 
        'Raw_Reads/{library}_2.fq.gz'
    output: 
        'Trimmed_Reads/{library}_1_val_1.fq.gz', 
        'Trimmed_Reads/{library}_2_val_2.fq.gz'
    conda:
        'env.yaml'
    shell: 
        'trim_galore -q 20 {input} -o Trimmed_Reads --paired'

rule star_align:
    input: 
        r1 = 'Trimmed_Reads/{library}_1_val_1.fq.gz', 
        r2 = 'Trimmed_Reads/{library}_2_val_2.fq.gz', 
        index = 'Genome_Index/chrName.txt'
    output: 
        'STAR_Output/{library}_Aligned.sortedByCoord.out.bam'
    conda:
        'env.yaml'
    threads: config['threads']
    shell: 
        'STAR --runThreadN {threads} --genomeDir Genome_Index --readFilesIn {input.r1} {input.r2} --readFilesCommand zcat --outFileNamePrefix STAR_Output/{wildcards.library}_ --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --outSAMtype BAM SortedByCoordinate'

rule filter_star:
    input: 
        rules.star_align.output
    output: 
        'Aligned_Reads/{library}.bam'
    conda:
        'env.yaml'
    shell: 
        'samtools view -b -q {config[mapq]} {input} > {output}'

rule make_alignment_report:
    input: 
        lambda wildcards: expand('Aligned_Reads/{library}.bam', library=RAW_FILES)
    output: 
        txt = 'Alignment_Summary.txt', 
        csv = 'Alignment_Summary.csv'
    params: 
        libraries = lambda wildcards: expand('{library}', library=RAW_FILES)
    threads: config['threads']
    script: 
        'Make_Alignment_Summary.py'

rule htseq_count:
    input:
        rules.filter_star.output
    output:
        'HTSeq/{library}.txt'
    conda:
        'env.yaml'
    shell:
        'htseq-count -f bam --stranded=reverse {input} {config[annotation]} > {output}'

rule make_group_files:
    input:
        expand('HTSeq/{library}.txt', library=RAW_FILES) #TODO: Only require files for group
    output:
        'Groups/{group}.txt'
    run:
        import itertools
        import os
        try:
            os.mkdir('Groups')
        except FileExistsError:
            pass
        dirlist = os.listdir('HTSeq')
        for group in itertools.product(*config['factors'].values()):
            files = ['HTSeq/' + x for x in dirlist]
            for factor in group:
                files = [x for x in files if factor in x]
            outfile = open('Groups/' + '_'.join(group) + '.txt', 'w')
            outfile.write('\n'.join(files) + '\n')
            outfile.close()
            
rule deseq:
    input:
        unpack(deseq_inputs)
    output:
        'DESeq/{condition}_vs_{control}_{fixed}.tsv'
    wildcard_constraints:
        control = "[^_]+"
    conda:
        'env.yaml'
    script:
        'deg.R'

