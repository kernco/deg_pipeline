#Creates a table summarizing the alignment results.
#
#Author: Colin Kern (https://github.com/kernco)
#License: MIT

import subprocess
import multiprocessing
import os

def get_data(library):
    if os.path.isfile("Trimmed_Reads/{}_2.fq.gz_trimming_report.txt".format(library)):
        reportname = "Trimmed_Reads/{}_2.fq.gz_trimming_report.txt".format(library)
        divfactor = 2
    elif os.path.isfile("Trimmed_Reads/{}.fq.gz_trimming_report.txt".format(library)):
        reportname = "Trimmed_Reads/{}.fq.gz_trimming_report.txt".format(library)
        divfactor = 1
    with open(reportname) as f:
        for line in f:
            if line.startswith("Total reads processed"):
                rawreads = int(line.split()[3].replace(',',''))
    trimreads = subprocess.check_output(["grep", "Number of input reads", "STAR_Output/{}_Log.final.out".format(library)])
    trimreads = int(trimreads.split(b'|')[1])
    uniquely = subprocess.check_output(["grep", "Uniquely mapped reads number", "STAR_Output/{}_Log.final.out".format(library)])
    multimap = subprocess.check_output(["grep", "Number of reads mapped to multiple loci", "STAR_Output/{}_Log.final.out".format(library)])
    alignpairs = int(uniquely.split(b'|')[1]) + int(multimap.split(b'|')[1])
    filtered = int(int(subprocess.check_output(["samtools", "view", "-c", "Aligned_Reads/{}.bam".format(library)])) / divfactor)
    return [library, rawreads, trimreads, (float(trimreads) / rawreads)*100, alignpairs, (float(alignpairs) / trimreads)*100, filtered, (float(filtered) / alignpairs)*100, (float(filtered) / rawreads)*100]

with multiprocessing.Pool(processes=snakemake.threads) as pool:
    data = list(pool.imap(get_data, snakemake.params.libraries))

with open(snakemake.output.txt, 'w') as f:
    f.write("{: >15s} {: >14s} {: >23s} {: >23s} {: >23s} {: >9s}\n".format("Library", "Raw Reads", "Trimmed Reads", "Aligned Reads", "Filtered Reads", "Total %"))
    for line in data:
        print(line)
        f.write("{: >15s} {: >14,d} {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >14,d} ({: >5.1f}%) {: >8.1f}%\n".format(*line))

with open(snakemake.output.csv, 'w') as f:
    headers = ["Library", "Raw Reads", "Trimmed Reads", "%", "Aligned Reads", "%", "Filtered Reads", "%", "Total %"]
    f.write(','.join(headers) + '\n')
    for line in data:
        f.write(','.join([str(x) for x in line]) + '\n')
