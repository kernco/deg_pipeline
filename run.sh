#!/bin/bash
# Runs the pipeline on a cluster that uses the SLURM scheduler.
# Check cluster.yaml for cluster configuration options.

module load bio3

mkdir -p Logs
snakemake -j 999 --cluster-config /home/ckern/DEG_Pipeline/cluster.yaml --cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -n {cluster.cpus} -J {cluster.name} -o {cluster.output} -e {cluster.output}" -s /home/ckern/DEG_Pipeline/Snakefile --configfile config.yaml --latency-wait 120 -p -k --use-conda $@

