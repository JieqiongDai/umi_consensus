#!/bin/bash
mkdir -p log
sbatch wrapper/sbatch_snakemake.batch 2>log/sbatch.err
