# Fastq - QC, subsampling, and consensus

##### Author: Jieqiong Dai

## Description
This snakemake pipeline is for Next Generation Sequencing reads QC, trimming, subsampling, and UMI based consensus. The pipeline may be run on an HPC or in a local environment.

Major functions in this pipeline include:
* Quality trimming and QC
* Optional pre-trimming QC
* Optional subsampling
* Optional consensus

The final deliverables include:
* A multiQC report of all QC metrics with all samples
* Processed fastq files

## Software requirement
* [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Fastp](https://github.com/OpenGene/fastp)
* [FastqQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [STAR](https://github.com/alexdobin/STAR)
* [BWA](http://bio-bwa.sourceforge.net/)
* [Fgbio](https://github.com/fulcrumgenomics/fgbio)
* [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
* [Seqtk](https://github.com/lh3/seqtk)
* [MultiQC](https://multiqc.info/)
* [Samtools](http://www.htslib.org/)

## User's guide
### I. Input requirement
* Edited config/config.yaml
* Demultiplexed fastq files
* Reference genome file and annotation files (optional, required for consensus)
* STAR index (optional, required for RNAseq data consensus)
* BWA index (optional, required for DNAseq data consensus)
* GATK4 reference genome sequence dictionary (optional, required for DNAseq data consensus)

### II. Editing the config.yaml
#### Basic (need to be modified case by case)
* `run_id`: Run ID
* `build`: Reference genome build, eg. hg19, hg38
* `data_type`: Data Type, DNA or RNA
* `reference_bwa`: Path to reference genome fasta file (within the BWA index directory), required for DNAseq data consensus
* `reference_gatk`: Path to reference genome fasta file (within the GATK sequence dictionary directory), required for DNAseq data consensus
* `gtf`: Path to reference genome annotation GTF file, required for RNAseq data consensus
* `star_indice`: Path to STAR index, required for RNAseq data consensus
* `fastq`: Path to the merged fastq files stored directory.
* `name_split`: Sample ID delimitation. eg. _R, _sub_R
* `fastq_name`: Naming format of fastq files. eg. {sample}, {sample}_sub
* `out`: Path to the output directory. Default output/ in working directory
* `umi`: If used KAPA UMI , Y or N
* `umi_structure`: UMI read structure, eg, "3M3S+T 3M3S+T" for KAPA UMI, leave it blank if consensus is not required
* `subsample`: If needs subsampling, Y or N
* `read_number`: Number of read pairs in subsampling, leave it blank if subsampling is not required
* `extract_umi`: Fgbio ExtractUmisFromBam parameters, required for consensus
* `sam2fastq`: GATK4 SamToFastq parameters, required for consensus
* `trim`: Fastp trimming parameters
* `star`: STAR alingment parameters, required for RNA consensus
* `bwa`: BWA mem alignment parameters, required for DNA consensus
* `consensus`: Fgbio CallMolecularConsensusReads parameters, required for consensus

#### Optional
* `pre_qc`: Optional pre-trimming QC, Y or N, default N
* `fastq2sam`: GATK4 FastqToSam parameters
* `groupumi`: Fgbio GroupReadsByUmi parameters

### III. To run
#### Clone repository
* Clone the repository to your working directory
```bash
git clone git@github.com:JieqiongDai/umi_consensus.git
```

#### Environment
* Install Conda if not pre-installed
* If not pre-installed, create a conda environment in the appropriate directory using the provided yaml file from the git cloned directory `workflow/envs/` and activate it after installation:
```bash
conda env create -f tidyfastq_env.yaml
conda activate tidyfastq
```

#### Running the code
* Edit and save `config/config.yaml`
* To run on an HPC using slurm job scheduler:
  Edit `config/cluster_slurm_config.yaml` according to your HPC information;
  Run `wrapper/run_sbatch.sh` to initiate running of the pipeline
```bash
bash wrapper/run_sbatch.sh
```
* To run on an HPC using SGE/UGE job scheduler:
  Edit `config/cluster_SGE_config.yaml` according to your HPC information;
  Run `wrapper/run_qsub.sh` to initiate running of the pipeline
```bash
bash wrapper/run_qsub.sh
```
* To run on a local server:
  Run `wrapper/snakemake.batch` to initiate running of the pipeline
```bash
bash wrapper/snakemake.batch
```
* Look in log directory for logs for each rule
* To view the snakemake rule graph:
```bash
snakemake --rulegraph | dot -T png > tidyfastq.png
```

### IV. Example output
```bash
.  user/defined/output_dir
├── align (optional)
│   ├── {sample_A}
│   ├── {sample_B}
│   ├── {sample_C}
├── archive
│   ├── config.yaml
│   └── {run_id}_{build}_multiqc_report.html
├── consensus (optional)
│   ├── fastq
│   └── ...
├── extract_umi (optional)
├── final_fastq (the final fastq files for downstream processing)
├── multiqc
├── qc
├── subsample_fastq (optional)
└── trimmed
```

### Special Notes:
* In UMI processing, the consensus reads in the final fastq files are filtered and does NOT contain the UMI information. The UMI information is stored in the bam file: `{output_dir}/consensus/{sample}/{sample}_umi_consensus_unmapped.bam`
