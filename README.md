# Temporal dynamics of population structure and genetic diversity in rye

This repository contains the code, analysis outputs, and supporting files generated for the Master's thesis:

**Temporal dynamics of population structure and genetic diversity in rye in response to nutrient deficiency – Insights from a three-generation experiment**

## Project overview

This project investigates short-term population-genomic changes in rye (*Secale cereale*) across three successive generations under contrasting fertilization conditions.

Generations 1–3 were grown in the Eternal Rye field experiment in Halle (Saale) under two treatments:
- **ND:** unfertilized conditions on long-term nutrient-depleted soil
- **control:** mineral nitrogen, phosphorus, and potassium fertilization
An independently established Gen0 population was used as an additional reference population.

The main aim was to investigate whether genomic changes can be detected over only a few generations and whether these changes are primarily associated with generation or fertilization treatment. Temporal allele-frequency changes were additionally compared with neutral Wright–Fisher drift expectations.

The analyses include:

- NGS quality control, read trimming, alignment, variant calling, and variant filtering
- population structure analysis using PCA and ADMIXTURE
- nucleotide diversity, Tajima's D, and FST analysis
- genome-wide temporal allele-frequency analysis
- Wright–Fisher drift simulations
- SNP-wise modelling of allele-frequency trajectories
- treatment-specific tests of candidate SNP trajectories
- functional annotation and GO enrichment analysis

---

## Repository structure

```
population_dynamics_in_rye/
│
├── code/
│   ├── ngs_main.sbatch
│   │
│   ├── ngs_pipeline/
│   │   ├── fastqc_module.sbatch
│   │   ├── trimming_module.sbatch
│   │   ├── alignment_module.sbatch
│   │   ├── variant_calling_module.sbatch
│   │   ├── variant_filtering_module.sbatch
│   │   └── workers/
│   │
│   └── R_scripts/
│
├── additional_data/
│   └── id_data/
│       └── samples_id_gen.csv
│
├── general_analysis/
│   ├── NGS_QC/
│   ├── PCA/
│   │   ├── indiv_pcadapt_PCA/
│   │   └── gen0_diff_PCA/
│   ├── Admixture/
│   ├── Fst_StAMPP/
│   ├── vcftools/
│   │   ├── windowed_pi/
│   │   ├── tajimaD/
│   │   ├── Fst/
│   │   └── freq/
│   │       └── logistic_trend/
│   └── annotation_and_GO/
├── R_packages.txt
└── README.md
```

### `code/`

Contains the scripts used for data processing and statistical analysis.

#### `code/ngs_main.sbatch`

Central configuration and submission script for the SLURM-based NGS processing pipeline. Dataset-specific paths, enabled modules, analysis parameters, and computational resources are configured in this file.

#### `code/ngs_pipeline/`

Contains the individual pipeline modules for:

- initial FastQC
- trimming
- alignment
- variant calling
- variant filtering

The `workers/` subdirectory contains the SLURM worker scripts performing the individual computational steps.

#### `code/R_scripts/`

Contains the R scripts used for the population-genomic analyses and visualisations described in the thesis.

These include scripts for PCA, ADMIXTURE visualisation, diversity statistics, FST analysis, allele-frequency analysis, Wright–Fisher simulations, logistic modelling of SNP trajectories, and functional annotation.

### `additional_data/`

Contains additional files required by the analysis scripts.

The `id_data/` directory contains sample ID lists and the assignment of individuals to generations and treatments.

### `general_analysis/`

Contains selected analysis outputs, result tables, and figures.

Important subdirectories include:

| Directory | Content |
| --- | --- |
| `NGS_QC/` | MultiQC reports and NGS preprocessing summaries |
| `PCA/indiv_pcadapt_PCA/` | Individual-level PCA results |
| `PCA/gen0_diff_PCA/` | PCA of allele-frequency differences relative to Gen0 |
| `Admixture/` | ADMIXTURE results and visualisations |
| `Fst_StAMPP/` | Pairwise FST analysis using StAMPP |
| `vcftools/windowed_pi/` | Nucleotide diversity results |
| `vcftools/tajimaD/` | Tajima's D results |
| `vcftools/Fst/` | Window-based FST results |
| `vcftools/freq/` | Allele-frequency and MAF analyses |
| `vcftools/freq/logistic_trend/` | SNP-wise allele-frequency trajectory modelling |
| `annotation_and_GO/` | Functional annotation and GO enrichment results |

### `R_packages.txt`

Lists the R packages and package versions loaded across the downstream analysis scripts. The analyses were performed and tested using R 4.3.3.

---

## External data required for individual R scripts

Most intermediate and result files required by the R analyses are included in this repository. Some large or externally hosted files are not included.

### Large VCF file available on request

The following VCF file is not stored in the repository because of its size:
```
LD-pruned_minDP4_maxDP100_maf0.002_maxMiss0.9_minQ40_noIndels_biallelic_noPopVar.vcf.gz
```

It is required by:
```
code/R_scripts/Fst_analysis_StAMPP.R
```

The file is available on request.

When running the script, provide the local path to the VCF file when requested.

### Lo7 v3 annotation files available from Zenodo

Additional Lo7 v3 genome annotation files can be downloaded from:

https://zenodo.org/records/15801361

The following files are required:

| File | Required by |
| --- | --- |
| `SECALE.CEREALE.Lo7V3.pgsb.r2.Feb2025.featuretable.csv` | `annotation_sig_SNPs.R`, `GO_term_enrichment.R` |
| `SECALE.CEREALE.Lo7V3.pgsb.r2.Feb2025_all.aa.fa` | `find_missense_prot_seqs.R` |

After downloading the files, provide their local paths when requested by the respective R scripts.

---


## NGS processing pipeline

This repository contains a reusable SLURM-based pipeline developed for the processing and quality control of single-end sequencing data on the **IPK Gatersleben SLURM cluster**.
The pipeline configuration, SLURM directives, resource settings, and job-submission logic were designed for the IPK computing environment. The pipeline can in principle be adapted to other SLURM-based clusters, but cluster-specific settings and software availability need to be reviewed before use outside the IPK environment.

The central configuration and submission file is:

```
code/ngs_main.sbatch
```

This script does not perform the computationally intensive processing itself. Instead, it configures and submits the corresponding module and worker jobs.

### Available modules

Modules are activated in `ngs_main.sbatch` using `0` or `1`:
```
INITIAL_FASTQC=0
INITIAL_MULTIQC=0
TRIMMING=0
ALIGNMENT=0
VARIANT_CALLING=0
VARIANT_FILTERING=0
```

| Module | Function |
| --- | --- |
| `INITIAL_FASTQC` | Run FastQC on input FASTQ files and summarize reports with MultiQC |
| `INITIAL_MULTIQC` | Generate a MultiQC report from existing QC reports |
| `TRIMMING` | Trim reads with Trim Galore and generate QC reports |
| `ALIGNMENT` | Align reads with minimap2, filter/sort/index BAM files, and generate mapping statistics |
| `VARIANT_CALLING` | Perform variant calling using bcftools |
| `VARIANT_FILTERING` | Filter an existing VCF using VCFtools and optionally generate additional QC statistics |

### Module dependencies

The main processing modules can be chained automatically in the following order:

```
TRIMMING
   ↓
ALIGNMENT
   ↓
VARIANT_CALLING
```

If several of these modules are enabled in the same run, SLURM `afterok` dependencies are used. A downstream step is therefore submitted only after successful completion of the preceding module.

The output of the previous module is automatically used as input for the next step.

`INITIAL_FASTQC`, `INITIAL_MULTIQC`, and `VARIANT_FILTERING` are independent of this dependency chain.

For inspection of intermediate results and quality control between processing stages, running the main modules separately is recommended.

---

## Running the NGS pipeline

### 1. Configure the repository path

Set `CODE_DIR` to the `ngs_pipeline` directory of the cloned repository:

```
CODE_DIR="/path/to/population_dynamics_in_rye/code/ngs_pipeline"
```

### 2. Configure SLURM settings

At the beginning of `ngs_main.sbatch`, adjust the cluster-specific SLURM settings, particularly:

```
#SBATCH --output=/path/to/job/logs/submit_ngs_%j.out
#SBATCH --error=/path/to/job/logs/submit_ngs_%j.err
#SBATCH --mail-user=your@email.address
```

Also set:

```
MAIL_USER="your.email@example.org"
```

`MAIL_USER` is passed to the module jobs submitted by the main script.

### 3. Select the modules

Enable the required modules by changing the corresponding flags from `0` to `1`.
For example:
```
TRIMMING=1
ALIGNMENT=1
VARIANT_CALLING=1
```
This performs a chained trimming → alignment → variant-calling run.

### 4. Set dataset-specific paths

The following settings should be checked for every run:

```
JOB_NAME="my_dataset"

ALL_OUT="/path/to/output"

INPUT="/path/to/input"

REFERENCE_GENOME="/path/to/reference/genome.fasta"

INDEX="/path/to/reference/genome.fasta.mmi"
```

The file-selection patterns should also be adapted where necessary:
```
GLOB_PATTERN_FASTQ=(
    "${INPUT}/*.fastq"
    "${INPUT}/*.fastq.gz"
    "${INPUT}/*.fq"
    "${INPUT}/*.fq.gz"
)
```

For variant calling from existing BAM files:
```
GLOB_PATTERN_BAM=(
    "${INPUT}/*.bam"
)
```

For standalone variant filtering:
```
VCF_IN="${INPUT}/input.vcf.gz"
```

### 5. Adjust analysis parameters

#### Trimming

| Parameter | Default | Description |
| --- | ---: | --- |
| `QUALITY` | `20` | Minimum Phred quality used for read-end trimming |
| `MIN_LENGTH` | `20` | Minimum retained read length |
| `CLIP3` | `0` | Number of bases removed from the 3' end |
| `CLIP5` | `0` | Number of bases removed from the 5' end |
| `ADAPTER` | `""` | Adapter sequence; an empty value allows automatic adapter detection |

#### Alignment

| Parameter | Default | Description |
| --- | ---: | --- |
| `MAP_QUAL` | `30` | Minimum mapping quality of retained alignments |
| `MERGING_FLAG` | `0` | Merge generated BAM files (`0 = no`, `1 = yes`) |

#### Variant calling

| Parameter | Default | Description |
| --- | ---: | --- |
| `VAR_MAP_QUAL` | `30` | Minimum MAPQ retained after filtering |
| `BASE_QUAL` | `30` | Minimum base quality |
| `ANNOT_TAGS` | `"DP,AD"` | INFO/FORMAT annotations generated during `bcftools mpileup` |
| `VARIANT_QUAL` | `40` | Minimum variant quality |
| `MAX_BASE_DEPTH` | `260000` | Maximum total depth allowed for a variant |
| `MIN_BASE_DEPTH` | `10` | Minimum total depth required for a variant |
| `MIN_ALT_DEPTH` | `0` | Minimum alternative-allele depth |
| `REGION` | `1` | Run variant calling separately by chromosome/region (`0 = no`, `1 = yes`) |

When region-wise calling is enabled, variant calling is parallelized and the resulting VCF files are concatenated automatically.

#### Variant filtering

| Parameter | Default | Description |
| --- | ---: | --- |
| `MIN_DP` | `4` | Minimum genotype depth |
| `MAX_DP` | `100` | Maximum genotype depth |
| `MAF` | `0` | Minimum minor allele frequency threshold |
| `MAX_MISSING` | `0` | Minimum proportion of non-missing genotypes required per site |
| `MIN_Q` | `40` | Minimum variant quality |
| `RM_INDELS` | `1` | Remove indels (`0 = no`, `1 = yes`) |
| `RM_MULTIALLELIC` | `1` | Remove multiallelic variants (`0 = no`, `1 = yes`) |

Additional VCFtools QC statistics can be enabled with the following settings:

| Parameter | Default | Description |
| --- | ---: | --- |
| `QC_STATS_FLAG` | `1` | Generate additional VCFtools QC statistics (`0 = no`, `1 = yes`) |
| `TAJIMAD_WINDOW` | `10000` | Window size in base pairs for Tajima's D |
| `PI_WINDOW` | `10000` | Window size in base pairs for nucleotide diversity |

### 6. Configure computational resources

CPU allocation, memory, and runtime can be adjusted separately for each processing step in the **Optional technical settings** section of `ngs_main.sbatch`.
The maximum number of simultaneously active array tasks is controlled by:
```
ARRAY_CAP=25
```
This can be adjusted according to the available cluster resources.

### 7. Submit the pipeline

After configuration, submit the main script with:
```
sbatch code/ngs_main.sbatch
```
Output directories are generated from `JOB_NAME` and the selected module, for example:
```
<ALL_OUT>/<JOB_NAME>_trimming/
<ALL_OUT>/<JOB_NAME>_alignment/
<ALL_OUT>/<JOB_NAME>_variant_calling/
```

SLURM log files are stored in:
```
<ALL_OUT>/job/
```

---

## Software

The NGS pipeline uses the following main software:
- FastQC
- MultiQC
- Trim Galore / Cutadapt
- minimap2
- samtools
- bcftools
- VCFtools
- SLURM

The downstream analyses were performed and tested using **R 4.3.3**. The R packages and package versions loaded across the analysis scripts are listed in [`R_packages.txt`](R_packages.txt).

Software versions and detailed parameters of the NGS processing pipeline and downstream analyses are documented in the thesis.

---

## Reproducibility and data availability

This repository contains the analysis code and selected intermediate and final results required to document the analyses performed in the thesis. Large primary sequencing files are not stored in the repository.

The large LD-pruned VCF required for `Fst_analysis_StAMPP.R` is available on request. Publicly available Lo7 v3 annotation resources required for the functional analyses can be obtained from the Zenodo record listed above.

---

## Citation

If you use code or results from this repository, please cite the associated Master's thesis:

> Weise, Annemarie Bianka (2026). *Temporal dynamics of population structure and genetic diversity in rye in response to nutrient deficiency – Insights from a three-generation experiment*. Master's thesis, Martin Luther University Halle-Wittenberg.

For questions regarding the repository or access to raw data, please contact the repository author.
