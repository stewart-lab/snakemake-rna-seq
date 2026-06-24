# snakemake-rna-seq

A [Snakemake](https://snakemake.github.io/) pipeline for bulk paired-end RNA-seq:
from raw FASTQs to differential expression, gene-set enrichment, QC reports, and
an auto-generated methods paragraph. Each step runs in its own pinned conda
environment.

## What it does

Trim → align/quantify → QC → build count matrix → filter → explore (PCA) →
differential expression → gene-set analysis (GSVA **and** GSEA) → reporting.

| Step | Output folder | Tool | Description |
|------|---------------|------|-------------|
| 1 | `001_fastp_trimmed` | fastp | Adapter/quality trimming + polyX removal |
| 2 | `002_fastqc` | FastQC | QC of the trimmed reads |
| – | `<genome_dir>/star_reference` | RSEM + STAR | Build the alignment/quantification index |
| – | `<genome_dir>/picard_reference` | Picard | Build the sequence dictionary, refFlat, and rRNA intervals |
| 3 | `003_rsem` | RSEM + STAR | Align and quantify gene-level expression |
| 4 | `004_qc_reports` | Picard + RSEM/STAR | Per-sample QC "report card" (mapping %, reads, genes detected, rRNA, exonic, intergenic, 5'/3' bias) |
| 5 | `005_multiqc` | MultiQC | Aggregate fastp/FastQC/RSEM/Picard QC into one report |
| 6 | `006_raw_counts` | tximport | Assemble the per-sample results into a count matrix |
| 7 | `007_raw_counts_samples_dropped` | pandas | Drop excluded samples (see `exclude.csv`) |
| 8 | `008_raw_counts_filtered` | pandas | Low-count gene filter (all-samples **and** dropped matrices) |
| 9 | `009_pca` | DESeq2/Plotly | VST → PCA + sample clustermap, all samples |
| 10 | `010_pca_samples_dropped` | DESeq2/Plotly | VST → PCA + sample clustermap, samples dropped |
| 11 | `011_deseq2` | DESeq2 | Differential expression per contrast |
| 12 | `012_deseq2_gene_symbols` | gffutils | Map gene IDs → gene names using the GTF |
| 13 | `013_gsva` | GSVA + limma | Per-sample gene-set scores + differential enrichment |
| 14 | `014_gsea` | fgsea | Ranked-list gene-set enrichment per contrast |
| 15 | `015_prose` | – | Auto-generated methods paragraph + tool versions |
| – | `reports` | Quarto | Render any user-supplied `.qmd` reports |

## Running

The pipeline is packaged as a Docker image. `app/main.py` stages the workflow
next to your config (it creates a sibling `workflow/` directory containing the
`Snakefile`, `envs/`, `rules/`, and `scripts/`) and invokes Snakemake with
`--use-conda`.

```bash
# Mount a working directory that contains your config + inputs (and your
# reference directory), then point main.py at the config file.
# Tip: give the mounted directory the same path it has on your filesystem, so
# the paths in your config stay valid for collaborators.
docker run --rm \
  -v /w5home/myname/mystuff:/w5home/myname/mystuff \
  snakemake-rna-seq \
  python main.py \
  /w5home/myname/mystuff/myproject/config/config.yaml
```

Snakemake creates and caches each step's conda environment on first run, then
builds the targets defined in `rule all`. Re-running only rebuilds what changed.

### Suggested project directory structure

`main.py` treats the **parent of the config file's directory** as the project
root and stages `workflow/` there. A tidy layout that keeps inputs, config, and
outputs separate looks like this:

```text
myproject/
├── config/                # the run configuration + sample/design sheets
│   ├── config.yaml         #   ← pass THIS path to main.py
│   ├── samples.csv
│   ├── design.csv
│   ├── contrasts.csv
│   ├── gene_sets.csv
│   └── exclude.csv
├── resources/             # all inputs the pipeline reads
│   ├── reads/              #   raw FASTQs referenced by samples.csv
│   ├── reference/          #   genome_dir: exactly one FASTA + one GTF
│   │                       #   (the star_reference/ and picard_reference/
│   │                       #    indexes are built in here on first run)
│   └── gene_sets/          #   .gmt files referenced by gene_sets.csv
├── reports/               # reports_dir: your Quarto reports, one folder each
│   └── summary/summary.qmd
├── workflow/              # auto-staged by main.py (Snakefile + envs/rules/scripts)
├── results/              # results_dir: all pipeline outputs
└── logs/                 # logs_dir: per-rule logs
```

You are free to name and arrange these however you like — the pipeline only
cares about the paths you set in `config.yaml`. **Absolute paths are
recommended** (see below), since the example above is just a convention.

## Inputs

All inputs are referenced from a single config file (see
[`example/config.yaml`](example/config.yaml)). The CSV examples below live in
[`example/`](example/) and can be used as templates.

### `config.yaml`

| Key | Meaning |
|-----|---------|
| `samples` | Path to the sample sheet CSV |
| `design` | Path to the design CSV |
| `contrasts` | Path to the contrasts CSV |
| `gene_sets` | Path to the gene-sets CSV (lists GMT files) |
| `exclude` | Path to the exclusions CSV (optional; remove the key to disable) |
| `reports_dir` | Directory of user-supplied Quarto `.qmd` reports to render (optional; omit to skip — see [Reports](#reports--quarto-optional)) |
| `genome_dir` | Directory holding exactly one FASTA and one GTF |
| `filter.min_count` / `filter.min_samples` | Low-count gene filter thresholds |
| `cores` | Number of CPU cores for the pipeline (RSEM/STAR threads, etc.) |
| `results_dir` / `logs_dir` | Output locations |

Paths may be absolute or relative to the directory the pipeline runs in.
Absolute paths are recommended so your colleagues can find the data
if you send them your config files.

### Sample sheet — [`example/samples.csv`](example/samples.csv)

One row per sample. FASTQs may be `.fastq` or `.fastq.gz`.

```csv
sample,r1,r2
WT_rep1,reads/WT_rep1_R1.fastq.gz,reads/WT_rep1_R2.fastq.gz
KO_rep1,reads/KO_rep1_R1.fastq.gz,reads/KO_rep1_R2.fastq.gz
```

### Design — [`example/design.csv`](example/design.csv)

Maps each sample to a condition (used as the DESeq2 design `~ condition`).

```csv
sample,condition
WT_rep1,WT
KO_rep1,KO
```

### Contrasts — [`example/contrasts.csv`](example/contrasts.csv)

Each row is one comparison; results are named `<treatment>_vs_<control>`.

```csv
control,treatment
WT,KO
```

### Gene sets — [`example/gene_sets.csv`](example/gene_sets.csv)

A single `gmt` column listing GMT files. All listed GMTs are merged for GSVA and
GSEA. Gene sets must use **gene symbols** (matching the GTF-derived names).
Download collections such as the [MSigDB Hallmark sets](https://www.gsea-msigdb.org/gsea/msigdb/)
as `.gmt` files. A tiny illustrative GMT is in
[`example/gene_sets/example_sets.gmt`](example/gene_sets/example_sets.gmt).

```csv
gmt
example/gene_sets/example_sets.gmt
```

### Exclusions — [`example/exclude.csv`](example/exclude.csv) (optional)

Samples to drop before differential expression (e.g. PCA outliers). Every listed
sample must exist in the sample sheet. The `reason` is
included in the methods paragraph.

```csv
sample,reason
KO_rep3,it had a low number of reads
```

### Reference (`genome_dir`)

Point `genome_dir` at a directory that contains **exactly one** genome FASTA
(`.fa`/`.fasta`/`.fna`) and **exactly one** annotation GTF (`.gtf`) — both are
auto-detected. The pipeline builds the RSEM/STAR index (`star_reference/`) and
the Picard reference artifacts (`picard_reference/`: sequence dictionary,
refFlat, rRNA intervals) inside this directory. Gene names are read from the
GTF (`gene_name`/`gene`/`gene_symbol`/`Name` attributes), so the pipeline is not
tied to any single organism.

### Reports — Quarto (optional)

Point `reports_dir` at a directory of your own Quarto reports. Each `.qmd` is
expected one level down, i.e. `<reports_dir>/<name>/<name>.qmd`. They are
rendered **last**, after the rest of the pipeline finishes (so a report can read
any pipeline output), and the rendered HTML lands under
`<results_dir>/reports/<name>/<name>.html`. This step is skipped entirely if
`reports_dir` is omitted from the config or contains no `.qmd` files. A minimal
starter report is in
[`example/reports/summary/summary.qmd`](example/reports/summary/summary.qmd).

## Outputs

Everything lands under `results_dir` (per-step folders in the table above);
per-rule logs under `logs_dir`. Key results:

- `011_deseq2/<contrast>.csv` — differential expression (DESeq2)
- `012_deseq2_gene_symbols/<contrast>.csv` — the same, with gene symbols
- `013_gsva/<contrast>.csv` and `gsva_scores.csv` — GSVA enrichment
- `014_gsea/<contrast>.csv` — GSEA (fgsea) enrichment
- `009_pca/` and `010_pca_samples_dropped/` — interactive PCA (HTML/PNG) + clustermaps
- `004_qc_reports/report_cards.txt` — per-sample text QC report cards
- `005_multiqc/multiqc_report.html` — aggregated QC
- `015_prose/methods.md` — auto-generated methods paragraph + references
- `reports/<name>/<name>.html` — rendered user-supplied Quarto reports

## Development

**TODO:**

- Unit tests
- Support single-end read samples
- Handle paired study design in DESeq2
- Add package parameters to prose (e.g., GSVA kcdf, tau)
- Autogenerate exploratory plots
- Put all Python code in the Snakemake rules into Python script files for easier unit testing
- Pipeline integration test
