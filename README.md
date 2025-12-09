# VCFperfcomp

Compares the variants from a VCF file (e.g. produced by a bioinformatics pipeline) against an external "truth" VCF, computing accuracy metrics.

Designed to validate variant-calling pipelines that use tools like BWA-MEM + bcftools mpileup/bcftools call, but it works with any pair of VCFs that use consistent coordinates.

The script loads both VCFs, converts them to sets of (CHROM, POS, REF, ALT) tuples, and performs exact matching. It also supports restricting evaluation to specific regions using a BED file.


## HOW TO USE

### Quick install

Recommended: Set up a virtual machine and activate it

```bash
python -m venv .venv
source .venv/bin/activate
```

Install required python libraries
```bash
pip install -r requirements.txt
```


### Basic usage (required arguments)

```bash
python compare_vcf_performance.py \
    --truth true_variants.vcf \
    --pred pipeline_variants.vcf
```

### Restricting evaluation to high-confidence regions using a BED file

```bash
python compare_vcf_performance.py \
    --truth true_variants.vcf \
    --pred pipeline_variants.vcf \
    --bed highconf_regions.bed
```

### Saving global results to JSON and CSV

```bash
python compare_vcf_performance.py \
    --truth true_variants.vcf \
    --pred pipeline_variants.vcf \
    --save-json output.json \
    --save-csv output.csv
```

### Saving per-chromosome results

```bash
python compare_vcf_performance.py \
    --truth true_variants.vcf \
    --pred pipeline_variants.vcf \
    --per-chrom \
    --save-perchrom-csv output_per_chrom.csv
```

### Full example: truth + predictions + BED + all output files

```bash
python compare_vcf_performance.py \
    --truth true_variants.vcf \
    --pred pipeline_variants.vcf \
    --bed highconf_regions.bed \
    --save-json output.json \
    --save-csv output.csv \
    --per-chrom \
    --save-perchrom-csv output_per_chrom.csv
```


## OUTPUT FILES

```--save-json <file>```

Writes a JSON file with overall statistics:

```json
{
    "TP": ...,
    "FP": ...,
    "FN": ...,
    "precision": ...,
    "recall": ...,
    "f1": ...,
    "truth_variants": ...,
    "predicted_variants": ...
}
```

```--save-csv <file>```

Writes a CSV table with two columns:

| metric             | value |
|--------------------|-------|
| TP                 | ...   |
| FP                 | ...   |
| FN                 | ...   |
| precision          | ...   |
| recall             | ...   |
| f1                 | ...   |
| truth_variants     | ...   |
| predicted_variants | ...   |

```--save-perchrom-csv <file>```

Writes one line per chromosome (or contig):

| chrom | TP | FP | FN | precision | recall | f1 |
|-------|----|----|----|-----------|--------|----|
| chr1  | …  | …  | …  | …         | …      | …  |
| chr2  | …  | …  | …  | …         | …      | …  |
| …     | …  | …  | …  | …         | …      | …  |

Useful for detecting chromosomes with abnormal error patterns.


## NOTES


* The script uses exact matching of ```(CHROM, POS, REF, ALT)```.
* Skips ALT alleles that are not standard SNVs or INDELs.
* REF and ALT fields are automatically uppercased before comparison.
* BED files must follow standard BED format: 0-based start, 1-based end.
* VCFs must use the same reference genome version for meaningful comparisons.

