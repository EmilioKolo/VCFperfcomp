#!/usr/bin/env bash

set -euo pipefail




python compare_vcf_performance.py \
        --truth ./data/ecoli_breseq/SRR2983178/output/output.vcf \
        --pred ./data/ecoli/SRR2983178_variants.vcf.gz \
        --save-json ./data/comparisons/SRR2983178_variants.json \
        --save-csv ./data/comparisons/SRR2983178_variants.csv \
        --per-chrom \
        --save-perchrom-csv ./data/comparisons/SRR2983178_variants_per_chrom.csv


python compare_vcf_performance.py \
        --truth ./data/ecoli_breseq/SRR2983178/output/output.vcf \
        --pred ./data/ecoli/SRR2984607_snpeff.vcf.gz \
        --save-json ./data/comparisons/SRR2984607_snpeff.json \
        --save-csv ./data/comparisons/SRR2984607_snpeff.csv \
        --per-chrom \
        --save-perchrom-csv ./data/comparisons/SRR2984607_snpeff_per_chrom.csv
