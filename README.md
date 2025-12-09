# RNA-seq Analysis: Control vs Experimental Groups

## Overview

This repository contains a comprehensive RNA-seq analysis pipeline comparing control and experimental groups, with integration of literature data for validation and context. The analysis aims to identify differential gene expression patterns and provide functional insights into the biological processes affected by the experimental condition.

## Project Structure

```
RNAsq/
├── README.md                 # This file
├── data/
│   ├── raw/                 # Raw sequencing data
│   ├── processed/           # Quality-controlled and aligned reads
│   ├── literature/          # Reference data from literature
│   └── metadata/            # Sample metadata and annotations
├── scripts/
│   ├── qc/                  # Quality control scripts
│   ├── alignment/           # Read alignment scripts
│   ├── quantification/      # Gene expression quantification
│   ├── differential/        # Differential expression analysis
│   └── visualization/       # Plotting and visualization scripts
├── results/
│   ├── qc_reports/         # QC summary reports
│   ├── alignment_stats/     # Alignment statistics
│   ├── gene_counts/         # Raw and normalized gene counts
│   ├── deseq2/             # DESeq2 analysis results
│   ├── functional/         # GO/KEGG enrichment analyses
│   └── figures/            # Publication-quality figures
├── notebooks/               # Jupyter notebooks for exploratory analysis
├── docs/                   # Detailed documentation
├── environment.yml         # Conda environment specification
└── config/                 # Configuration files for analysis
```

## Analysis Workflow

### 1. Quality Control (QC)
- **FastQC**: Assessment of raw sequencing read quality
- **MultiQC**: Aggregate quality control reports across all samples
- **Adapter Detection**: Identification and documentation of sequencing adapters
- **Quality Metrics**:
  - Per-base sequence quality
  - Per-sequence GC content
  - Sequence length distribution
  - Duplication levels

### 2. Read Alignment
- **Reference Genome**: [Specify genome version]
- **Alignment Tool**: HISAT2/STAR for read alignment
- **Alignment Quality Metrics**:
  - Overall alignment rate (target: >85%)
  - Uniquely mapped reads percentage
  - Multi-mapped reads assessment
  - Strand specificity validation

### 3. Gene Expression Quantification
- **Quantification Method**: HTSeq-count or featureCounts
- **Annotation**: GENCODE/RefSeq gene annotations
- **Output**: Gene-level read counts matrix
- **Normalization**: TMM (Trimmed Mean of M-values) or RLE (Relative Log Expression)

### 4. Differential Expression Analysis
- **Statistical Framework**: DESeq2
- **Comparison**: Control vs Experimental groups
- **Filtering**:
  - Minimum count threshold: 10 reads (across all samples)
  - FDR adjustment: Benjamini-Hochberg (α = 0.05)
  - Log2 fold-change threshold: |log2FC| ≥ 1
- **Output Metrics**:
  - Log2 fold-change (log2FC)
  - P-values and adjusted p-values
  - Base mean expression levels
  - Effect sizes

### 5. Functional Enrichment Analysis
- **Gene Ontology (GO)**: Biological processes, molecular functions, cellular components
- **KEGG Pathways**: Metabolic and signaling pathway analysis
- **Tools**: clusterProfiler, enrichR, or g:Profiler
- **Filtering**: Adjusted p-value < 0.05, minimum gene set size

### 6. Literature Data Integration
- **Comparative Analysis**: Compare identified DEGs with published studies
- **Validation**: Cross-reference results with:
  - Published RNA-seq studies on similar conditions
  - Microarray datasets from GEO database
  - TCGA or other public repositories
- **Meta-analysis**: Assess consistency of findings across studies

## Sample Information

### Control Group
- **Description**: [Describe control samples]
- **Number of Replicates**: [N]
- **Sample IDs**: [List or reference]
- **Biological Conditions**: [Specify conditions]

### Experimental Group
- **Description**: [Describe experimental samples]
- **Number of Replicates**: [N]
- **Sample IDs**: [List or reference]
- **Biological Conditions**: [Specify conditions]

### Key Metadata
- Sequencing Platform: [e.g., Illumina NovaSeq, NextSeq]
- Read Type: Single-end / Paired-end
- Read Length: [e.g., 100 bp, 150 bp]
- Estimated Library Size: [M reads]
- Sequencing Depth: [M reads per sample]

## Key Results Summary

### Differential Expression Statistics
- **Total Genes Analyzed**: [Number]
- **Significantly Expressed Genes** (FDR < 0.05):
  - Upregulated in Experimental: [Number]
  - Downregulated in Experimental: [Number]
- **Validation Rate** (literature concordance): [Percentage]

### Top Differentially Expressed Genes

| Gene Symbol | Gene ID | log2FC | p-value | adj.p-value |
|------------|---------|--------|---------|-------------|
| [Gene1] | [ID1] | [FC] | [p] | [padj] |
| [Gene2] | [ID2] | [FC] | [p] | [padj] |
| [Gene3] | [ID3] | [FC] | [p] | [padj] |

*Full results available in `results/deseq2/`*

### Top Enriched Pathways (KEGG)

| Pathway ID | Pathway Name | Gene Count | p-value | Regulation |
|-----------|--------------|-----------|---------|-----------|
| [Path1] | [Name1] | [N] | [p] | Up/Down |
| [Path2] | [Name2] | [N] | [p] | Up/Down |

*Full enrichment results available in `results/functional/`*

### Top Gene Ontology Terms (GO)

| Category | Term ID | Term Name | Gene Count | p-value |
|----------|---------|-----------|-----------|---------|
| BP | [GO1] | [Term1] | [N] | [p] |
| MF | [GO2] | [Term2] | [N] | [p] |
| CC | [GO3] | [Term3] | [N] | [p] |

*BP = Biological Process, MF = Molecular Function, CC = Cellular Component*

## Validation Against Literature

### Published Studies - Consistency Analysis

1. **Study Citation [Year]**
   - Overlap: [X]% of genes in common
   - Direction Concordance: [Y]%
   - Key Overlapping Genes: [List]

2. **Study Citation [Year]**
   - Overlap: [X]% of genes in common
   - Direction Concordance: [Y]%
   - Key Overlapping Genes: [List]

### Known Markers Validation
- [Marker Gene 1]: [Expression Status] (Expected: [Expected], Found: [Found])
- [Marker Gene 2]: [Expression Status] (Expected: [Expected], Found: [Found])

## Methods

### Data Processing Details

#### Raw Data QC
```bash
fastqc *.fastq.gz
multiqc .
```

#### Alignment
```bash
hisat2 -x reference_index -1 control_R1.fastq.gz -2 control_R2.fastq.gz | samtools sort -o control.bam
samtools index control.bam
```

#### Quantification
```bash
featureCounts -p -T 8 -a genes.gtf -o counts.txt *.bam
```

#### Differential Expression
- Statistical test: Wald test (DESeq2)
- Dispersion estimation: Maximum likelihood estimation
- Multiple testing correction: Benjamini-Hochberg

### Statistical Parameters
- Significance threshold: FDR < 0.05
- Fold-change threshold: |log2FC| ≥ 1
- Minimum expression: >0 counts in ≥3 samples
- Normalization method: DESeq2 median of ratios

## Software and Dependencies

### Core Tools
- **FastQC** (v0.11.x): Quality control
- **MultiQC** (v1.12+): Aggregate QC reports
- **HISAT2** (v2.2.x) or **STAR** (v2.7.x): Read alignment
- **Samtools** (v1.13+): BAM file processing
- **featureCounts** (v2.0+): Read quantification
- **R** (v4.0+): Statistical analysis
- **DESeq2** (v1.32+): Differential expression
- **clusterProfiler** (v4.0+): Functional enrichment

### Python Packages
```
pandas >= 1.3.0
numpy >= 1.20.0
matplotlib >= 3.4.0
seaborn >= 0.11.0
scipy >= 1.7.0
scikit-learn >= 0.24.0
```

### Installation

**Using Conda:**
```bash
conda env create -f environment.yml
conda activate rnaseq_analysis
```

**Using Docker:**
```bash
docker pull [your-docker-image]
docker run -it [your-docker-image] /bin/bash
```

## Running the Analysis

### Step-by-step Execution

1. **Quality Control**
   ```bash
   bash scripts/qc/run_fastqc.sh
   bash scripts/qc/run_multiqc.sh
   ```

2. **Alignment**
   ```bash
   bash scripts/alignment/align_reads.sh
   bash scripts/alignment/bam_processing.sh
   ```

3. **Quantification**
   ```bash
   bash scripts/quantification/feature_count.sh
   ```

4. **Differential Expression**
   ```bash
   Rscript scripts/differential/deseq2_analysis.R
   ```

5. **Functional Enrichment**
   ```bash
   Rscript scripts/functional/enrichment_analysis.R
   ```

6. **Visualization**
   ```bash
   python scripts/visualization/generate_plots.py
   ```

### Using Snakemake (Optional)

```bash
snakemake --cores 16 --use-conda
```

## Output Files Description

### Quality Control
- `results/qc_reports/multiqc_report.html`: Interactive QC summary
- `results/qc_reports/fastqc_*.html`: Individual sample QC reports

### Alignment Statistics
- `results/alignment_stats/alignment_summary.txt`: Overall mapping statistics
- `results/alignment_stats/per_sample_metrics.csv`: Per-sample statistics

### Gene Counts
- `results/gene_counts/raw_counts.csv`: Raw read counts matrix
- `results/gene_counts/normalized_counts.csv`: TMM-normalized counts
- `results/gene_counts/tpm_values.csv`: Transcripts per million

### Differential Expression
- `results/deseq2/deseq2_results.csv`: Complete DESeq2 output
- `results/deseq2/significant_genes.csv`: Filtered significant genes (FDR < 0.05)
- `results/deseq2/upregulated_genes.csv`: Upregulated genes in experimental group
- `results/deseq2/downregulated_genes.csv`: Downregulated genes in experimental group

### Functional Analysis
- `results/functional/go_enrichment.csv`: Gene Ontology enrichment results
- `results/functional/kegg_enrichment.csv`: KEGG pathway enrichment results
- `results/functional/enrichment_plots.pdf`: Visualization of enriched terms

### Figures
- `results/figures/pca_plot.pdf`: Principal component analysis
- `results/figures/volcano_plot.pdf`: Volcano plot of DE genes
- `results/figures/heatmap_top_genes.pdf`: Heatmap of top DE genes
- `results/figures/expression_boxplots.pdf`: Box plots of key genes

## Data Visualization

### Principal Component Analysis (PCA)
PCA plot visualizing overall sample clustering and reproducibility between replicates. Should show clear separation between control and experimental groups.

### Volcano Plot
Visualization of log2 fold-change vs. -log10 adjusted p-values. Highlights significantly up- and downregulated genes.

### Heatmap
Hierarchical clustering heatmap of top differentially expressed genes, showing expression patterns across all samples.

### MA Plot
Log2 fold-change vs. mean log-scale expression, identifying genes with significant changes in expression.

## Quality Control Checklist

- [ ] All samples pass FastQC quality assessment
- [ ] Alignment rates > 85%
- [ ] Strand specificity confirmed
- [ ] Replicate samples cluster together in PCA
- [ ] No outlier samples detected
- [ ] Library sizes comparable across samples
- [ ] Batch effects assessed and documented
- [ ] GC bias evaluated
- [ ] Duplication levels acceptable

## Troubleshooting

### Low Alignment Rates
- Check reference genome version compatibility
- Verify adapter trimming effectiveness
- Assess sample quality (may require resequencing)

### Inconsistent Results Between Replicates
- Investigate technical variability (sequencing run, library prep)
- Check for batch effects
- Consider removing outlier samples
- Verify sample identity

### No Significant DE Genes Found
- Verify biological condition differences
- Check sample grouping assignments
- Assess sequencing depth adequacy
- Consider reducing fold-change or p-value thresholds

### High False Discovery Rate
- Increase sample size/replicates
- Improve experimental design
- Check for contamination or mix-ups
- Verify batch correction is applied

## Reproducibility

### Random Seed
All analyses use fixed random seeds for reproducibility:
- DESeq2 analysis: `set.seed(12345)`
- Python scripts: `np.random.seed(12345)`

### Session Information
Complete session information for R and Python environments is documented in:
- `results/session_info_r.txt`
- `results/session_info_python.txt`

### Data Availability
Raw sequencing data availability and accession numbers:
- [SRA Accession Numbers or Data Repository Links]
- [GEO Accession Numbers if applicable]

## Citation and References

### Key References

1. [Reference 1] - Methods for RNA-seq analysis
2. [Reference 2] - DESeq2 methodology
3. [Reference 3] - Related biological findings
4. [Reference 4] - Enrichment analysis approaches

### How to Cite This Work
If using this analysis pipeline or results, please cite:

```
Author(s). (2025). RNA-seq Analysis: Control vs Experimental Groups. 
GitHub Repository: https://github.com/rincowang-rgb/RNAsq
```

## Contributing

### Guidelines for Contributors
1. Create a feature branch for new analyses
2. Document any new scripts or methods
3. Update results summary and figures
4. Add appropriate citations for new reference data
5. Submit pull request with detailed description

### Reporting Issues
Please report bugs, questions, or suggestions through the GitHub Issues page.

## License

[Specify your license - e.g., MIT, GPL-3.0, CC-BY-4.0]

## Contact and Support

**Project Lead**: [Your Name]
**Email**: [Your Email]
**GitHub**: [@rincowang-rgb](https://github.com/rincowang-rgb)

For questions, suggestions, or collaboration opportunities, please open an issue or contact the maintainers.

## Changelog

### Version 1.0 (2025-12-09)
- Initial comprehensive documentation
- Complete workflow documentation
- Methods and validation framework
- Software requirements and installation instructions
- Quality control and troubleshooting guides

---

**Last Updated**: 2025-12-09
**Status**: Active Development
