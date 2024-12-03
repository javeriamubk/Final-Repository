# Introduction
The goal of this workflow is to analyze a gene family through alignment, phylogenetic tree construction, domain prediction, and expression analysis. This pipeline provides a comprehensive analysis, including midpoint rooting and protein domain visualization.

There are:

Scripts for each step of the workflow.
Outputs for reference and reuse.
Documentation for reproducibility.
Workflow Structure
The repository is organized as follows:

```
├── 01_Raw_Sequences

├── 02_Alignment

├── 03_Phylogenetic_Analysis

├── 04_Domain_Prediction

├── 05_Expression_Analysis

├── scripts

```

### 1. Raw Data Preparation

Goal: Organize and prepare raw sequence data.


##### Steps:
Clone the repository:
```
Copy code: git clone https://github.com/username/GeneFamilyAnalysis.git
cd GeneFamilyAnalysis/01_Raw_Sequences
```
Copy your raw FASTA sequences to the 01_Raw_Sequences directory.
Scripts:

prepare_sequences.sh: Remove stop codons and preprocess sequences.
example_raw.fasta: Example input.

Output:
processed_sequences.fasta: Preprocessed sequences.

### 2. Quality Control and Alignment
Goal: Align sequences using MUSCLE for further analysis.

Steps:
Run the alignment:

```
Copy code: muscle -align processed_sequences.fasta -output aligned_sequences.aln
Check alignment quality using alv or other tools.
```
MUSCLE aligns the cleaned sequences (processed_sequences.fasta) and produces an alignment file (aligned_sequences.aln), which is important for making phylogenetic trees.




Scripts:
run_alignment.sh: Automates the alignment step.

Output:
aligned_sequences.aln: Alignment file for tree construction.

### 3. Phylogenetic Tree Construction

Goal: Build a phylogenetic tree and visualize it.

Steps:
Build the tree using IQ-TREE:

```
Copy code: iqtree -s aligned_sequences.aln -bb 1000 -nt 2
Midpoint rooting:
```
```
Copy code: gotree reroot midpoint -i treefile -o rooted_treefile
```
IQ-TREE makes a phylogenetic tree from the aligned sequences with bootstrap analysis (-bb 1000). "gotree" re-roots the tree at its midpoint to provide a more meaningful evolutionary perspective.


Visualization: Use plotUnrooted.R for unrooted trees:

```
Copy code: Rscript plotUnrooted.R treefile tree_unrooted.pdf
Scripts:
```
run_tree.sh: Automates tree construction and rooting.
Output:

treefile: Unrooted tree.

rooted_treefile: Rooted tree.

tree_unrooted.pdf: Visual representation.

### 4. Protein Domain Prediction

Goal: Identify protein domains using RPS-BLAST.

Steps:

Run RPS-BLAST with Pfam:

```
Copy code: rpsblast -query processed_sequences.fasta -db Pfam -out domains.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue 1e-10
Visualize domains on the phylogenetic tree:
```
RPS-BLAST identifies protein domains in the sequences using the Pfam database. The R script plots these domains next to the phylogenetic tree, creating visualizations of their distribution across evolutionary branches.

This command removes stop codons (*) from the input FASTA file and saves the cleaned sequences to a new file (processed_sequences.fasta). This step ensures compatibility with downstream tools that may not handle stop codons properly.

```
Copy code: Rscript plotTreeAndDomains.r rooted_treefile domains.rps-blast.out tree_with_domains.pdf
```
Scripts: run_rpsblast.sh: Automates RPS-BLAST.
plotTreeAndDomains.r: R script for visualization.

Output: domains.rps-blast.out: 

Domain prediction results.

tree_with_domains.pdf: Tree annotated with domains.

### 5. Gene Expression Analysis

Goal: Quantify and analyze differential expression.

Steps:
Index the reference transcriptome:

```
Copy code: kallisto index -i reference.index processed_sequences.fasta
```
Quantify gene expression:
```
Copy code: kallisto quant -i reference.index -o output_sample sample_R1.fastq sample_R2.fastq
```
Combine counts for analysis:

Use R scripts for DESeq2 analysis.

Scripts:
run_kallisto.sh: Automates quantification.
analyze_expression.R: Differential expression in R.

Output:

abundance.tsv: Counts and TPM values.
differential_expression_results.tsv: DE analysis output.

### 6. Writing Results and Figures

Goal: Prepare figures and summarize results.

Steps:
Generate plots and summary statistics:
R

### Example in R

library(ggplot2)

```
ggplot(data, aes(x = gene, y = expression)) + geom_bar()
```
This R code creates a bar plot of gene expression levels using ggplot2, producing a visualization of differential expression results. theme_minimal() applies a clean, publication-ready styling.


Combine results into a report:

Use markdown or LaTeX for formatting.

Scripts:

generate_plots.R: Create figures from DE results.

Output:
results_report.md: Summary of results.
figures: Plots and figures for publication.


### References

```
Pfam Database: Pfam
IQ-TREE: IQ-TREE Website
RPS-BLAST: NCBI RPS-BLAST
```
