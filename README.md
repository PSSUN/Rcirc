# Rcirc
Rcirc helps to identify the coding ability of circRNA and visualize the feature at various aspects. Rcirc provides general visualization for both single circRNA and meta-feature of thousands of circRNAs, Rcirc was designed as a user-friendly tool which covers lots of high automatically functions without running many complicated processes by user. 

## Introduction
The main goal of our study was to develop a user-friendly tool which covers the main demands for circRNA research. Rcirc is a capable and user-friendly package based on R language. The package provides numerous analyses for both upstream and downstream research include circRNA detection, coding ability identification, single feature analyses and visualization of meta-feature. Furthermore, the users can visualize the reads mapping situation for each junction site of circRNA by Rcirc with sequencing data. With growing attention on circRNA, Rcirc will become an auxiliary tool to encourage researchers to proceed with further analyses in circRNA world. All the detail of usage is included by the documents of Rcirc in GitHub's online pages.

## Documents
https://rcirc-doc.readthedocs.io/en/latest/index.html

## Usage overview

The functions under [`R/`](R) expose the main workflows in the package.  After
loading the project in R (for example with
`devtools::load_all("/path/to/Rcirc")`), you can call the exported helpers
directly.

### Detecting circRNAs

* `PredictCirc(sam, fa, out, gtf)` wraps the bundled **CIRI2** script and takes
  a mapped SAM file together with the reference genome FASTA and annotation GTF
  to produce a circRNA prediction report.【F:R/PredictCirc.R†L5-L16】
* `TranslateCirc(out, fastq, adapter, trimmomatic, genome, tmp, circgenome, rRNA)`
  orchestrates the preprocessing of Ribo-seq reads: it builds Bowtie/STAR
  indices, trims adapters with Trimmomatic, removes rRNA and linear reads, and
  finally aligns the cleaned reads to the circRNA genome.【F:R/TranslateCirc.R†L1-L73】

### Summarising junction support

* `circ_summary(bamfile, gff)` reads Ribo-seq alignments together with a circRNA
  annotation GFF, filters junction-spanning reads, and returns a list with the
  trimmed annotations, supporting reads, coverage summaries, and per-position
  views that drive the downstream plots.【F:R/rCirc.R†L23-L112】

### Visualising results

Rcirc contains multiple high-level plotting functions to explore the detected
circRNAs:

* `showOverview(circbed, gff, genomefasta, ribo, rna)` classifies circRNAs by
  genomic context (intergenic, exon, intron) and combines annotation, Ribo-seq,
  and RNA-seq summaries into a whole-genome overview plot.【F:R/showOverview.R†L1-L89】
* `showDistribution(gbed, cbed)` produces category counts for circRNAs based on
  their annotated start and end regions, helping compare intergenic, exonic, and
  intronic origins.【F:R/showDistribution.R†L1-L83】
* `showJunction(data, max = 20, title = 0)` inspects circRNA FASTA files and
  plots the most frequent splice signal combinations observed at back-splice
  junctions.【F:R/showJunction.R†L1-L37】
* `showMapping(summary, circ_index, genomefile, ...)` (exported as
  `showMapping`) visualises read density and sequence context around a selected
  circRNA junction using the summary list returned by `circ_summary`.【F:R/mappingPlot.R†L1-L107】
* Additional helpers such as `showCodon`, `showLength`, `classByType`, and
  `stemRing` focus on codon usage, length distribution, category summaries, and
  stem-loop structures respectively; see the source files in [`R/`](R) for their
  parameters and plotting options.

These functions can be combined in scripts or notebooks to build a complete
analysis workflow: predict circRNAs, process Ribo-seq support, summarise junction
evidence, and visualise the resulting catalogue.
