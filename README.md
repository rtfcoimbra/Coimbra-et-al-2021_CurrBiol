# Code used in Coimbra *et al.* (2021) - Current Biology

Code used to analyse whole-genome sequencing data of giraffe in [Coimbra *et al.* (2021)](https://doi.org/10.1016/j.cub.2021.04.033).

### Citation

    Coimbra RTF, Winter S, Kumar V, Koepfli K-P, Gooley RM, Dobrynin P, Fennessy J, Janke A
    Whole-genome analysis of giraffe supports four distinct species
    Current Biology, 2021, in press. https://doi.org/10.1016/j.cub.2021.04.033

### Kordofan giraffe genome assembly and annotation

- `assemble_chromium_genome.sh`: *de novo* genome assembly with [Supernova](https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/welcome).
- `check_assembly_quality`: genome assembly quality assessment with [QUAST](http://quast.sourceforge.net/index.html), [BUSCO](https://busco.ezlab.org/), and [JupiterPlot](https://github.com/JustinChu/JupiterPlot/tree/1.0).
- `process_assembly.sh`: index assembly, annotate repeats with [RepeatMasker](http://www.repeatmasker.org/) and [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/), mask FASTA, extract scaffolds >= 1 Mb from masked FASTA, and create a BED file of non-repetitive regions in scaffolds >= 1 Mb.

### Quality control and read mapping

- `sra2fastq.sh`: retrieve data from the Sequence Read Archive.
- `remove_10x_barcodes.sh`: remove 10X Chromium barcodes from sample ENP11 with [proc10xG](https://github.com/ucdavis-bioinformatics/proc10xG).
- `trim_reads.sh`: trim paired-end reads with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
- `map_reads.sh`: map reads against a reference assembly with [BWA](https://github.com/lh3/bwa) and sort the output BAMs with [SAMtools](https://www.htslib.org/).
- `merge_bams.sh`: merge same sample BAMs from different lanes with [SAMtools](https://www.htslib.org/).
- `edit_read_id.sh`: remove 10X Chromium barcode sequences from read ID in the BAM file of sample ENP11.
- `mark_duplicates.sh`: mark PCR/optical duplicate reads with [Picard](https://broadinstitute.github.io/picard/).
- `check_mapping_quality.sh`: assess mapping quality with [QualiMap](http://qualimap.bioinfo.cipf.es/).
- `realigner_target_creator.sh`: create list of target intervals for indel realignment with [GATK](https://software.broadinstitute.org/gatk/).
- `indel_realigner.sh`: perform local realignment around indels with [GATK](https://software.broadinstitute.org/gatk/).
- `clean_bams.sh`: remove bad reads (flags 4, 256, 512, 1024 or 2048) from BAM files and keep only properly paired reads (flag 2) mapped to non-repetitive regions in scaffolds >= 1 Mb with [SAMtools](https://www.htslib.org/).
-

### SNP calling and linkage pruning

- `snp_calling.sh` and `site_depth_stats.py`: calculate the global site depth for multiple BAMs with [sambamba](https://github.com/biod/sambamba), compute summary statistics from the global site depth distribution (i.e. 5th and 95th percentiles, median, and median absolute deviation) with Python (NumPy and SciPy), and estimate genotype likelihoods with [ANGSD](https://github.com/ANGSD/angsd).
- `ld_pruning.sh`: calculate pairwise linkage disequilibrium with [ngsLD](https://github.com/fgvieira/ngsLD), plot LD decay curve with `fit_LDdecay.R`, prune linked sites with `prune_graph.pl`, and extract unlinked sites from ANGSD's genotype likelihoods file.

### Population structure and admixture analyses

- `pcangsd_hwe.sh`: estimate covariance matrix and perform a Hardy-Weinberg equilibrium test with [PCAngsd](https://github.com/Rosemeis/pcangsd).
- `plot_3dPCA.R`: perform PCA and make a 3-dimensional plot.
- `ngsadmix.sh`: estimate admixture proportions with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix).
- `plot_admixture.R`: plot admixture proportions with a custom R script and [ggplot2](https://ggplot2.tidyverse.org/index.html).
- `plot_likelihoods.R`: plot run likelihoods for each K.

### Genetic distances and neighbor-joining tree

### Nuclear phylogenomic inference

### Assembly and phylogeny of mitochondrial genomes

### Divergence time estimates

### Demographic reconstruction

### Heterozygosity

### ROH and Inbreeding
