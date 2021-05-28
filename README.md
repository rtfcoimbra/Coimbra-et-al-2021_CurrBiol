# Code from: Whole-genome analysis of giraffe supports four distinct species

Code used to analyze whole-genome sequencing data of giraffe in [Coimbra *et al.* (2021)](https://doi.org/10.1016/j.cub.2021.04.033):

- Coimbra RTF, Winter S, Kumar V, Koepfli K-P, Gooley RM, Dobrynin P, Fennessy J, Janke A (2021) Whole-genome analysis of giraffe supports four distinct species. *Current Biology*, in press. https://doi.org/10.1016/j.cub.2021.04.033

**Note:** The scripts used to generate plots do not necessarily reproduce the figures exactly as shown in the paper. In many cases, I used a free image editing software, namely [Krita](https://krita.org/en/), to assemble independent plots and add or correct plot annotations. I tried to use code as much as I could for reproducibility, but my skills are still not quite in the level of coding everything.

### Kordofan giraffe genome assembly and annotation

- `denovo_assembly.sh`: generate a *de novo* pseudohaploid genome assembly from a Chromium-prepared library using [Supernova](https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/welcome).
- `assembly_qc`: assess genome assembly quality with [QUAST](http://quast.sourceforge.net/index.html), [BUSCO](https://busco-archive.ezlab.org/v3/), and [JupiterPlot](https://github.com/JustinChu/JupiterPlot/tree/1.0).
- `process_assembly.sh`: index assembly, annotate repeats with [RepeatMasker](http://www.repeatmasker.org/) and [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/), mask FASTA, extract scaffolds >= 1 Mbp from masked FASTA, and create a BED file of non-repetitive regions in scaffolds >= 1 Mbp.
- `annotate_coding.sh`: soft-mask assembly and annotate protein-coding genes with [BRAKER](https://github.com/Gaius-Augustus/BRAKER) using protein sequences from *Bos taurus* ([GCA_002263795.2](https://www.ncbi.nlm.nih.gov/genome/?term=txid9913[orgn])) as extrinsic evidence.

### Quality control and read mapping

- `remove_10Xbarcodes.sh`: remove 10X Chromium barcodes from sample ENP11 with [proc10xG](https://github.com/ucdavis-bioinformatics/proc10xG).
- `trim_reads.sh`: trim paired-end reads with [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic).
- `map_reads.sh`: map reads against the Kordofan giraffe assembly with [BWA](https://github.com/lh3/bwa) and sort the output BAMs with [samtools](https://www.htslib.org/).
- `merge_bams.sh`: merge same sample BAMs from different lanes with [samtools](https://www.htslib.org/).
- `edit_readID.sh`: remove 10X Chromium barcode sequences from read ID in the BAM file of sample ENP11.
- `deduplicate_bams.sh`: mark PCR/optical duplicate reads with [Picard MarkDuplicates](https://broadinstitute.github.io/picard/).
- `mapping_qc.sh`: assess mapping quality with [QualiMap](http://qualimap.bioinfo.cipf.es/).
- `realigner_target_creator.sh`: create list of target intervals for indel realignment with [GATK](https://software.broadinstitute.org/gatk/).
- `indel_realigner.sh`: perform local realignment around indels with [GATK](https://software.broadinstitute.org/gatk/).
- `clean_bams.sh`: remove bad reads (flags 4, 256, 512, 1024 or 2048) from BAM files and keep only properly paired reads (flag 2) mapped to non-repetitive regions in scaffolds >= 1 Mbp with [samtools](https://www.htslib.org/).

### SNP calling and linkage pruning

- `snp_calling.sh` and `site_depth_stats.py`: calculate the global site depth for multiple giraffe BAMs with [sambamba](https://github.com/biod/sambamba), compute summary statistics from the global site depth distribution (i.e. 5th and 95th percentiles, median, and median absolute deviation) with Python (NumPy and SciPy), and estimate genotype likelihoods with [ANGSD](https://github.com/ANGSD/angsd).
- `ld_pruning.sh`: calculate pairwise linkage disequilibrium with [ngsLD](https://github.com/fgvieira/ngsLD), plot LD decay curve with [fit_LDdecay.R](https://github.com/fgvieira/ngsLD/blob/master/scripts/fit_LDdecay.R), prune linked sites with [prune_graph.pl](https://github.com/fgvieira/ngsLD/blob/master/scripts/prune_graph.pl), and extract unlinked sites from ANGSD's genotype likelihoods file.

### Population structure and admixture analyses

- `pcangsd_hwe.sh`: estimate covariance matrix and perform a Hardy-Weinberg equilibrium test with [PCAngsd](https://github.com/Rosemeis/pcangsd).
- `plot_3dPCA.R`: perform PCA and make a 3-dimensional plot with a custom R script.
- `ngsadmix.sh`: estimate admixture proportions with [NGSadmix](http://www.popgen.dk/software/index.php/NgsAdmix).
- `plot_admixture.R`: plot admixture proportions with a custom R script.
- `plot_likelihoods_deltaK.R`: plot run likelihoods for each K, generate input file for the delta K analysis ([Evanno *et al.* 2005](https://doi.org/10.1111/j.1365-294X.2005.02553.x)) implemented in the [CLUMPAK](http://clumpak.tau.ac.il/bestK.html) webserver, and plot delta K values (after running CLUMPAK separately).

### Genetic distances and neighbor-joining tree

- `snp_calling_gendist.sh` and `site_depth_stats.py`: similar to `snp_calling.sh` but including the okapi BAM as an outgroup.
- `gendist_tree.sh` and `make_pops_label.py`: calculate bootstrapped genetic distance matrices with [ngsDist](https://github.com/fgvieira/ngsDist) and generate a BioNJ tree with [FastME](http://www.atgc-montpellier.fr/fastme/).
- `plot_bionj_tree.R`: plot BioNJ tree with [ggtree](https://guangchuangyu.github.io/software/ggtree/).

### Nuclear phylogenomic inference

The scripts for phylogenomic analyses with genome fragments presented here were based on and modified from Fritjof Lammers' [original pipeline](https://github.com/mobilegenome/phylogenomics/tree/whales) used in [Árnason *et al.* (2018)](https://doi.org/10.1126/sciadv.aap9873).

- `generate_fasta.sh` and `site_depth_stats.py`: calculate per-base depth from cleaned BAMs using [sambamba](https://lomereiter.github.io/sambamba/) and generate consensus sequences in FASTA format with [ANGSD](https://github.com/ANGSD/angsd).
- `clean_fasta.sh` and `concatenate_fasta.py`: remove repetitive regions from the consensus sequences with [bedtools](https://bedtools.readthedocs.io/en/latest/index.html) and a custom Python script and calculate the percentage of ambiguous sites in each sequence.
- `gf_phylogenomics.txt`: workflow for the phylogenomic analyses with genome fragments (GFs) - i.e. generate GFs with a custom Python script, perform approximately unbiased (AU) tree topology test with [IQ-TREE](http://www.iqtree.org/), infer multispecies coalescent tree with [ASTRAL](https://github.com/smirarab/ASTRAL), analyze quartet frequencies with [DiscoVista](https://github.com/esayyari/DiscoVista).
- `genome_fragments.py`: generate GF alignments from consensus genome sequences.
- `check_gfs.py`: check the percentage of ambiguous sites (Ns) per sequence in each GF alignment and save a list of "good" and "bad" GFs (check code comments for description). Only "good" GFs were used for analyses.
- `iqtree_parallel.sh` and `iqtree_au_test.sh`: perform AU tree topology test with [IQ-TREE](http://www.iqtree.org/) and collect the inferred trees and p-AU values.
- `parse_iqlog.py`: extract *p*<sub>AU</sub> values from IQ-TREE's log file.
- `plot_pAU.R`: plot the *p*<sub>AU</sub> values of each topology for each GF size (with confidence intervals) and plot the alternative topologies analyzed in the AU test.

### Assembly and phylogeny of mitochondrial genomes

- `mitogenomes_workflow.txt`: subsample raw reads with [seqtk](https://github.com/lh3/seqtk), assemble and annotate mitochondrial genomes with [MitoZ](https://github.com/linzhi2013/MitoZ/tree/master/version_2.3), extract protein coding sequences, and construct maximum likelihood phylogeny based on protein coding sequences with [IQ-TREE](http://www.iqtree.org/). **Note:** some steps must be performed manually and/or with the use of independent GUI programs.
- `plot_mtdna_tree.R`: plot maximum likelihood mitochondrial tree with [ggtree](https://guangchuangyu.github.io/software/ggtree/).

### Demographic reconstruction

- `psmc_bootstrap.sh` and `site_depth_stats.py`: calculate per-base depth from cleaned BAMs using [sambamba](https://lomereiter.github.io/sambamba/), generate consensus genome sequences in FASTQ format with [bcftools](https://samtools.github.io/bcftools/), and infer giraffe effective population size changes through time with the [PSMC](https://github.com/lh3/psmc) model.

### Heterozygosity

- `generate_foldedSFS.sh`: generate the per sample folded site frequency spectrum (SFS) with [ANGSD](https://github.com/ANGSD/angsd).
- `plot_heterozygosity.R`: calculate and plot global heterozygosity per sample in R with [tidyverse](https://www.tidyverse.org/).

### ROH and Inbreeding

- `genotype_calling`: call genotypes with [bcftools](https://samtools.github.io/bcftools/) using the SNPs called with ANGSD and convert BCF to Oxford GEN format.
- `estimate_roh`: identify runs of homozygosity (ROH) and estimate realized inbreeding coefficients (F<sub>ROH</sub>) with [RZooRoH](https://doi.org/10.1111/2041-210X.13167).
- `plot_froh.R`: plot F<sub>ROH</sub> and the number vs. the accumulated length of ROH with [tidyverse](https://www.tidyverse.org/) and [patchwork](https://patchwork.data-imaginist.com/).
