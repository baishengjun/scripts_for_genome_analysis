# scripts_for_genome_analysis
## Genome assembly
```
canu -p ONT -d populus_deltoides -t 40 genomeSize=450m -nanopore-raw total_populus_deltiodes.fastq.gz
```
## polish genome for three rounds with racon and one round with Medaka using raw Nanopore reads
```
racon -t 40 total_populus_deltiodes.fastq racon01.paf ont.contigs.fasta >round01.fasta
racon -t 40 total_populus_deltiodes.fastq racon02.paf round01.fasta > round02.fasta
racon -t 50 total_populus_deltiodes.fastq racon03.paf round02.fasta > racon03.fasta
medaka_consensus -i total_populus_deltiodes.fastq.gz -d round03.fasta -o medaka_polished -t 40 -m r941_min_fast_g303
```
## polish genome using NGS data 
`nextpolish run.cfg`
### run.cfg
```
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 4
multithread_jobs = 10
genome =./consensus.fasta
genome_size = auto
workdir = ./nextpolished
polish_options = -p {multithread_jobs}

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100
```
##  identify and remove redundant heterozygous contigs
```
minimap2 -t 40 -ax map-ont nextpolished.fasta total_populus_deltiodes.fastq.gz --secondary=no  | samtools sort -m 5G -o aligned.bam -T tmp.ali
purge_haplotigs hist -b aligned.bam -g nextpolished.fasta -t 40
purge_haplotigs cov -i aligned.bam.gencov -l 10 -m 72 -h 180 
purge_haplotigs purge -g nextpolished.fasta -a 95 -c coverage_stats.csv -t 20 -d -b aligned.bam
```
## Scaffold contigs based on HiC data
### align the Hi-C reads to draft genome and computed the Hi-C contact frequency using Juicer pipeline
```
bwa index curated.fasta
python ~/juicer/misc/generate_site_positions.py MboI genome curated.fasta
awk 'BEGIN{OFS="\t"}{print $1, tre $NF}' genome_MboI.txt > MZHY.chrom.size
~/juicer/scripts/juicer.sh -g genome -s MboI -z reference/curated.fasta -p reference/MZHY.chrom.size -y reference/genome_MboI.txt -D ~/juicer -t 40 
```
### correct mis-joins, order, orient and anchor contigs into a candidate chromosome-length assembly using the 3D-DNA pipeline
```
~/bash-4.0/bash ~/3d-dna/run-asm-pipeline.sh -r 2 ../reference/curated.fasta ../aligned/merged_nodups.txt
```
### contig misassemblies and scaffold mis-joins were detected and corrected using JBAT softwarte, and then generate the curated_final.review.assembly file
### generate final assembly file after JBAT review
```
~/bash-4.0/bash ~/3d-dna/run-asm-pipeline-post-review.sh -r curated_final.review.assembly  ../../reference/curated.fasta ../../aligned/merged_nodups.txt
```
## BUSCO analysis
```
busco -m genome -i populus_deltoides.fasta -o run_busco -l embryophyta_odb10 -c 20
```
## Repeat annotation
## Noncoding RNA anntation
## Gene prediction
## Gene annatation
