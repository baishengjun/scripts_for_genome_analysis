# scripts_for_genome_analysis
## Genome assembly
```
canu -p ONT -d populus_deltoides -t 40 genomeSize=450m -nanopore-raw total_populus_deltiodes.fastq.gz
```
## Polish genome for three rounds with racon and one round with Medaka using raw Nanopore reads
```
racon -t 40 total_populus_deltiodes.fastq racon01.paf ont.contigs.fasta >round01.fasta
racon -t 40 total_populus_deltiodes.fastq racon02.paf round01.fasta > round02.fasta
racon -t 50 total_populus_deltiodes.fastq racon03.paf round02.fasta > racon03.fasta
medaka_consensus -i total_populus_deltiodes.fastq.gz -d round03.fasta -o medaka_polished -t 40 -m r941_min_fast_g303
```
## Polish genome using NGS data 
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
##  Identify and remove redundant heterozygous contigs
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
```
BuildDatabase -name populus_deltoides -engine ncbi populus_deltoides.fasta
RepeatModeler -engine ncbi -pa 30 -database populus_deltoides
queryRepeatDatabase.pl -species "embryophyta" > repeatmasker_embryophyta.fa
cat consensi.fa.classified repeatmasker_embryophyta.fa > merge.lib
RepeatMasker -e ncbi -pa 10 -dir . -gff -xsmall -nolow -lib merge.lib  populus_deltoides.fasta
```
## Gene prediction
### homolog-based
```
funannotate util prot2genome -g populus_deltoides.fasta.masked  -p merge_protein.fa -f tblastn  --cpus 30 -o homolog_deltoides.gff3 --logfile homlog.log
```
### RNA-seq-assisted using Trinity and PASA
```
Trinity --seqType fq --left totalRNAseq_1.fq --right totalRNAseq_2.fq --CPU 40 --max_memory 200G 
seqclean Trinity.fasta  -v UniVec
Launch_PASA_pipeline.pl -c alignAssembly.config  -C -R -g populus_deltoides.fasta.masked -t Trinity.fasta.clean -T -u 
pasa_asmbls_to_training_set.dbi --pasa_transcripts_fasta sql_lite.assemblies.fasta --pasa_transcripts_gff3 sql_lite.pasa_assemblies.gff3
pasa_asmbls_to_training_set.extract_reference_orfs.pl sql_lite.assemblies.fasta.transdecoder.genome.gff3 300 > best_candidates.gff3
```
### ab initio prediction
#### AUGUSTUS
```
gff2gbSmallDNA.pl best_candidates.gff3 
populus_deltoides.fasta 1000 genes.raw.gb
new_species.pl  --species=for_bad_genes_removing --AUGUSTUS_CONFIG_PATH=/mnt/sdf1/boshengjun/miniconda3/envs/augustus/config
etraining --species=for_bad_genes_removing --stopCodonExcludedFromCDS=false  genes.raw.gb 2> train.err
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
filterGenes.pl badgenes.lst genes.raw.gb > genes.db
randomSplit.pl genes.db 200
new_species.pl --species=populus_deltiodes_train
etraining --species=populus_deltiodes_train genes.db.train > train.out
augustus --species=populus_deltiodes_train genes.db.test | tee firsttest.out
optimize_augustus.pl --species=populus_deltiodes_train --cpus=8 --rounds=5 genes.db.train
augustus --species=deltiodes_train --gff3=on
populus_deltoides.fasta.masked > deltoides_augustus.gff3
```
#### SNAP
```
maker2zff best_candidates.gff3
perl chomp.pl genome.dna > deltiodes.dna
perl sequence_by_id.pl deltiodes.dna populus_deltoides.fasta.masked > genome.dna
fathom genome.ann genome.dna  -gene-stats &> gene-stats.log
fathom genome.ann genome.dna  -validate &> validate.log
fathom genome.ann genome.dna -categorize 1000 &
fathom uni.ann uni.dna -export 1000 -plus
mkdir params && cd params
forge ../export.ann ../export.dna
hmm-assembler.pl deltiodes params > deltiodes.hmm
snap -gff deltiodes.hmm populus_deltoides.fasta.masked  > deltoides_snap.gff
```
#### GeneMark-ES
```
gmes_petap.pl --ES --core 10 --sequence populus_deltoides.fasta.masked
gtf2gff.pl --gff3 <genemark.gtf --out=genemarkES_deltoides.gff3
```
### Combine results using EVM
```
partition_EVM_inputs.pl --genome populus_deltoides.fasta.masked  --gene_predictions denovo.gff3 --protein_alignments homolog_deltoides.gff3 --transcript_alignments sql_lite.assemblies.fasta.transdecoder.genome.gff3 --repeats repeat.gff3  --segmentSize 100000 --overlapSize 10000 --partition_listing partitions_list.out
write_EVM_commands.pl --genome populus_deltoides.fasta --repeats repeat.gff3 --gene_predictions denovo.gff3 --transcript_alignments sql_lite.assemblies.fasta.transdecoder.genome.gff3  --weights /mnt/sdf1/boshengjun/worka/populus_deltoides/Gene_prediction/abinitio/EVM/weights.txt --output_file_name evm.out --partitions partitions_list.out > commands.list
execute_EVM_commands.pl commands.list
convert_EVM_outputs_to_GFF3.pl --partitions partitions_list.out --output evm.out --genome populus_deltoides.fasta
find . -regex ".*evm.out.gff3" -exec cat {} \; > EVM.all.gff3
```
### PASA update 2 rounds
```
Load_Current_Gene_Annotations.dbi -c ../run_pasa/alignAssembly.config -g populus_deltoides.fasta.masked -P EVM.all.gff3 
Launch_PASA_pipeline.pl -c annotationCompare.config -A -g populus_deltoides.fasta.masked -t Trinity.fasta.clean
Load_Current_Gene_Annotations.dbi -c ../run_pasa/alignAssembly.config -g populus_deltoides.fasta.masked -P sql_lite.gene_structures_post_PASA_updates.50691.gff3 
Launch_PASA_pipeline.pl -c annotationCompare.config -A -g populus_deltoides.fasta.masked -t Trinity.fasta.clean
```



