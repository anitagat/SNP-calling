# SNP-calling
Pipeline for aligning samples to reference genomes, calling genomic SNPs and annotating them

## Part 1: Alignment 
### Part 1a: Alignment of single-end data
Assuming the existence of a folder containing sequencing data for all sequenced samples. 
In this experiments there are **2 replicates**.
Here we also assume that sequencing was **single-end**.

```
cd ~/data

# build genome index
cd ~/data/genome/

bowtie2-build wildtype.fna.fasta wildtype_index

# inspect sequences

cd ~/data/sequences/
head mutant_R1.fq
head mutant_R2.fq

wc -l mutant_R1.fq
wc -l mutant_R2.fq

# make folders
cd ~/data/
mkdir alignments
mkdir trimmed_sequences
mkdir fastp
cd fastp
mkdir mutant_strainX_R1
mkdir mutant_strainX_R2

## trim sequences for each replicate
cd ~/data/sequences

# trim replicate 1 and 
fastp -i mutant_strainX_R1.fq -o mutant_strainX_R1_trim.fq

mv ~/data/sequences/fastp.html ~/data/fastp/mutant_strainX_R1/fastp.html

mv ~/data/sequences/fastp.json ~/data/fastp/mutant_strainX_R1/fastp.json

mv ~/data/sequences/mutant_strainX_R1_trim.fq ~/data/trimmed_sequences/mutant_strainX_R1_trim.fq

# trim replicate 2
fastp -i mutant_strainX_R2.fq -o mutant_strainX_R2_trim.fq 
mv ~/data/sequences/fastp.html ~/data/fastp/mutant_strainX_R2/fastp.html

mv ~/data/sequences/fastp.json ~/data/fastp/mutant_strainX_R2/fastp.json

mv ~/data/sequences/mutant_strainX_R2_trim.fq ~/data/trimmed_sequences/mutant_strainX_R2_trim.fq

# align trimmed fq to genome index, outputs a sam file
# repeat for each replicate 

cd ~/data/trimmed_sequences/
bowtie2 -x ~/data/genome/wildtype_index -U mutant_strainX_R1_trim.fq -S mutant_strainX_R1_trim.sam
bowtie2 -x ~/data/genome/wildtype_index -U mutant_strainX_R2_trim.fq -S mutant_strainX_R2_trim.sam

mv ~/data/trimmed_sequences/mutant_strainX_R1_trim.sam ~/data/alignments/mutant_strainX_R1_trim.sam
mv ~/data/trimmed_sequences/mutant_strainX_R2_trim.sam ~/data/alignments/mutant_strainX_R2_trim.sam

# convert sam to bam, sort and index
cd ~/data/alignments
samtools view -bS mutant_strainX_R1_trim.sam > mutant_strainX_R1_trim.bam
samtools view -bS mutant_strainX_R2_trim.sam > mutant_strainX_R2_trim.bam

samtools sort mutant_strainX_R1_trim.bam -o mutant_strainX_R1_trim_sorted.bam
samtools sort mutant_strainX_R2_trim.bam -o mutant_strainX_R2_trim_sorted.bam

samtools index mutant_strainX_R1_trim_sorted.bam
samtools index mutant_strainX_R2_trim_sorted.bam

# tidy
rm mutant_strainX_R1_trim.sam
rm mutant_strainX_R2_trim.sam
rm mutant_strainX_R1_trim.bam
rm mutant_strainX_R2_trim.bam`

```

### Part 1b: Alignment of paired-end data 
The following code assumes the existence of a folder called paired.end where fq files containing the sequences are.

```
# make folders
cd ~/data/paired.end/
mkdir alignments
mkdir trimmed_sequences
mkdir fastp

# iterate through samples
for sample in sample_rep1 sample_rep2
do

  # run pipeline for sample (paired end)

  echo $sample
  # trim 
  fastp -i sequences/$sample"_R1.fq" -I sequences/$sample"_R2.fq" -o trimmed_sequences/$sample"_R1.trim.fq" -O trimmed_sequences/$sample"_R2.trim.fq"

  # align (paired end)
  bowtie2 -x genome/wildtype_index -1 trimmed_sequences/$sample"_R1.trim.fq" -2 trimmed_sequences/$sample"_R2.trim.fq" -S alignments/$sample".sam"

  # sam into bam 
  samtools view -bS alignments/$sample".sam" > alignments/$sample".bam"
  
  # sort
  samtools sort alignments/$sample".bam" -o alignments/$sample".sorted.bam"

  # index
  samtools index alignments/$sample".sorted.bam"
  
  # tidy up the sample folder
  mkdir fastp/$sample
  mv fastp.html fastp/$sample/fastp.html
  mv fastp.json fastp/$sample/fastp.json
  rm alignments/$sample".sam"
  rm alignments/$sample".bam" 

done

```

## Part 2: Call SNPs
Here, we use the **free bayes** tool to find SNPs. 

Freebayes outputs a VCF (variant call format) file. 

```
cd ~/data
mkdir vcf
freebayes -f genome.fasta sample_sorted.bam -C 5 > output.vcf
```
The -C 5 option is the minimum number of supporting observations to make the SNP call.

## Part 3: Annotate SNPs and predict effects
First, we run SnpEff, i.e. *a genetic variant annotation, and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes) [source: SnpEff](http://pcingola.github.io/SnpEff/).*

```
# find the reference genome of interest using grep from snpEff database
java -jar snpEff.jar databases | grep -i musculus

# run SnpEff
java -jar snpEff.jar GRCm38.99 sample.vcf -v > sample.annotated.vcf

# grep mutations that are predicted to have a high effect on the protein
# by specifying the flags -o -P, it is possible to only grep a defined range of characters around the word to be grepped from the file

grep -o -P '.{0,30}LOW.{0,100}' vcf/sample_rep1.sorted.annotated.vcf
grep -o -P '.{0,30}MODERATE.{0,100}' vcf/sample_rep1.sorted.annotated.vcf
grep -o -P '.{0,30}HIGH.{0,100}' vcf/sample_rep1.sorted.annotated.vcf

```
## Part 4: Complete the pipeline 
Here's a comprehensive and automated pipeline that aligns samples, calls SNPs and annotates them. 

```
# move to base directory
cd ~/data/

# variables (these are examples)
bowtie_index="wildtype_index"
snpEff_organism="Staphylococcus_aureus"

# make folders
mkdir alignments
mkdir trimmed_sequences
mkdir fastp
mkdir vcf
mkdir annotated_vcf

# iterate through samples
for sample in sample_rep1 sample_rep2
do

  # run pipeline for sample
  echo $sample
  
  # trim
  fastp -i sequences/$sample"_R1.fq" -I sequences/$sample"_R2.fq" -o trimmed_sequences/$sample"_R1.trim.fq" -O trimmed_sequences/$sample"_R2.trim.fq"
  # align
  bowtie2 -x genome/$bowtie_index -1 trimmed_sequences/$sample"_R1.trim.fq" -2 trimmed_sequences/$sample"_R2.trim.fq" -S alignments/$sample".sam"
  # sam to bam 
  samtools view -bS alignments/$sample".sam" > alignments/$sample".bam"
  samtools sort alignments/$sample".bam" -o alignments/$sample".sorted.bam"
  # index
  samtools index alignments/$sample".sorted.bam"
  # call SNPs 
  freebayes -f genome/wildtype.fna.fasta alignments/$sample".sorted.bam" -C 5 > vcf/$sample".sorted.vcf"
  # annotate SNPs
  java -jar /home/john/Documents/snpEff_latest_core/snpEff/snpEff.jar $snpEff_organism vcf/$sample".sorted.vcf" -v > annotated_vcf/s$sample".sorted.annotated.vcf"
  
  # tidy up the sample folders
  mkdir fastp/$sample
  mv fastp.html fastp/$sample/fastp.html
  mv fastp.json fastp/$sample/fastp.json
  mkdir annotated_vcf/$sample
  mv snpEff_genes.txt annotated_vcf/$sample/snpEff_genes.txt
  mv snpEff_summary.html annotated_vcf/$sample/snpEff_summary.html
  rm alignments/$sample".sam"
  rm alignments/$sample".bam" 

done

```