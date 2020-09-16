# BachelorThesis
A metagenomic workflow.

For this documentation the pyr_d60 samples are used as examples for easier understanding of the command lines.

Starting the pipeline:

1. Quality check with FastQC
Using the Illumina Sequences obtained from environmental samples a quality check with FastQC is to be conducted. This analysis will show where the sequences need to be refined, such as trimming out adapters and potentially contaminated parts (mostly located at the beginning and end of the sequences). 

fastqc -o OUTPUT_PATH -f fastq PATH-TO-RAW-READS/pyr_d60_all_1.fq PATH-TO-RAW-READS/pyr_d60_all_2.fq -t num_threads

2. Trimming with BBDuk
After analysing the reads with FastQC, determine where trimming is necessary. Trimming is performed in 3 steps: Adapter trimming left, adapter trimming right and quality trimming. Make sure to always update the input file with the most recent file (i.e. using the righttrimmed file to perform the adapter trim on the left side of the sequence and then the rl trimmed file to perform the quality trim on).

  2.1 Right trim:

bbduk.sh t=num_threads in1=PATH-TO-RAW-READS/pyr_d60_all_1.fq in2=PATH-TO-RAW-READS/pyr_d60_all_2.fq out1=OUTPUT-PATH/pyr_d60_all_1_rtrim.fq out2=OUTPUT-PATH/pyr_d60_all_2_rtrim.fq ref=PATH-TO-ADAPTER-FILE/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo  

  2.2 Left trim:

bbduk.sh t=num_threads in1=PATH-TO-RIGHT-TRIMMED-FILE/pyr_d60_all_1_rtrim.fq in2=PATH-TO-RIGHTTRIMMED-FILE/pyr_d60_all_2_rtrim.fq out1=OUTPUT-PATH/pyr_d60_all_1_rltrimmed.fq out2=OUTPUT-PATH/pyr_d60_all_2_rltrimmed.fq ref=PATH-TO-ADAPTER-FILE/adapters.fa ktrim=l k=23 mink=11 hdist=1 tpe tbo  

  2.3 Quality trim:

The parameters ftl, ftr, trimq, maq and minlen have to be set according to the FastQC results.
bbduk.sh t=num_threads in1=PATH-TO-RL-TRIMMED-FILE/pyr_d60_all_1_clean.fq in2=PATH-TO-RL-TRIMMED-FILE/pyr_d60_all_2_clean.fq out1=OUTPUT-PATH-FOR-CLEAN-READS/pyr_d60_all_1_crisp.fq out2=OUTPUT-PATH-FOR-CLEAN-READS/pyr_d60_all_2_crisp.fq qtrim=rl trimq=20 ftl=6 ftr=144 maq=20 minlen=100

3. Repeat quality check with FastQC

4. Abundance Estimation with Kraken2 and Bracken
Abundance estimation gives a first overview of species identified in the sample. It is useful to check abundance before continuing the pipeline in case the target organism is not very abundant. The higher the abundance of an organism the higher the chance of forming high quality bins. When installing the packages for Kraken2 and Bracken check the packages in the environment using conda list (if using anaconda). Bracken can install kraken1 which will cause problems. If listed, remove before continuing. 

  4.1 Kraken2
Kraken2 is the base for Bracken to run on. This takes a while, go read some nice papers and have a coffee.

kraken2 --db PATH-TO-KRAKEN-DATABASE/kraken2_db --paired --classified-out pyr_d60#.fq PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq --threads num_threads --output PATH-TO-KRAKEN-OUTPUT/pyr_d60_Kraken.out --report Pyr_d60.report --confidence 0.05

  4.2 Bracken
Now that the initial work is done (thx Kraken<3) the estimation can be conducted. Bracken will default to estimate organisms on species level. Make sure to set the level (-l) for each iteration (levels: D=domain, P=phylum, C=class, O=order, F=family, G=genus, S=species (default)). The parameter -r sets the readlength (use what you set as minlen in the quality trim).

bracken -d PATH-TO-KRAKEN-DATABASE/kraken2_db -i PATH-TO-KRAKEN-REPORT-FILE/Pyr_d60_.report -o PATH-TO-BRACKEN-OUTPUT/Pyr_d60.bracken -r 100 -l D -t num_threads
bracken -d PATH-TO-KRAKEN-DATABASE/kraken2_db -i PATH-TO-KRAKEN-REPORT-FILE/Pyr_d60_.report -o PATH-TO-BRACKEN-OUTPUT/Pyr_d60.bracken -r 100 -l P -t num_threads
...

5. <i>De novo</i> Assembly using Megahit
The next step to take with the clean reads is to assemble them. This step tries to align sequences to lengthen reads. The goal is to achieve a sequence of a (nearly) complete genome/gene. This also runs about forever depending on you num_threads so go make some food at this point? -o creates a new directory. Set -m (usable memory) to go easy on the server and enable everyone else to also work.

megahit -1 PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq -2 PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq -o PATH-TO-ASSEMBLY-OUTPUT/pyr_d60_all.megahit_asm -m 0.2 -t num_threads

6. Readmapping using Bowtie2
Now that the reads have been assembled they have to be mapped for quantification. There are plenty mapping tools, in this step Bowtie2 was used.

  6.1 Building a Bowtie2 Index
bowtie2-build --threads num_threads PATH-TO-ASSEMBLED-CONTIGS/pyr_d60_all.megahit_asm/final.contigs.fa OUTPUT-PATH-FOR-MAPPED-CONTIGS/contigs

  6.2 Running the mapping
bowtie2 --threads num_threads -x PATH-TO-NEWLY-CREATED-INDEX-DIRECTORY/contigsindex -1 PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq -2 PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq -S OUTPUT-PATH-FOR-MAPPED-FILE/pyr_d60_all_mapped.sam

7. .sam > .bam
Various further needed programs require sorted and indexed bam files to work. Sorting and indexing can be done using only samtools (samtools view, samtools sort and samtools index) or using a combination of samtools and anvi'o. I use the combination method.

  7.1 creating the "raw" (not sorted or indexed) .bam file.

samtools view --threads num_threads -F 4 -bS PATH-TO-MAPPED-FILE/pyr_d60_all_reform_mapped.sam > OUTPUT-PATH/pyr_d60_all_reform_mapped_RAW.bam

7.2 sorting and indexing the file (r_bam > .bam and .bam.bai)
Since anvi'o and samtools require different versions of Python, switch to the environment that has the right Python version (3.0 I think) installed (The installed version of Python can be checked with <i> "conda list" </i> or <i> "conda version python" </i> in the desired environment). 

