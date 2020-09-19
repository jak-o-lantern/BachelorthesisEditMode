# A metagenomic workflow.

For this documentation the pyr_d60 samples are used as examples for easier understanding of the command lines.

<center> Starting the pipeline: </center>

<h3> <center> 1. Quality check with FastQC </center> </h3> 
Using the Illumina Sequences obtained from environmental samples a quality check with FastQC is to be conducted. This analysis will show where the sequences need to be refined, such as trimming out adapters and potentially contaminated parts (mostly located at the beginning and end of the sequences). 

<code> fastqc -o OUTPUT_PATH -f fastq PATH-TO-RAW-READS/pyr_d60_all_1.fq PATH-TO-RAW-READS/pyr_d60_all_2.fq -t num_threads </code>

<h3> <center> 2. Trimming with BBDuk </center> </h3> 
After analysing the reads with FastQC, determine where trimming is necessary. Trimming is performed in 3 steps: Adapter trimming left, adapter trimming right and quality trimming. Make sure to always update the input file with the most recent file (i.e. using the righttrimmed file to perform the adapter trim on the left side of the sequence and then the rl trimmed file to perform the quality trim on).

<h3> <center> 2.1 Right trim: </center> </h3>

<code> bbduk.sh t=num_threads in1=PATH-TO-RAW-READS/pyr_d60_all_1.fq in2=PATH-TO-RAW-READS/pyr_d60_all_2.fq out1=OUTPUT-PATH/pyr_d60_all_1_rtrim.fq out2=OUTPUT-PATH/pyr_d60_all_2_rtrim.fq ref=PATH-TO-ADAPTER-FILE/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo  </code>

<h3> <center> 2.2 Left trim: </center> </h3>

<code> bbduk.sh t=num_threads in1=PATH-TO-RIGHT-TRIMMED-FILE/pyr_d60_all_1_rtrim.fq in2=PATH-TO-RIGHTTRIMMED-FILE/pyr_d60_all_2_rtrim.fq out1=OUTPUT-PATH/pyr_d60_all_1_rltrimmed.fq out2=OUTPUT-PATH/pyr_d60_all_2_rltrimmed.fq ref=PATH-TO-ADAPTER-FILE/adapters.fa ktrim=l k=23 mink=11 hdist=1 tpe tbo </code>  

<h3> <center> 2.3 Quality trim: </center> </h3>
The parameters ftl, ftr, trimq, maq and minlen have to be set according to the FastQC results.

<code> bbduk.sh t=num_threads in1=PATH-TO-RL-TRIMMED-FILE/pyr_d60_all_1_clean.fq in2=PATH-TO-RL-TRIMMED-FILE/pyr_d60_all_2_clean.fq out1=OUTPUT-PATH-FOR-CLEAN-READS/pyr_d60_all_1_crisp.fq out2=OUTPUT-PATH-FOR-CLEAN-READS/pyr_d60_all_2_crisp.fq qtrim=rl trimq=20 ftl=6 ftr=144 maq=20 minlen=100 </code>

<h3> <center> 3. Repeat quality check with FastQC </center> </h3>
You know how to do this

<h3> <center> 4. Abundance Estimation with Kraken2 and Bracken </center> </h3>
Abundance estimation gives a first overview of species identified in the sample. It is useful to check abundance before continuing the pipeline in case the target organism is not very abundant. The higher the abundance of an organism the higher the chance of forming high quality bins. When installing the packages for Kraken2 and Bracken check the packages in the environment using conda list (if using anaconda). Bracken can install kraken1 which will cause problems. If listed, remove before continuing. 

<h3> <center> 4.1 Kraken2 </center> </h3>
Kraken2 is the base for Bracken to run on. This takes a while, go read some nice papers and have a coffee.

<code> kraken2 --db PATH-TO-KRAKEN-DATABASE/kraken2_db --paired --classified-out pyr_d60#.fq PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq --threads num_threads --output PATH-TO-KRAKEN-OUTPUT/pyr_d60_Kraken.out --report Pyr_d60.report --confidence 0.05 </code>

<h3> <center> 4.2 Bracken </center> </h3>
Now that the initial work is done (thx Kraken<3) the estimation can be conducted. Bracken will default to estimate organisms on species level. Make sure to set the level (-l) for each iteration (levels: D=domain, P=phylum, C=class, O=order, F=family, G=genus, S=species (default)). The parameter -r sets the readlength (use what you set as minlen in the quality trim). <br>

<code> bracken -d PATH-TO-KRAKEN-DATABASE/kraken2_db -i PATH-TO-KRAKEN-REPORT-FILE/Pyr_d60_.report -o PATH-TO-BRACKEN-OUTPUT/Pyr_d60.bracken -r 100 -l D -t num_threads </code>
<code> bracken -d PATH-TO-KRAKEN-DATABASE/kraken2_db -i PATH-TO-KRAKEN-REPORT-FILE/Pyr_d60_.report -o PATH-TO-BRACKEN-OUTPUT/Pyr_d60.bracken -r 100 -l P -t num_threads </code>
...

<h3> <center> 5. <i>De novo</i> Assembly using Megahit </center> </h3>
The next step to take with the clean reads is to assemble them. This step tries to align sequences to lengthen reads. The goal is to achieve a sequence of a (nearly) complete genome/gene. This also runs about forever depending on you num_threads so go make some food at this point? -o creates a new directory. Set -m (usable memory) to go easy on the server and enable everyone else to also work.

<code> megahit -1 PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq -2 PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq -o PATH-TO-ASSEMBLY-OUTPUT/pyr_d60_all.megahit_asm -m 0.2 -t num_threads </code>

<h3> <center> 6. Readmapping using Bowtie2 </center> </h3>
Now that the reads have been assembled they have to be mapped for quantification. There are plenty mapping tools, in this step Bowtie2 was used.

<h3> <center> 6.1 Building a Bowtie2 Index </center> </h3>

<code> bowtie2-build --threads num_threads PATH-TO-ASSEMBLED-CONTIGS/pyr_d60_all.megahit_asm/final.contigs.fa OUTPUT-PATH-FOR-MAPPED-CONTIGS/contigs </code>

<h3> <center> 6.2 Running the mapping </center> </h3>

<code> bowtie2 --threads num_threads -x PATH-TO-NEWLY-CREATED-INDEX-DIRECTORY/contigsindex -1 PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq -2 PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq -S OUTPUT-PATH-FOR-MAPPED-FILE/pyr_d60_all_mapped.sam </code>

<h3> <center> 7. .sam > .bam </center> </h3>
Various further needed programs require sorted and indexed bam files to work. Sorting and indexing can be done using only samtools (samtools view, samtools sort and samtools index) or using a combination of samtools and anvi'o. I use the combination method.

<h3> <center> 7.1 creating the "raw" (not sorted or indexed) .bam file. </center> </h3>

<code> samtools view --threads num_threads -F 4 -bS PATH-TO-MAPPED-FILE/pyr_d60_all_reform_mapped.sam > OUTPUT-PATH/pyr_d60_all_reform_mapped_RAW.bam </code>

<h3> <center> 7.2 sorting and indexing the file (r_bam > .bam and .bam.bai) </center> </h3>
Since anvi'o and samtools require different versions of Python, switch to the environment that has the right Python version (3.0 I think) installed (The installed version of Python can be checked with <code> "conda list" </code> or <code> "conda version python" </code> in the desired environment). <br>

<code> anvi-init-bam -T num_threads PATH-TO-MAPPED-RAW-BAM-FILE/pyr_d60_all_mapped_RAW.bam -o OUTPUT-PATH/pyr_d60_all_mapped.bam </code>

<h3> <center> 8. Binning </center> </h3>
This step will sort the assembled and hopefully nice and long contigs into bins. This means it will group those contigs that belong to one organism or genome together, attempting to achieve the highest completion and lowest contamination per bin. Make sure to switch back to the right Python environment.

<h3> <center> 8.1 Metabat2 </center> </h3>
By running metabat with the bam file as reference the accuracy of binning is improved. (creates a new directory with default naming in the current directory)

<code> runMetaBat.sh PATH-TO-ASSEMBLED-CONTIGS/contigs.fasta PATH-TO-COMPLEMENTING-BAM-FILE/pyr_d60_all_ mapped.bam </code>

<h3> <center> 8.2 Maxbin2 </center> </h3>

<code> run_MaxBin.pl -contig PATH-TO-ASSEMBLED-CONTIGS/contigs.fasta -out OUTPUT-PATH/maxbin_pyr_d60_all -thread num_threads -reads PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq -reads2 PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq </code>

<h3> <center> 8.3 CONCOCT </center> </h3>
Running CONCOCT requires a few more steps than the other two binning programs. CONCOCT also requires Python 3.0(?), so make sure to switch environments again.

<h3> <center> 8.3.1 Cutting the contigs into smaller pieces. </center> </h3>
The recommended cut length as specified <here> is 10000 bp. If necessary (i.e. in later stages of the bin targeted reassembly) cut length can be reduced as appropriate.
  
<code> cut_up_fasta.py PATH-TO-ASSEMBLED-CONTIGS/contigs.fasta -c 10000 -o 0 --merge_last -b OUTPUT-PATH/contigs_10k.bed > OUTPUT-PATH/contigs_10k.fa </code>

<h3> <center> 8.3.2 Generating the coverage table. </center> </h3>

<code> concoct_coverage_table.py PATH-TO-BED-FILE/contigs_10k.bed PATH-TO-MAPPED-SORTED-AND-INDEXED-CONTIGS/ pyr_d60_all_ mapped.bam > OUTPUT-PATH/coverage_table.tsv </code>

<h3> <center> 8.3.3 Binning. </center> </h3>
-b will create a new directory as output for concoct with the name you specify.

<code> concoct --composition_file PATH-TO-CUT-UP-CONTIGS/contigs_10k.fa --coverage_file PATH-TO-COVERAGE-TABLE/coverage_table.tsv -b PATH-TO-AND-NAME-FOR-OUTPUT-DIRECTORY/ </code>

<h3> <center> 8.3.4 Merging the cut up contigs with the original contigs </center> </h3>

<code> merge_cutup_clustering.py PATH-TO-CONCOCT-OUTPUT/clustering_gt1000.csv > OUTPUT-PATH/clustering_merged.csv </code>

<h3> <center> 8.3.5 Extracting the bins as FASTA files </center> </h3>

<code> mkdir PATH-TO-CONCOCT-OUTPUT/fasta_bins </code>
<code> extract_fasta_bins.py PATH-TO-ORIGINAL-CONTIGS/contigs.fasta PATH-TO-CONCOCT-OUTPUT/clustering_merged.csv --output_path PATH-TO-CONCOCT-OUTPUT/fasta_bins </code>

<h3> <center> 9. Bin refinement using METAWRAP </center> </h3>
METAWRAP takes the bins created by Metabat2, Maxbin2 and CONCOCT and bins them, again. It provides a combination of all three binning approaches to create 4 new sets of bins. The combinations are AB, ABC, AC and BC. This optimises the binning. Remember to document which binning program is used as which input (A=Metabat, …). -c enables to set a threshold for completeness and -x a threshold for contamination. The default is fairly conservative with -c 70 and -x 10, which can be a good approach for the first rebinning, but should be increased in later iterations. -o creates a new directory to be used as output for the refinement using the name you specify. <br>

<code> metawrap bin_refinement -o PATH-TO-METAWRAP-OUTPUT/Bin_Refinement -t num_threads -A PATH-TO-METABAT2-BINS/fasta_bins/ -B PATH-TO-MAXBIN-BINS/fasta_bins/ -C PATH-TO-CONCOCT-BINS/fasta_bins/ -c 90 -x 5  </code>

<h3> <center> 10. Classification of assembled reads </center> </h3>
To begin the bin targeted reassembly, it is necessary to know which bin contains contigs for which organism/gene/genome. Therefore, classification needs to be run on the created bins. Here, classification using the GTDB is conducted utilising GTDBtk. Note that using this database means that the GTDB names will be used and might therefore differ from the NCBI database organism names (i.e. Sva1033 (NCBI) = BM103 (GTDB)). -x sets the file format. If the contigs are .fa format change to -x fa. <br>

<code> gtdbtk classify_wf --genome_dir PATH-TO-REFINED-BINS/binsAB/ --out_dir OUTPUT-PATH/classified_pyr_d60_AB --cpus num_threads --extension fasta </code>
<code> gtdbtk classify_wf --genome_dir PATH-TO-REFINED-BINS/binsABC/ --out_dir OUTPUT-PATH/classified_py_d60_ABC --cpus num_threads --extension fasta </code>
…

<h3> <center> 11. Analysing statistics for the created bins using CheckM </center> </h3>
To know how complete and contaminated a bin is, CheckM is run. Bins with a completion below 70% and a contamination above 10% should be disregarded or refined. The higher the completeness and the lower the contamination the better the bin. This will not display the GTDB names but just the bin numbers. Check back with your GTDBtk file to know which bin is which.

<h3> <center> 11.1 Lineage_wf </center> </h3>
-x defines the input format (set -x fasta if bins are in .fasta), --tab_table is optional. If not set CheckM will just give out the results in the console, setting the flag will print out a file named as in -f in the set output directory. <br>

<code> checkm lineage_wf PATH-TO-FASTA-BINS/fasta_bins/ BachelorAnna/checkm/concoct/ -x fa -t 4 --tab_table -f checkm_out </code>

<h3> <center> 11.2 CheckM QA </center> </h3>
This step gives out the file (if flag --tab_table is set) that has the statisctical analysis of all bins located in the specified input folder (completeness, contamination, N50 value, ...)

<code> checkm qa -o 2 -f checkm_qa_out --tab_table -t 4 PATH-TO-CHECKM-OUTPUT/lineage.ms output_folder/ </code>

<h3> <center> 12. Bin targeted reassembly </center> </h3>
After having checked the bins and selected the one(s) suitable for your cause, bin targeted reassembly starts. This process also has a couple steps necessary.

<h3> <center> 12.1 remapping the clean reads to the selected bin using BBMap </center> </h3>
The minimum identity was set to 0.98 for the first iteration of mapping, and increased to 0.99 for the following iterations. <br>

<code> bbmap.sh in1=PATH-TO-CLEAN-READS/pyr_d60_all_1_clean.fq in2=PATH-TO-CLEAN-READS/pyr_d60_all_2_clean.fq minid=0.98 threads=num_threads outm=READ_98.sam ref=PATH-TO-BIN/BinX.fasta </code>

<h3> <center> 12.2 .sam > unsorted and unindexed .bam </center> </h3> 
<code> samtools view -S -b PATH-TO-MAPPED-FILE/READ_98.sam > OUTPUT-PATH/READ_98_r.bam </code>

<h3> <center> 12.2 raw .bam > .bam </center> </h3> 
Switch back to Anvi’o environment <br>

<code> anvi-init-bam -T 8 PATH-TO-RAW-BAM/READ98_r.bam -o OUTPUT-PATH/READ98.bam </code>

<h3> <center> 12.3 .bam > .fq </center> </h3>
Changes the .bam format to fasta/fastq format <br>

<code> samtools bam2fq PATH-TO-BAM-FILE/READ_98.bam > OUPUT-PATH/READ_98.fastq </code>

<h3> <center> 12.4 divide into two fastq files </center> </h3> <br>

<code> cat READ_98.fastq | grep '^@.*/1$' -A 3 --no-group-separator > READ_98_1.fastq </code> <br>
<code> cat READ_98.fastq | grep '^@.*/2$' -A 3 --no-group-separator > READ_98_2.fastq </code>

<h3> <center>  13. Reassembly with SPAdes </center> </h3> 
Utilising metagenome and paired end flags --meta --pe <br>

<code> spades.py -o OUTPUT-PATH/reassembly --pe1-1 PATH-TO-READ-1/READ_98_1.fastq --pe1-2 PATH-TO-READ-2/READ_98_2.fastq --meta -t num_threads </code> <br>
Moving on, map the clean reads back to the contigs created in the new assembly. Repeat steps 12.1 through 12.4. 

<font size=5> Steps 8 through 13 are repeated until the N50 value does no longer increase. </font>

<h3> <center>  14. Extracting 16s RNA from bins using Barrnap </center> </h3> 
Barrnap can identify bacterial, archeal and prokaryotic RNA, the default is bacterial. 

<code> --kingdom bac --threads num_threads --outseq OUTPUT-PATH/BC_Bin24_RNA.fa < PATH-TO-DESIRED-BIN/Refined_24.fa </code>

<h3> <center> 15. Annotating genes with RAST </center> </h3>
Send the obtained sequence to RAST for annotation. Depending on their job load and your sample it might take a while. 







