# A metagenomic workflow.

This project describes the workflow utilised in the Bachelorthesis "Metagenomic Analysis Of Representatives Of The Sva1033 Cluster; Potential Iron- And Sulfate Reducers". The analysed metagenomes were Pyr_d60, Cell_F9 and SF_I_II.

For this documentation the Pyr_d60 samples are used as examples for easier understanding of the command lines. 

<h3> <center> Preparations </center> </h3>

In order to use the metagenomic pipeline as described below, the installation of different programs is necessary when using the conda method and when not working in a pre set-up environment. Therefore, different environments running with different versions of Python are required. This project uses one main environment with Python Version 3.8 that contains most programs. <p>Furthermore four different environments were created using Python Versions 2.7.15, 2.7.18, 3.6 and 3.7.8. Prokka was installed in a clean environment without other programs. Therefore, switching between environments to use certain programs is necessary. <p> To deactivate an environment <code> conda deactivate </code>  is used, to activate a new environment <code> conda activate environment_name </code>  is used. Installation was conducted with <code> conda install program_name </code>  unless otherwise specified, installing all dependencies if not yet installed.

<h4> Environment 1: "metagenomics" </h4>
This is the main environment used for this project running on Python 3.8.5. The installed programs are:
- FastQC (v0.11.9) <p>
- BBTools (BBduk, BBMap) (Ver. 38.18) <p>
- MEGAHIT (v1.2.9) <p>
- Bowtie2 (Ver. 2.4.1) <p>
- Samtools (Ver. 1.7) <p>
- Metabat2 (Ver. 2.15) <p>
- Maxbin2 (Ver. 2.2.7) <p>
- GTDBtk (Ver. 1.3.0) <p>
- CheckM (Ver. 1.1.3)[Installation with pip3, see www.x, unzipped using <code> tar xvzf </code>] <p>
- SPAdes (Ver. 3.14.1) <p>
- Barrnap (Ver. 0.9) <p>

<h4> Environment 2: "python2.7" </h4>
This environment runs on Python Version 2.7.18 and contains the following two programs:
- Kraken2 (Ver. 2.0.9beta) <p>
- Bracken (Ver. 2.6.0)

<h4> Environment 3: "anvio" </h4>
This environment runs on Python 3.6.10 and contains:
- Anvi'o (Ver. 6.2) <p>
- CONCOCT (Ver. 1.1.0)

<h4> Environment 4: "metawrap" </h4>
This environment runs on Python 2.7.15 and contains:
- metaWRAP (Ver. 1.3.0)

<h4> Environment 5: "prokka" </h4>
This environment runs on Python 3.7.8 and contains:
- Prokka  (Ver. 1.14.6) <p>

Databases required for Kraken2 and GTDBtk were not built. The databases used were built and kindly provided by another member of the workgroup.


<h3> <center> 1. Quality check with FastQC </center> </h3> 
Using the Illumina Sequences obtained from environmental samples a quality check with FastQC is to be conducted. This analysis will show where the sequences need to be refined, such as trimming out adapters and potentially contaminated parts (mostly located at the beginning and end of the sequences). <p>

<code>fastqc -o OUTPUT_PATH -f fastq PATH/TO/RAW/READS/pyr_d60_all_1.fq PATH/TO/RAW/READS/pyr_d60_all_2.fq -t num_threads </code>

<h3> <center> 2. Trimming with BBDuk </center> </h3> 
After analysing the reads with FastQC, it was determined where trimming is necessary. Trimming is performed in 3 steps: Adapter trimming left, adapter trimming right and quality trimming. Always update the input file with the most recent file (i.e. using the right-trimmed file to perform the adapter trim on the left side of the sequence and then the rl-trimmed file to perform the quality trim on).

<h4> <center> 2.1 Right trim: </center> </h4>
in sets the input reads, out the outputfiles. ktrim defines which side the file is trimmed on (r,l), k defines that all kmers of size x (here 23) are used, mink sets the minimum kmer size (meaning that it will look for all kmers between k and mink), hdist sets the hamming distance, tpe Trims Pairs Evenly and tbo Trims By Overlap. <p>

<code> bbduk.sh t=num_threads in1=PATH/TO/RAW/READS/pyr_d60_all_1.fq in2=PATH/TO/RAW/READS/pyr_d60_all_2.fq out1=OUTPUT/PATH/pyr_d60_all_1_rtrim.fq out2=OUTPUT/PATH/pyr_d60_all_2_rtrim.fq ref=PATH/TO/ADAPTER/FILE/adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo  </code>

<h4> <center> 2.2 Left trim: </center> </h4>

<code>bbduk.sh t=num_threads in1=PATH/TO/RIGHT/TRIMMED/FILE/pyr_d60_all_1_rtrim.fq in2=PATH/TO/RIGHT/TRIMMED/FILE/pyr_d60_all_2_rtrim.fq out1=OUTPUT/PATH/pyr_d60_all_1_rltrimmed.fq out2=OUTPUT/PATH/pyr_d60_all_2_rltrimmed.fq ref=PATH/TO/ADAPTER/FILE/adapters.fa ktrim=l k=23 mink=11 hdist=1 tpe tbo </code>  

<h4> <center> 2.3 Quality trim: </center> </h4>
The parameters ftl (trim on left side), ftr (trim on right side), trimq (trims regions below set number), maq (trims regions below specified minimum average quality) and minlen (regions shorter than specified here are trimmed) have to be set according to the FastQC results. <p>

<code>bbduk.sh t=num_threads in1=PATH/TO/RL/TRIMMED/FILE/pyr_d60_all_1_clean.fq in2=PATH/TO/RL/TRIMMED/FILE/pyr_d60_all_2_clean.fq out1=OUTPUT/PATH/FOR/CLEAN/READS/pyr_d60_all_1_crisp.fq out2=OUTPUT/PATH/FOR/CLEAN/READS/pyr_d60_all_2_crisp.fq qtrim=rl trimq=20 ftl=6 ftr=144 maq=20 minlen=100 </code>

<h3> <center> 3. Repeat quality check with FastQC </center> </h3>

<h3> <center> 4. Abundance Estimation with Kraken2 and Bracken </center> </h3>
Abundance estimation gives a first overview of species identified in the sample. It is useful to check abundance before continuing the pipeline. When installing the packages for Kraken2 and Bracken check the packages in the environment using conda list (if using anaconda). Bracken can install kraken1 which will cause problems. If listed, remove before continuing. 

<h4> <center> 4.1 Kraken2 </center> </h4>
Kraken2 is the base for Bracken to run on. <p>

<code>kraken2 --db PATH/TO/KRAKEN/DATABASE/kraken2_db --paired --classified-out pyr_d60#.fq PATH/TO/CLEAN/READS/pyr_d60_all_1_clean.fq PATH/TO/CLEAN/READS/pyr_d60_all_2_clean.fq --threads num_threads --output PATH/TO/KRAKEN/OUTPUT/pyr_d60_Kraken.out --report Pyr_d60.report --confidence 0.05 </code>

<h4> <center> 4.2 Bracken </center> </h4>
Bracken now uses the the report Kraken created to estimate organism abundances. Set the level (-l) for each iteration (levels: D=domain, P=phylum, C=class, O=order, F=family, G=genus, S=species (default)). The parameter -r sets the readlength (use what was set as minlen in the quality trim). <p>

<code> bracken -d PATH/TO/KRAKEN/DATABASE/kraken2_db -i PATH/TO/KRAKEN/REPORT/FILE/Pyr_d60_.report -o PATH/TO/BRACKEN/OUTPUT/Pyr_d60.bracken -r 100 -l D -t num_threads </code> <p>
<code> bracken -d PATH/TO/KRAKEN/DATABASE/kraken2_db -i PATH/TO/KRAKEN/REPORT/FILE/Pyr_d60_.report -o PATH/TO/BRACKEN/OUTPUT/Pyr_d60.bracken -r 100 -l P -t num_threads </code>
...

<h3> <center> 5. <i>De novo</i> Assembly using Megahit </center> </h3>
The next step to take with the clean (quality trimmed) reads is to assemble them. Using overlapping kmers contigs are built which are then fit together to build scaffolds. -o creates a new directory. Setting -m (usable memory) limits memory usage. <p>

<code> megahit -1 PATH/TO/CLEAN/READS/pyr_d60_all_1_clean.fq -2 PATH/TO/CLEAN/READS/pyr_d60_all_2_clean.fq -o PATH/TO/ASSEMBLY/OUTPUT/pyr_d60_all.megahit_asm -m 0.2 -t num_threads </code>

<h3> <center> 6. Readmapping using Bowtie2 </center> </h3>
After having assembled the reads, Bowtie2 was used to map the assembly back to the reads.

<h4> <center> 6.1 Building a Bowtie2 Index </center> </h4>

<code> bowtie2-build --threads num_threads PATH/TO/ASSEMBLED/CONTIGS/pyr_d60_all.megahit_asm/final.contigs.fa OUTPUT/PATH/FOR/MAPPED-CONTIGS/contigs </code>

<h4> <center> 6.2 Running the mapping </center> </h4>

<code> bowtie2 --threads num_threads -x PATH/TO/NEWLY/CREATED/INDEX/DIRECTORY/contigsindex -1 PATH/TO/CLEAN/READS/pyr_d60_all_1_clean.fq -2 PATH/TO/CLEAN/READS/pyr_d60_all_2_clean.fq -S OUTPUT/PATH/FOR/MAPPED/FILE/pyr_d60_all_mapped.sam </code>

<h3> <center> 7. .sam > .bam </center> </h3>
Various further needed programs require sorted and indexed bam files to work. Sorting and indexing can be done using only samtools (samtools view, samtools sort and samtools index) or using a combination of samtools and Anvi'o. Here, the combination method is used.

<h4> <center> 7.1 Creating the "raw" (not sorted or indexed) .bam file. </center> </h4>

<code> samtools view --threads num_threads -F 4 -bS PATH/TO/MAPPED/FILE/pyr_d60_all_reform_mapped.sam > OUTPUT/PATH/pyr_d60_all_reform_mapped_RAW.bam </code>

<h4> <center> 7.2 Sorting and indexing the file (r_bam > .bam and .bam.bai) </center> </h4>

<code> anvi-init-bam -T num_threads PATH/TO/MAPPED/RAW/BAM/FILE/pyr_d60_all_mapped_RAW.bam -o OUTPUT/PATH/pyr_d60_all_mapped.bam </code>

<h3> <center> 8. Binning </center> </h3>
This step will sort the assembled contigs into bins. This means it will group contigs together, attempting to achieve the highest completion and lowest contamination per bin. <b> EDIT ME PLS </b>

<h4> <center> 8.1 Metabat2 </center> </h4>
By running metabat with the bam file as reference the accuracy of binning is improved. (creates a new directory with default naming in the current directory) <p>

<code> runMetaBat.sh PATH/TO/ASSEMBLED/CONTIGS/contigs.fasta PATH/TO/COMPLEMENTING/BAM/FILE/pyr_d60_all_ mapped.bam </code>

<h4> <center> 8.2 Maxbin2 </center> </h4>

<code> run_MaxBin.pl -contig PATH/TO/ASSEMBLED/CONTIGS/contigs.fasta -out OUTPUT/PATH/maxbin_pyr_d60_all -thread num_threads -reads PATH/TO/CLEAN/READS/pyr_d60_all_1_clean.fq -reads2 PATH/TO/CLEAN/READS/pyr_d60_all_1_clean.fq </code>

<h4> <center> 8.3 CONCOCT </center> </h4>

<h5> <center> 8.3.1 Cutting contigs into smaller pieces. </center> </h5>
The recommended cut length as specified <here> is 10000 bp. If necessary (i.e. in later stages of the bin targeted reassembly) cut length can be reduced as deemed appropriate. <p>
  
<code> cut_up_fasta.py PATH/TO/ASSEMBLED/CONTIGS/contigs.fasta -c 10000 -o 0 --merge_last -b OUTPUT/PATH/contigs_10k.bed > OUTPUT/PATH/contigs_10k.fa </code>

<h5> <center> 8.3.2 Generating a coverage table. </center> </h5>

<code> concoct_coverage_table.py PATH/TO/BED/FILE/contigs_10k.bed PATH/TO/MAPPED/SORTED/AND/INDEXED/CONTIGS/pyr_d60_all_ mapped.bam > OUTPUT/PATH/coverage_table.tsv </code>

<h5> <center> 8.3.3 Binning. </center> </h5>
-b will create a new directory as output for concoct with the name specified. <p>

<code> concoct --composition_file PATH/TO/CUT/UP/CONTIGS/contigs_10k.fa --coverage_file PATH/TO/COVERAGE/TABLE/coverage_table.tsv -b PATH/TO/AND/NAME/FOR/OUTPUT/DIRECTORY/ </code>

<h5> <center> 8.3.4 Merging cut up contigs with the original contigs </center> </h5>

<code> merge_cutup_clustering.py PATH/TO/CONCOCT/OUTPUT/clustering_gt1000.csv > OUTPUT/PATH/clustering_merged.csv </code>

<h5> <center> 8.3.5 Extracting the bins as FASTA files </center> </h5>

<code> mkdir PATH/TO/CONCOCT/OUTPUT/fasta_bins </code> <p>
<code> extract_fasta_bins.py PATH/TO/ORIGINAL/CONTIGS/contigs.fasta PATH/TO/CONCOCT/OUTPUT/clustering_merged.csv --output_path PATH/TO/CONCOCT/OUTPUT/fasta_bins </code>

<h3> <center> 9. Bin refinement using METAWRAP bin_refinement </center> </h3>
METAWRAP bin_refinement takes the bins created by Metabat2, Maxbin2 and CONCOCT and bins them, again. It provides a combination of all three binning approaches to create 4 new sets of bins. The combinations are AB, ABC, AC and BC. This optimises the binning. For all samples, a consistent input was used. (A=Metabat2, B=Maxbin2, C=CONCOCT). -c set a threshold for completeness and -x a threshold for contamination. -o creates a new directory to be used as output for the refinement using the specified name. <p>

<code> metawrap bin_refinement -o PATH/TO/METAWRAP/OUTPUT/Bin_Refinement -t num_threads -A PATH/TO/METABAT2/BINS/fasta_bins/ -B PATH/TO/MAXBIN2/BINS/fasta_bins/ -C PATH/TO/CONCOCT-BINS/fasta_bins/ -c 90 -x 5  </code>

<h3> <center> 10. Classification of assembled reads </center> </h3>
To begin the bin targeted reassembly, it is necessary to know which bin is assigned to the target organism. Therefore, classification needs to be run on the created bins. Here, classification using the GTDB is conducted utilising GTDBtk. Note that using this database means that the GTDB classification will be used and might therefore differ from the NCBI classification (i.e. Sva1033 (NCBI) = BM103 (GTDB)). --extension (-x) sets the file format. If the contigs are .fa format change to -x fa. <p>

<code> gtdbtk classify_wf --genome_dir PATH/TO/REFINED/BINS/binsAB/ --out_dir OUTPUT/PATH/classified_pyr_d60_AB --cpus num_threads --extension fasta </code> <p>
<code> gtdbtk classify_wf --genome_dir PATH/TO/REFINED/BINS/binsABC/ --out_dir OUTPUT/PATH/classified_py_d60_ABC --cpus num_threads --extension fasta </code>
…

<h3> <center> 11. Analysing statistics for the created bins using CheckM </center> </h3>
To compute completeness and contamination of a bin, CheckM is used. Bins with a completion below 70% and a contamination above 10% should be disregarded or refined. The higher the completeness and the lower the contamination the better the bin. 

<h3> <center> 11.1 Lineage_wf </center> </h3>
-x defines the input format (set -x fasta if bins are in .fasta), --tab_table is optional. If not set CheckM will give out the results in the console, setting the flag will print out a file named as in -f in the set output directory. <p>

<code> checkm lineage_wf PATH/TO/FASTA/BINS/fasta_bins/ OUTPUT/PATH/checkm/concoct/ -x fa -t 4 --tab_table -f checkm_out </code>

<h3> <center> 11.2 CheckM QA </center> </h3>
This step gives out the file (if flag --tab_table is set) that has the statistical analysis of all bins located in the specified input folder (completeness, contamination, N50 values, ...) <p>

<code> checkm qa -o 2 -f checkm_qa_out --tab_table -t 4 PATH/TO/CHECKM/OUTPUT/lineage.ms output_folder/ </code>

<h3> <center>  12. Extracting 16s RNA from bins using Barrnap </center> </h3> 
Barrnap identifies bacterial, archeal and prokaryotic RNA. A phylogeny tree can be created using ARB with the retrieved 16s RNA genes.<p>

<code> --kingdom bac --threads num_threads --outseq OUTPUT/PATH/BC_Bin24_RNA.fa < PATH/TO/DESIRED/BIN/Refined_24.fa </code>

<h3> <center> 13. Bin targeted reassembly </center> </h3>
After having checked the bins and selected the <b> EDIT ME PLS suitable one(s) EDIT ME PLS </b>, bin targeted reassembly starts. This process includes multiple steps.

<h3> <center> 13.1 Remapping the clean reads to the selected bin using BBMap </center> </h3>
The minimum identity was set to 0.98 for the first iteration of mapping, and increased to 0.99 for the following iterations. <p>

<code> bbmap.sh in1=PATH/TO/CLEAN/READS/pyr_d60_all_1_clean.fq in2=PATH/TO/CLEAN/READS/pyr_d60_all_2_clean.fq minid=0.98 threads=num_threads outm=READ_98.sam ref=PATH/TO/BIN/BinX.fasta </code>

<h3> <center> 13.2 .sam > unsorted and unindexed .bam (RAW) </center> </h3> 
<code> samtools view -S -b PATH/TO/MAPPED/FILE/READ_98.sam > OUTPUT/PATH/READ_98_r.bam </code>

<h3> <center> 13.2 RAW .bam > .bam (sorted and indexed) </center> </h3> 

<code> anvi-init-bam -T 8 PATH/TO/RAW/BAM/READ98_r.bam -o OUTPUT/PATH/READ98.bam </code>

<h3> <center> 13.3 .bam > .fq </center> </h3>
Changes the .bam format to fasta/fastq format. <p>

<code> samtools bam2fq PATH/TO/BAM/FILE/READ_98.bam > OUPUT/PATH/READ_98.fastq </code>

<h3> <center> 13.4 Divide into two fastq files </center> </h3> 

<code> cat READ_98.fastq | grep '^@.*/1$' -A 3 --no-group-separator > READ_98_1.fastq </code> <p>
<code> cat READ_98.fastq | grep '^@.*/2$' -A 3 --no-group-separator > READ_98_2.fastq </code>

<h3> <center>  14. Reassembly with SPAdes </center> </h3> 
Utilising metagenome and paired end flags --meta --pe <p>

<code> spades.py -o OUTPUT-PATH/reassembly --pe1-1 PATH-TO-READ-1/READ_98_1.fastq --pe1-2 PATH-TO-READ-2/READ_98_2.fastq --meta -t num_threads </code> <br>
Moving on, the clean reads were mapped back to the contigs created in the new assembly. Repeat steps 12.1 through 12.4. <p>

<font size=5> Steps 8 through 13 are repeated until the N50 scaffold value does no longer increase. </font>

<h3> <center> 15. Annotating genes with RAST </center> </h3>
Send the obtained sequence to RAST for annotation. 

<h3> <center> 16. Annotating genes with Prokka </center> </h3>
Headers need to be simplified before running Prokka. This was conducted with anvi-script-reformat-fasta. <p>
<code> anvi-script-reformat-fasta PATH/TO/final.contigs.fa -o OUTPUT/PATH/contigs-fixed.fa -l 0 --simplify-names -r OUTPUT/PATH/pyr_defline_report </code> <p>
<code> prokka --kingdom Bacteria --outdir OUTPUT/PATH/prokka --genus M0040 --locustag GCA_006226895.1 --cpus 4 --rnammer PATH/TO/SIMPLIFIED/FILE/simplify/Pyr_d60_AC4_1_simplified.fa </code> <p>
Genus and locustag can be obtained from the GTDB for the desired organism.






