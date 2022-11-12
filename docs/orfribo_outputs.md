## Output folder architecture (/workdir/orfribo/):
Here is the folder architecture of the ORFribo output stored in /workdir/orfribo.   

The following output example is based on the example provided in the /ORFmine/examples/ directory with the fastq SRR1520313_17031088.


<pre>
/workdir/orfribo/
├── <i>dag_all.svg</i>
├── <i>dag_last_run.svg</i>
├── RESULTS/
	├── <i>config.yaml</i>
	├── <i>ORFribo_yeast_example.Analysis_Report.txt</i>
	├── adapter_lists/
		└── <i>SRR1520313_17031088.txt</i>
	├── Bam2Reads_genome_output/
		├── <b>all_samples_genome.25-35.mean70_median70_reads_concatenated.tab</b>
        └── <i>SRR1520313_17031088</i>
            ├── <b>genome.25-35.mean70_median70_reads_concatenated.tab</b>
            ├── <i>length_25</i>
                ├── <i>genome.25-35.mean70_median70_reads.tab</i>
                ├── <i>genome.25-35.mean70_median70_periodicity_all.tab</i>
                ├── <i>genome.25-35.mean70_median70_periodicity_start.tab</i>
                └── <i>genome.25-35.mean70_median70_periodicity_stop.tab</i>
            ├── <i>length_26</i>
                ├── <i>genome.25-35.mean70_median70_reads.tab</i>
            ├ ...
            └── <i>length_35</i>
                └── <i>genome.25-35.mean70_median70_reads.tab</i>
	├── annex_database/
		├── <i>NamedCDS_Scer.gff</i>
		├── <i>index_bowtie2.1.bt2</i>
        ├ ...
		├── <i>index_hisat2.1.ht2</i>
        ├ ...
		├── <i>outRNA_bowtie2.1.ht2</i>
        └ ...
	├── fastqc/
		├── <i>SRR1520313_17031088_fastqc.html</i>
		└── <i>SRR1520313_17031088_fastqc.zip</i>
    ├── selected_tables/
        └── <i>threshold_mean70_median70</i>
        	└── <i>SRR1520313_17031088.txt</i>
	├── BAM/
		├── <i>SRR1520313_17031088.25-35.bam</i>
		└── <i>SRR1520313_17031088.25-35.bam.bai</i>
	├── no-outRNA/
        └── <i>SRR1520313_17031088.25-35.no-outRNA.fastq.gz</i>
	└── ORFribo/
		├── database/
    		├── <i>exome.nfasta</i>
            ├── <i>exome_elongated.nfasta</i>
            ├── <i>exome_elongated.nfasta.fai</i>
        	├── <i>exome_elongated.gff</i>
        	└── <i>exome_elongated_with_gene_features.gff</i>
		├── annex_database/
    		├── <i>exome_elongated.exons_Scer.fna</i>
    		├── <i>exome_elongated.exons_Scer.gff.gtf</i>
    		├── <i>exome_elongated.exome_index_bowtie2.1.bt2</i>
    		├── ...
            ├── <i>exome_elongated.exome_index_hisat2.1.ht2</i>
    		└── ...
    	├── BAM_exome/
    		├── <i>exome_elongated.SRR1520313_17031088.25-35.bam</i>
    		└── <i>exome_elongated.SRR1520313_17031088.25-35.bam.bai</i>
		├── Bam2Reads_exome_output/
            ├── <i>SRR1520313_17031088_25</i>
                ├── <b>exome.25-35_reads.tab</b>
                ├── <i>exome.25-35_periodicity_all.tab</i>
                ├── <i>exome.25-35_periodicity_start.tab</i>
                ├── <i>exome.25-35_periodicity_stop.tab</i>
                └── <i>exome.25-35_reads.stats</i>
            ├── ...
            └── <i>SRR1520313_17031088_35</i>
                ├── <b>exome.25-35_reads.tab</b>
                ├── <i>exome.25-35_periodicity_all.tab</i>
                ├── <i>exome.25-35_periodicity_start.tab</i>
                ├── <i>exome.25-35_periodicity_stop.tab</i>
                └── <i>exome.25-35_reads.stats</i>
		└── riboWaltz/
		      ├── <i>psite_offset.csv</i>
		      ├── <i>best_offset.txt</i>
		      └── <i>exome_elongated.SRR1520313_17031088</i>
		            ├── <i>21.tiff</i>
		            ├── <i>22.tiff</i>
		            ├── ...
		            └── <i>35.tiff</i>
├── benchmarks/
    ├── <i>adapt_trimming</i>
        └── <i>SRR1520313_17031088.benchmark.txt</i>
    └── ...
├── logs/
	├── <i>RiboDoc_package_versions.txt</i>
	├── <i>adapt_trimming</i>
	   ├── <i>SRR1520313_17031088_cutadapt.log</i>
	   └── <i>SRR1520313_17031088_trim_value.log</i>
   └── ...
└── logsTmp/
	├── <i>SRR1520313_17031088_adapt_trimming.log</i>
	├── <i>SRR1520313_17031088_bowtie2_run_outRNA.log</i>
	├── <i>SRR1520313_17031088_run_mapping_bowtie2.log</i>
	└── <i>SRR1520313_17031088_run_mapping_hisat2.log</i>
</pre>



<a name="main-outable-table"></a>

## Main output table
The main output table can be found in the <i>RESULTS/Bam2Reads_genome_output/</i> folder and is named <b>all_samples_genome.25-35.meanXX_medianYY_reads_concatenated.tab</b> (with X and Y standing for the median and mean thresholds, here 70 and 70 respectively). This table summarizes the results for each ORF (i.e. numbers and fraction of reads in the three frames of the ORF)  calculated for all the retained kmers in each input dataset (in this example, there is only one dataset - SRR1520313_17031088.). If ORFribo has been performed on [multiple datasets](./orfribo_advanced.md), it also provides for each dataset a folder named according to the dataset that contains a table with the same format as the <b>all_samples_genome.25-35.meanXX_medianYY_reads_concatenated.tab</b> table but calculated only on the [retained kmers](./How_it_works_orfribo.md#Step2) of the corresponding dataset <i>dataset_XYZ/genome.25-35.mean70_median70_reads_concatenated.tab</i>. The directories of each dataset also contain the intermediate tables calculated for each retained kmer. The latter are stored in the <i>RESULTS/Bam2Reads_genome_output/dataset_xyz/length_i/</i> directory where i stands for the kmer sizes. When only one dataset has been provided to ORFribo (as in this example), the output table present in the dataset folder is the same as the <b>all_samples_genome.25-35.meanXX_medianYY_reads_concatenated.tab</b> present in the Bam2Reads_genome_output directory. 

The summary table <b>all_samples_genome.25-35.meanXX_medianYY_reads_concatenated.tab</b> contains for each ORF of the studied category(ies), the number and fraction of reads in its three frames (frame 0 (F0 or P0) and its two alternative frames +1 and +2 (P1/F1 and P2/F2 respectively)). This table has 8 columns:

<i>Seq_ID</i> : Identifier of the ORF <br>
<i>Num_reads</i> : Number of reads having a P-site aligned on the ORF <br>
<i>Num_p0</i> : Number of reads with their P-site in phase 0 of the ORF  <br>
<i>Num_p1</i> : Number of reads with their P-site in phase 1 of the ORF <br>
<i>Num_p2</i> : Number of reads with their P-site in phase 2 of the ORF <br>
<i>Perc_p0</i> : Percentage of reads aligned on the ORF with their P-site in phase 0 (i.e. Num_p0 divided by the total number of reads mapping on the ORF) <br>
<i>Perc_p1</i> : Percentage of reads aligned on the ORF with their P-site in phase 1 <br>
<i>Perc_p2</i> : Percentage of reads aligned on the ORF with their P-site in phase 2 <br>


## Output folder details :

* The <i>dag files</i> which represents the analysis steps with your datasets.  

* The <i>logs/</i> folder groups together all the error output messages from tools used in ORFribo analysis pipeline. Thus, in the event of an error, it allows you to identify the problematic step (and give us feedback if needed).

* The <i>RESULTS/</i> folder contains these files and folders: <br>
* i) <i>PROJECT_NAME.Analysis_report.txt</i> gathers standard output of each analysis pipeline tool. It allows to know how many reads are present at each step of the analysis :  a)raw reads b)reads after trimming and length selection c)after out RNA depletion d)after double alignment on the reference genome.

* *ii) <i>config.yaml</i> to have a backup of the parameters.

* *I) <b>Bam2Reads_genome_output/</b>: Contains the final outputs of the analysis (explained above).

* *II) <i>annex_database/</i>: It contains the indexes for the genome alignment and the gff with all CDS named.

* III) <i>fastqc/</i>: It contains data quality controls.

* IV) <i>adapter_lists/</i>: It contains a text file with the adapters list for each dataset that were found in the <i>config.yaml</i> file or determined from data if the user did not provide any adapter sequence in the configuration file.

* V) <i>selected_length_tables/</i>: It contains one subfolder with a file for every threshold the user chose for the alignment on CDSs step (median or mean of P0 proportions) as multiple thresholds may be tried. The file keeps the information of which read length passed the threshold and was kept for the alignment on all ORFs.

* VI) <i>BAM/</i>: It contains a BAM file for each dataset (allows visualization on tools such as IGV).

* VII) <i>no-outRNA/</i>: It contains fastq files trimmed and after removal of the reads aligned on unwaned sequences.

* VIII) <i>ORFribo/</i>:

* * I) <b>Bam2Reads_exome_output/</b>: It contains one subfolder per kmer (dataset_XYZ_i with i standing for the kmer size) that itsefl contains a table <b>exome.25-35_reads.tab</b>. The information in this table is the same as that contained in the <b>all_samples_genome.25-35.meanXX_medianYY_reads_concatenated.tab</b> table described above but only for the CDSs annotated in the original gff file. Each CDS is associated with the fractions of reads that map on its coding frame (i.e. in-frame reads named also P0 reads) or in its alternative frames. These tables are used for the [detection of good quality kmers (Step 2)](./How_it_works_orfribo.md#Step2) but can also be used to probe the translation activity of CDSs. In this case, you just need to merge the tables of all retained kmers (the list of the retained kmers can be found in the file /workdir/orfribo/RESULTS/selected_tables/threshold_meanXX_medianXX/dataset_XYZ.txt) into a global table (see [here](./orfribo_advanced.md#CDS) for more details). This final output table will contain for each CDS, the number and fraction of reads mapping in-frame or in the +1 and +2 frames of each CDS.

* * II) database/: It contains the re-formatted fasta and gff with artificial elongated CDSs to avoid the missing of reads which align on the borders of CDSs (i.e. on the start and stop codons).

* * III) annex_database/: It contains the indexes for the exome alignments (alignments on all exons of CDSs as indicated in the original gff file) and the gtf for riboWaltz.

* * IV) BAM_exome/: It contains a BAM file for each dataset corresponding to the alignment on the fasta with transcript by transcript feature.

* * V) riboWaltz/: It contains the P-site offsets file.


## Particular case of multiple inputs

ORFribo can be launched on multiple datasets (fastqs) as long as the input files concern the same organism - i.e. share the same reference fasta and gff files (more details [here](./orfribo_advanced.md)). For this, you just have to put all the fastq files in the /fastq/ folder and launch ORFribo as you usually do with:
``` bash
orfribo CPU MEMORY
```

The output architecture will be the same as the one obtained for a single input, except that all intermediate files will be stored in subdirectories corresponding to each input dataset as follows: 
<pre>
orfribo/
├── <i>dag_all.svg</i>
├── <i>dag_last_run.svg</i>
├── RESULTS/
	├── <i>config.yaml</i>
	├── adapter_lists/
		└── <i>one_file_by_sample.txt</i>
	├── Bam2Reads_genome_output/
		├── all_samples_genome.25-35.mean70_median70_reads_concatenated.tab"
        └── <i>one_folder_by_sample</i>
            ├── <i>concatenated_results_table.tab</i>
            └── <i>one_folder_by_length</i>
                └── <i>Bam2Reads_results_for_all_ORFs_alignment</i>
    ├── annex_database/
		├── <i>gff_files_with_named_CDSs.gff</i>
		├── <i>indexes_for_bowtie2_alignments.bt2</i>
		├── <i>indexes_for_hisat2_alignments.ht2</i>
    ├── fastqc/
        ├── <i>one_html_by_sample.html</i>
		└── <i>one_zip_by_sample.zip</i>
    ├── selected_tables/
        └── <i>one_folder_by_sample</i>
	├── BAM/
		├── <i>one_bam_by_sample.bam</i>
		└── <i>one_bai_by_bam.bai</i>
	├── no-outRNA/
        └── <i>one_file_by_sample.fastq.gz</i>
	└── ORFribo/
		├── database/
    		├── <i>intermediate_fasta_files.fa</i>
        	├── <i>intermediate_gff_files.gff</i>
		├── annex_database/
    		├── <i>intermediate_fasta_files.fa</i>
    		├── <i>intermediate_gtf_file_for_riboWaltz.gtf</i>
    		├── <i>indexes_for_bowtie2_alignments.bt2</i>
    		└── <i>indexes_for_hisat2_alignments.ht2</i>
    	├── BAM_exome_output/
    		├── <i>one_bam_by_sample.bam</i>
    		└── <i>one_bai_by_bam.bai</i>
		├── Bam2Reads_exome_output/
            └── <i>one_folder_by_sample_and_length</i>
                └── <i>Bam2Reads_results_for_CDS_alignment</i>
		└── riboWaltz/
		      └── <i>riboWaltz's qualitative analysis results</i>
├── benchmarks/
    ├── <i>one_benchmark_folder_by_rule</i>
        └── <i>one_benchmark_file_by_job</i>
├── logs/
	└── <i>one_log_folder_by_rule</i>
        └── <i>one_log_file_by_job_and_command</i>
├── logsTmp/
	└── <i>one_file_by_steps_of_interest_for_alignment_stats</i>
</pre>


