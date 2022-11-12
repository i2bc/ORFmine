## Run on multiple datasets

ORFribo can handle multiple Ribosome Profiling datasets provided they correspond to the same genome. To do so, the user only need to store all the corresponding fastq files in the /fastq/ directory. ORFribo will treat all the datasets independently and will generate as [many output tables](./orfribo_outputs.md#main-output-table) as there are input fastq files. The generated main output tables are then merged together into a single table written in the /workdir/orfribo/RESULTS/Bam2Reads_genome_output/ and is named "all_samples_genome_meanXX_medianYY_reads_concatenated.tab".

We recall that the config.yaml file must be [completed](./orfribo_configuration.md) and stored in the /workdir/ directory of the container. One should note that the same parameters will be applied to all datasets.

If you want to add some files to your experiment when ORFribo is already done for some fastqs, just put the new fastqs in the /fastq/ folder and launch the analysis again. ORFribo is smart enough (!) to only execute the needed steps for your analysis. The same logic applies if you change some parameters in the configuration file <i>config.yaml</i> or change the input files. If a file is modified in any way by the user or if a newer one takes its place, the downstream steps of the pipeline based on the modified file will be executed but no more (based on files' timestamp).
For example, you can also wait for the analysis to be finished, then change the selection thresholds in your configuration file and start ORFribo again. The new outputs will be created only after the necessary steps will be made, without trimming the adapters and aligning the reads on the genome's sequence once again.


<a name="CDS"></a>
## Probing the translation activity of CDSs
ORFmine is ORF-centered (i.e. ORFs defined from STOP-to-STOP), thus the ORFtrack output, that is used by most of ORFmine tools including ORFribo for the read mapping, only contains ORFs as defined by ORFtrack. Consequently, the ORFtrack output contains ORFs that include a CDS or CDS exons but not the real intronless CDSs (ATG-STOP) as defined in the original annotation gff file. If the user wants to probe the translation of CDSs, he/she therefore needs to extract the information that is stored in the outputs of step 2 (the kmer quality control performed on CDSs). To do so, the user has to identify for the studied dataset (e.g. dataset_XYZ), the retained kmers which are stored in this file: /workdir/orfribo/RESULTS/selected_tables/threshold_meanXX_medianXX/dataset_XYZ.txt. This file lists all the retained kmers. The CDS read counts are written in the corresponding tables stored in /workdir/orfribo/RESULTS/ORFribo/Bam2Reads_exome_output/dataset_XYZ_kmer_i/ with kmer_i standing for the retained kmers. The user thus needs to pool the corresponding count tables into a single table of same format as [the final output](./orfibo_outputs.md#main_table) through our script as described below.

If for example, for the dataset_XYZ for which the retained kmers were the 27mers, 28mers and 31mers, the user can pool all the count tables based on the mapping of the 27-, 28- and 31mers as follows:

``` python
python3 /ORFmine/orfribo/RiboDoc_BAM2Reads/tools/Bam2Reads_function/concatenate.py -tables /workdir/orfribo/RESULTS/ORFribo/Bam2Reads_exome_output/dataset_XYZ_27/exome.25-35_reads.tab  /workdir/orfribo/RESULTS/ORFribo/Bam2Reads_exome_output/dataset_XYZ_28/exome.25-35_reads.tab  /workdir/orfribo/RESULTS/ORFribo/Bam2Reads_exome_output/dataset_XYZ_31/exome.25-35_reads.tab -outpath /workdir/orfribo/ -outname CDS_translation_XYZ.tab
```

The resulting table contains for each CDS, the numbers and frequencies of reads in its F0, F1 anf F2 frames.
