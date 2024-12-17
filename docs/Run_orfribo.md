## Basic Run ORFribo

ORFribo is very user-friendly and offers two configuration options:

### 1. Using Command Line:

You can configure ORFribo directly via the command line using the following example:

```bash
orfribo --fna examples/database/Scer.fna --gff examples/database/Scer.gff --gff-intergenic examples/database/mapping_orf_Scer.gff --fastq /data/work/I2BC/fadwa.elkhaddar/BIM/fastq/ --not-trimmed
```

Make sure to communicate if your data are trimmed or not, if they are trimmed, replace --not-trimmed by trimmed. 
By Default, ORFribo will consider your data as not-trimmed. 


---

### 2. Using the Configuration File:

The `config.yaml` file contains parameters that users can modify. Below is an explanation of how to complete this file. (Example datasets are available in `ORFmine/examples/`.)

**Important**: We strongly recommend reading the [How it works?](./How_it_works_orfribo.md) page before editing the configuration file.

---

### Configuration File Details

- **Project Name (Optional):**
```yaml
project_name: ""
```
*(Provide a project name without spaces or special characters.)*

- **Genome/Transcriptome FASTA File Path (Required):**
```yaml
fna: ""
```
*(Example: `reference_genome_sequences.fa`: A FASTA file containing the nucleotide sequence of the complete genome.)*

- **Reference Annotation GFF File Path (Required):**
```yaml
gff: ""
```
*(Example: `ORFtrack_output.gff`: ORFtrack output GFF file with ORF coordinates for translation activity analysis.)*

- **Intergenic Annotation GFF File Path (Required):**
```yaml
gff_intergenic: ""
```
*(Example: `mapping_orf_Scer.gff`: GFF file mapping intergenic regions from ORFtrack.)*

- **FASTA File with Sequences to Exclude (Optional):**
```yaml
rna_to_exclude: ""
```
*(Example: `NA_sequences_to_remove.fa`: File containing sequences to exclude, like rRNAs.)*

- **Directory Containing FASTQ Files (Required):**
```yaml
fastq: ""
```

- **Output Directory Path (Optional):**
```yaml
out: ""
```
*(Default: `./orfribo_%datetime`.)*

---

Remember to save your changes in the `config.yaml` file before running ORFribo. For further guidance.
if you want to run orfribo using the config file, you can execute this command line.

```bash
orfribo  --config path/to/config.yml
```

