#!/usr/bin/env python3
import os
import argparse
import subprocess
from pathlib import Path

def create_output_file(output, fastq):
    with open(output, "w") as data_report:
        data_report.write("# Report Summary\n\n")
        l_files = [f for f in os.listdir(fastq) if not f.startswith('.') and f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]
        for file in l_files:
            basename = Path(file).stem.split('.')[0]
    return output

def extract_number_of_reads(file_path):
    result = subprocess.run(f"zcat {file_path} | wc -l", shell=True, capture_output=True, text=True)
    num_reads = int(result.stdout.strip()) // 4
    return num_reads

def format_log_content(log_file):
    try:
        with open(log_file, "r") as log:
            lines = log.readlines()
        filtered_lines = [line for line in lines if "perl: warning" not in line]
        return "".join(filtered_lines)
    except FileNotFoundError:
        return f"Log file {log_file} not found.\n"

def add_log_content(log_files, basename, log_type, data_report):
    data_report.write(f"## {log_type} Statistics\n")
    found_log = False
    for log_file in log_files:
        if basename in Path(log_file).stem.split('.')[0]:  
            found_log = True
            content = format_log_content(log_file)
            data_report.write(content)
    
    if not found_log:
        data_report.write(f"No {log_type} log found for sample {basename}.\n")
    data_report.write("\n")

def process_samples(fastq, fastq_trimmed, unwanted, exome, genome, output):
    with open(output, "a") as data_report:
        l_files = [f for f in os.listdir(fastq) if not f.startswith('.') and f.endswith(('.fastq', '.fq', '.fastq.gz', '.fq.gz'))]

        for file in l_files:
            basename = Path(file).stem.split('.')[0]  
            file_path = os.path.join(fastq, file)
            num_reads_original = extract_number_of_reads(file_path)

            trimmed_reads_count = 0
            for trimmed_file in fastq_trimmed:
                trimmed_basename = Path(trimmed_file).stem.split('.')[0]
                if basename in trimmed_basename:
                    trimmed_reads_count = extract_number_of_reads(trimmed_file)
                    break
            
            data_report.write(f"### Sample: {basename}\n")
            data_report.write(f"- Reads in FASTQ: {num_reads_original}\n")
            data_report.write(f"- Reads after trimming: {trimmed_reads_count}\n")

            if unwanted:
                add_log_content(unwanted, basename, "Unmapped Sequences and Filtering", data_report)
            add_log_content(exome, basename, "Mapping Exome", data_report)
            add_log_content(genome, basename, "Mapping Genome", data_report)
            data_report.write(f"{'-'*40}\n")

def main():
    parser = argparse.ArgumentParser(description="Report Analysis")
    parser.add_argument('--fastq', required=True, help="Path to the directory containing original FASTQ files")
    parser.add_argument('--fastq_trimmed', required=True, nargs="*", help="Paths to the trimmed FASTQ files")
    parser.add_argument('--unwanted_sequence', nargs="*", default=[], help="Paths to mapping unwanted sequences log (optional)")
    parser.add_argument('--exome', nargs="*", help="Paths to exome mapping logs")
    parser.add_argument('--genome', nargs="*", help="Paths to genome mapping logs")
    parser.add_argument('--output', required=True, help="Path to the output report file")
    args = parser.parse_args()

    fastq = args.fastq
    fastq_trimmed = args.fastq_trimmed
    unwanted = args.unwanted_sequence
    exome = args.exome
    genome = args.genome
    output = args.output

    create_output_file(output, fastq)
    process_samples(fastq, fastq_trimmed, unwanted, exome, genome, output)

if __name__ == "__main__":
    main()

