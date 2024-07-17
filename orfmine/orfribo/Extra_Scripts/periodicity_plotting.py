import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
import argparse

def read_file(file):
    """Reads each file and splits each line into parts using tab characters as separators.
    
    Args: 
        file (str): The file path to read.
        
    Returns:
        list: A list which stores the parts.
    """
    data = []
    with open(file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            data.append(parts)
    return data

def generate_full_gene_name(gene):
    """Generate a full gene name in the format 'CDS:<gene>_CDS'.
    
    Args:
        gene (str): The name of the gene.
        
    Returns:
        str: The full gene name in the format 'CDS:<gene>_CDS'.
    """
    return f"CDS:{gene}_CDS"

def aggregate_data(gene, files, cds=None):
    """
    Aggregate data for a specified gene from a list of files.
    
    Args:
        gene (str): The name of the gene.
        files (list of str): List of file paths to read data from.
        cds (int, optional): Number of coding sequences (CDS) to process. If None, process all CDS.
        
    Returns:
        tuple: A dictionary with aggregated data and a list of all processed data.
    """
    full_gene_name = generate_full_gene_name(gene)
    aggregate = {}
    all_data = []
    for file in files:
        data = read_file(file)
        count = 0
        for line in data:
            gene_name_data, cds_num, P0, P1, P2 = line
            if gene_name_data == full_gene_name:
                cds_num = int(cds_num)
                P0 = int(P0)
                P1 = int(P1)
                P2 = int(P2)
                if cds is not None and count >= cds:
                    break
                if cds_num not in aggregate:
                    aggregate[cds_num] = [P0, P1, P2]
                else:
                    aggregate[cds_num][0] += P0
                    aggregate[cds_num][1] += P1
                    aggregate[cds_num][2] += P2
                all_data.append((cds_num, P0, P1, P2))
                count += 1
    return aggregate, all_data

def main():
    parser = argparse.ArgumentParser(description="Plot periodicity tables.")
    parser.add_argument("--gene", type=str, help="The gene name to filter by (e.g., YAL003W).")
    parser.add_argument("--files", nargs='+', help="The list of files to process (.tab)")
    parser.add_argument("--cds", type=int, help="The number of first lines to filter by (optional), if not defined, it will use all the CDS in that file")
    parser.add_argument("--output", type=str, help="The output file to save the plot (e.g., plot.png).", required=True)

    args = parser.parse_args()

    gene = args.gene
    files = args.files
    cds = args.cds
    output_file = args.output

    aggregate, all_data = aggregate_data(gene, files, cds)

    cds_nums = [data[0] for data in all_data]
    P0_vals = [data[1] for data in all_data]
    P1_vals = [data[2] for data in all_data]
    P2_vals = [data[3] for data in all_data]

    bar_width = 0.15  
    plt.figure(figsize=(36, 18))
    plt.bar(np.array(cds_nums) - bar_width, P0_vals, width=bar_width, align='center', label=f'{gene} in-frame reads', color='orange')
    plt.bar(cds_nums, P1_vals, width=bar_width, align='center', label=f'{gene} out-of-frame reads', color='blue')
    plt.bar(np.array(cds_nums) + bar_width, P2_vals, width=bar_width, align='center', color='blue')

    plt.xlabel('Positions', fontsize=14)
    plt.ylabel('Reads Number', fontsize=14)
    plt.legend(fontsize='large')
    plt.title(f'{gene}', fontsize=14)

    max_cds_num = max(cds_nums)
    tick_positions = np.arange(0, (max_cds_num // 5 + 1) * 5, 5)
    plt.xticks(tick_positions, rotation=90)

    plt.gca().set_xlim(left=0)

    plt.savefig(output_file)
    print(f"Plot saved to {output_file}")

if __name__ == "__main__":
    main()

