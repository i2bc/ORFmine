import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import argparse


def read_file(tab):
    """Reads a tab-separated file and returns its lines as a list of parts."""
    data = []
    with open(tab, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            data.append(parts)
    return data


def read_files_from_directory(directory):
    """Reads all files in a directory and returns their data in a dictionary."""
    all_data = {}
    for name_file in os.listdir(directory):
        path_file = os.path.join(directory, name_file)
        if os.path.isfile(path_file):
            data = []
            with open(path_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    data.append(parts)
            all_data[path_file] = data
    return all_data


def generate_full_cds_name(cds_name):
    """Generates a full CDS name in the format 'CDS:<cds_name>_CDS'."""
    if cds_name.startswith("CDS:") and cds_name.endswith("_CDS"):
        return cds_name
    return f"CDS:{cds_name}_CDS"


def process_cds_from_file(cds_file):
    """Reads CDS names from a file and generates full CDS names."""
    full_cds_names = []
    with open(cds_file, 'r') as file:
        for line in file:
            cds = line.strip()
            full_cds_names.append(generate_full_cds_name(cds))
    return full_cds_names


def aggregate_data(cds_name, files, nb_cds=None):
    full_gene_name = generate_full_cds_name(cds_name)
    print(f"Looking for: {full_gene_name}")  # Debug print
    aggregate = {}
    all_data = []
    count = 0  # Initialize count here

    for file in files:
        print(f"Processing file: {file}")  # Debug print
        data = read_file(file)
        for line in data:
            if len(line) < 5:  # Ensure there are enough columns
                print(f"Skipping line (too few columns): {line}")
                continue

            gene_name_data, cds_num, P0, P1, P2 = line
            if gene_name_data == full_gene_name:
                cds_num = int(cds_num)
                P0 = int(P0)
                P1 = int(P1)
                P2 = int(P2)
                if nb_cds is not None and count >= nb_cds:
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


def plot_data(cds_name, all_data, pooled=False):
    """
    Generates a plot from the aggregated data for a specific CDS.
    
    Args:
        cds_name (str): Name of the CDS.
        all_data (list): Data to plot.
        pooled (bool): Whether the plot is generated in pooled mode.
    """
    cds_nums = [data[0] for data in all_data]
    P0_vals = [data[1] for data in all_data]
    P1_vals = [data[2] for data in all_data]
    P2_vals = [data[3] for data in all_data]

    bar_width = 0.15
    plt.figure(figsize=(36, 18))
    plt.bar(np.array(cds_nums) - bar_width, P0_vals, width=bar_width, align='center', label=f'{cds_name} in-frame reads', color='orange')
    plt.bar(cds_nums, P1_vals, width=bar_width, align='center', label=f'{cds_name} out-of-frame reads', color='blue')
    plt.bar(np.array(cds_nums) + bar_width, P2_vals, width=bar_width, align='center', color='blue')

    plt.xlabel('Positions', fontsize=14)
    plt.ylabel('Reads Number', fontsize=14)
    plt.legend(fontsize='large')
    plt.title(f'{cds_name}', fontsize=14)

    max_cds_num = max(cds_nums)
    tick_positions = np.arange(0, (max_cds_num // 5 + 1) * 5, 5)
    plt.xticks(tick_positions, rotation=90)

    plt.gca().set_xlim(left=0)
    clean_cds_name = cds_name.replace("CDS:", "").replace("_CDS", "")
    if pooled:
        clean_cds_name += "_pooled"
    
    output_filename = f"{clean_cds_name}.png"
    plt.savefig(output_filename)
    print(f"Plot saved to {output_filename}")


def main():
    parser = argparse.ArgumentParser(description="Script to generate periodicity plots based on CDS data.")
    parser.add_argument("--tab", type=str, help="Path to a tab file if using individual mode.")
    parser.add_argument("--directory", type=str, help="Path to the directory containing tab files if using pooled mode.")
    parser.add_argument("--nb_cds", type=int, help="Number of CDS to process (default: all).")
    parser.add_argument("--cds_file", type=str, help="Path to a file containing CDS names (in .txt format).")
    parser.add_argument("--cds_name", type=str, help="Name of the CDS you want to use.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--pooled", action='store_true', help="Use all the files in the directory and generate one plot for all data.")
    mode.add_argument("--individual", action='store_true', help="Use only the file provided in the command line to generate a plot.")

    args = parser.parse_args()

    # Validate input based on selected mode
    if args.pooled and not args.directory:
        print("Error: Directory path is required for pooled mode.")
        sys.exit(1)
    if args.individual and not args.tab:
        print("Error: A tab file is required for individual mode.")
        sys.exit(1)

    if args.pooled:
        print("Pooled mode selected...")
        all_files = read_files_from_directory(args.directory)
        files = list(all_files.keys())
    else:
        print("Individual mode selected...")
        files = [args.tab]

    # Process CDS names
    if args.cds_file:
        cds_names = process_cds_from_file(args.cds_file)
    elif args.cds_name:
        cds_names = [generate_full_cds_name(args.cds_name)]
    else:
        print("Error: You must provide either a CDS name or a file with CDS names.")
        sys.exit(1)

    for cds in cds_names:
        aggregate, all_data = aggregate_data(cds, files, args.nb_cds)
        if all_data:
            plot_data(cds, all_data, pooled=args.pooled)
        else:
            print(f"No data found for {cds}")


if __name__ == "__main__":
    main()
