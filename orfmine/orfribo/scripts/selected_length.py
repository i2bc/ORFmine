import pandas as pd
import re
import argparse
import glob
import os

def main():
    parser = argparse.ArgumentParser(description="Select maximum mean or median value from statistics files.")
    parser.add_argument('--directory', type=str, help="The directory containing statistics files.")
    parser.add_argument('--mean', action='store_true', help="Select maximum mean value.")
    parser.add_argument('--median', action='store_true', help="Select maximum median value.")
    parser.add_argument('--both', action='store_true', help="Select maximum mean and median value.")
    parser.add_argument('--threshold', type=float, help="Threshold for selection.")
    parser.add_argument('--output', type=str, help="Path to output file.")
    args = parser.parse_args()

    directory = args.directory
    mean = args.mean
    median = args.median
    both = args.both
    threshold = args.threshold
    output = args.output

    # Find all .stats files in the specified directory and its subdirectories
    files = glob.glob(os.path.join(directory, '**', '*.stats'), recursive=True)

    def reading_file(files):
        for filename in files:
            df = pd.read_csv(filename, sep="\s+", header=0)
            num = re.search(r"Exome_(\d+)_reads\.stats", filename).group(1)
            yield df, num

    def extract_max_mean(files, threshold):
        max_files_mean = []
        for df, num in reading_file(files):
            mean_PO = df[df["Metric"] == "Mean"]["P0"].values[0]
            if mean_PO > threshold:
                max_files_mean.append((num, mean_PO))
        return max_files_mean

    def extract_max_median(files, threshold):
        max_files_median = []
        for df, num in reading_file(files):
            median_PO = df[df["Metric"] == "Median"]["P0"].values[0]
            if median_PO > threshold:
                max_files_median.append((num, median_PO))
        return max_files_median

    if mean:
        max_files_mean = [file_mean[0] for file_mean in extract_max_mean(files, threshold)]
        with open(output, "w") as file:
            for file_name in max_files_mean:
                file.write(file_name + "\n")

    if median:
        max_files_median = [file_median[0] for file_median in extract_max_median(files, threshold)]
        with open(output, "w") as file:
            for file_name in max_files_median:
                file.write(file_name + "\n")
    
    if both:
        max_files_mean = set([file_mean[0] for file_mean in extract_max_mean(files, threshold)])
        max_files_median = set([file_median[0] for file_median in extract_max_median(files, threshold)])
        combined_files = max_files_mean.intersection(max_files_median)
        with open(output, "w") as file:
            for file_name in combined_files:
                file.write(file_name + "\n")

if __name__ == "__main__":
    main()
