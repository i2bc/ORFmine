import argparse
import sys
from typing import List


def get_args() -> argparse.Namespace:
    """Return command line parameters

    Returns:
        argparse.Namespace: command line parameters
    """
    parser = get_parser()
    args = parser.parse_args()

    if args.config:
        args.config = { x.split("=")[0]:x.split("=")[1] for x in args.config }

    args.mem_mb = args.mem_mb * 1000
    
    return args


def get_parser():
    parser = argparse.ArgumentParser(description='Run ORFribo pipeline.')

    parser.add_argument("-c", "--config", help="Path to a YAML configuration file for ORFribo.")
    parser.add_argument("--fna", help="Path to the genome/transcriptome fasta file")
    parser.add_argument("--gff", help="Path to the GFF annotation file")
    parser.add_argument("--gff_intergenic", help="Path to the GFF annotation mapping file")
    parser.add_argument("--path_to_fastq", help="Path to the directory containing .fastq.gz files")
    parser.add_argument("--project_name", help="Name for the experiment", default="")
    parser.add_argument("--out_base", help="Base directory location for orfribo outputs", default="")
    parser.add_argument("--fasta_outRNA", help="Path to a fasta file with sequences to exclude", default="")
    parser.add_argument("--already_trimmed", choices=["yes", "no"], help="Are sequencing adapters already removed")
    parser.add_argument("--adapt_sequence", help="Adapter sequence (e.g. 'AGATCGGAAGAGCACACGTCT')", default="")
    parser.add_argument("--readsLength_min", help="Minimum read length for ribosome profiling", type=int, default=25)
    parser.add_argument("--readsLength_max", help="Maximum read length for ribosome profiling", type=int, default=35)
    parser.add_argument("--gff_cds_feature", help="Elements to select during counting", default="CDS")
    parser.add_argument("--gff_name_attribute", help="Attribute containing gene names in GFF", default="Name")
    parser.add_argument("--orfstats_mean_threshold", help="Minimum mean of in-frame reads", type=int, default=70)
    parser.add_argument("--orfstats_median_threshold", help="Minimum median of in-frame reads", type=int, default=70)
    parser.add_argument("--final_counts", help="List of features in intergenic GFF", default="nc_intergenic")
    parser.add_argument("--mem_mb", help="Maximum allowed RAM to use in Mb (default: 2000Mb).", type=int, default=2000)
    parser.add_argument("--threads", help="Maximum number of threads to use", type=int, default=3)
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores. (default: 1)")
    parser.add_argument("-n", "--dry-run", action='store_true', default=False, help="Only dry-run the workflow (default False)")
    parser.add_argument('--dag', '-D', action='store_true', default=False, help='Generate a DAG image of the worfklow ("dag.svg")')
    parser.add_argument("-F", "--forceall", action='store_true', default=False, help="Force all output files to be re-created (default False)")
    parser.add_argument("-D, --dry-run", action='store_true', default=False, help="Flag used to show the docker command line. Must be used in conjonction with --docker")

    container_group = parser.add_mutually_exclusive_group()
    container_group.add_argument("--docker", action='store_true', default=False, help="Flag used to run computations on a docker container")
    container_group.add_argument("--singularity", action='store_true', default=False, help="Flag used to run computations on a singularity container")

    return parser


def get_provided_args(parser: argparse.ArgumentParser, args: argparse.Namespace, ignore_args: List=[]):
    provided_args = {}
    for action in parser._actions:
        if any(arg in sys.argv for arg in action.option_strings if arg not in ignore_args):
            provided_args[action.dest] = getattr(args, action.dest)
            
    return provided_args


def validate_required_args(config, required_args):
    missing_args = []
    for arg in required_args:
        parsed_arg = arg
        while str(parsed_arg).startswith("-"):
            parsed_arg = parsed_arg[1:]
        
        if not config[parsed_arg]:
            missing_args.append(arg)

    if missing_args:
        get_parser().print_usage()
        print(f"\nMissing required arguments: {', '.join(missing_args)}")
        sys.exit(1)

