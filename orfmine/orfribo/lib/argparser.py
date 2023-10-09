import argparse
import pkg_resources
import sys
from typing import List

from yaml import safe_load as yaml_safe_load


def get_args() -> argparse.Namespace:
    """Return command line parameters

    Returns:
        argparse.Namespace: command line parameters
    """
    parser = get_parser()
    args = parser.parse_args()

    if args.config:
        args.config = { x.split("=")[0]:x.split("=")[1] for x in args.config }

    args.ram = args.ram * 1000
    
    return args


def get_parser():
    parser = argparse.ArgumentParser(description='Run ORFribo pipeline.')

    parser.add_argument("-c", "--config", help="Path to a YAML configuration file for ORFribo.")
    parser.add_argument("--fna", help="Path to the genome/transcriptome fasta file")
    parser.add_argument("--gff", help="Path to the GFF annotation file")
    parser.add_argument("--gff-intergenic", help="Path to the GFF annotation mapping file")
    parser.add_argument("--fastq", help="Path to the directory containing .fastq.gz files")
    parser.add_argument("--project-name", help="Name for the experiment", default="")
    parser.add_argument("--out", type=str, help="Base directory location for orfribo outputs", default=".")
    parser.add_argument("--rna-to-exclude", "-X", help="Path to a fasta file with nucleotide sequences to exclude", default="")
    parser.add_argument("--adapter", help="Adapter sequence (e.g. 'AGATCGGAAGAGCACACGTCT')", type=str, default="")
    parser.add_argument("--min-read-length", help="Minimum read length for ribosome profiling", type=int, default=25)
    parser.add_argument("--max-read-length", help="Maximum read length for ribosome profiling", type=int, default=35)
    parser.add_argument("--gff-feature", help="Feature element to select during counting", default="CDS")
    parser.add_argument("--gff-attribute", help="Attribute containing gene names in GFF", default="Name")
    parser.add_argument("--mean-threshold", help="Minimum mean of in-frame reads", type=float, default=70)
    parser.add_argument("--median-threshold", help="Minimum median of in-frame reads", type=float, default=70)
    parser.add_argument("--intergenic-features", help="List of features in the intergenic GFF", nargs="*", default=["nc_intergenic"])
    parser.add_argument("--ram", help="Maximum allowed RAM to use (Mb). Defaults to 2000).", type=int, default=2000)
    parser.add_argument("--cores", help="Number of provided cores. Defaults to 1.", type=int, default=1)
    parser.add_argument("--threads", help="Maximum number of threads to use", type=int, default=3)
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores. (default: 1)")
    parser.add_argument("-P", "--preview", action='store_true', default=False, help="Only dry-run the workflow (default False)")
    parser.add_argument('--dag', action='store_true', default=False, help='Generate a DAG image of the worfklow ("dag.svg")')
    parser.add_argument("-F", "--forceall", action='store_true', default=False, help="Force all output files to be re-created (default False)")
    parser.add_argument("--debug", action="store_true", default=False, help="Allow to debug rules with e.g. PDB. This flag allows to set breakpoints in run blocks.")

    parser.add_argument("-D", "--dry-run", action='store_true', default=False, help="Flag used to show the docker command line. Must be used in conjonction with '--docker' or '--singularity'")

    trim_group = parser.add_mutually_exclusive_group(required=False)
    trim_group.add_argument("--trimmed", action='store_true', help="Flag indicating that the sequence adapters are already removed.")
    trim_group.add_argument("--not-trimmed", action='store_true', help="Flag indicating that the sequence adapters are not removed.")

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
        
        if parsed_arg not in config:
            missing_args.append(arg)

    if missing_args:
        get_parser().print_usage()
        print(f"\nMissing required arguments: {', '.join(missing_args)}")
        sys.exit(1)


def get_default_config():
    """Load the default snakefile configuration

    Returns:
        dict: default snakefile configuration
    """
    # 
    default_config_path = pkg_resources.resource_filename("orfmine.orfribo", 'config.yaml')
    default_config = None
    try:
        with open(default_config_path, "r") as f:
                default_config = yaml_safe_load(f)
    except:
        print("Error occured while trying to load the default snakefile configuration")

    return default_config


def update_config_from_file(default_config: dict, configfile: str):
    """ Update the default config from a yaml file

    Args:
        default_config (dict): dict to update
        configfile (str): yaml file containing key:value pairs to update
    """
     
    with open(configfile, "r") as f:
        custom_config = yaml_safe_load(f)
        
    # ensure that only valid keys will be used for the update
    for key in default_config:
         if key not in custom_config:
              print(f"Warning, key {key} in {configfile} is not allowed. It will not be considered.")
              _ = custom_config.pop(key)

    default_config.update(custom_config)


def update_config_from_args(args: argparse.Namespace, config: dict):

    # get provided arguments into a dictionary
    provided_args = get_provided_args(parser=get_parser(), args=args)

    # set --trimmed args value. "--not-trimmed" flag is used as a proxy to set --trimmed value to False if --not-trimmed is in provided_args, True otherwise
    if "not_trimmed" in provided_args:
        is_trimmed = not provided_args.pop("not_trimmed")
        provided_args["trimmed"] = is_trimmed     

    # update default config with given command line args
    for key, value in provided_args.items():
        config[key] = value

    return config


def load_config(args: argparse.Namespace):
    # load default yaml config file
    config = get_default_config()

    # update default config with config file, if given
    if args.config:
        update_config_from_file(default_config=config, configfile=args.config)

    # update default config with given command line args
    config = update_config_from_args(args=args, config=config)

    # check that mandatory arguments are given and valid
    required_args = ["--fna", "--gff", "--gff_intergenic", "--fastq", "--trimmed"]
    validate_required_args(config, required_args)

    return config
