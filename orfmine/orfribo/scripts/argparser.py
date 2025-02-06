import argparse
import pkg_resources
import sys
from typing import List
from yaml import safe_load as yaml_safe_load

from orfmine.utilities.container import add_container_args


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
    parser.add_argument("--project-name", help="Name for the experiment (Defaults '') ", default="")
    parser.add_argument("--out", type=str, help="Base directory location for orfribo outputs (Defaults '.')", default=".")
    parser.add_argument("--rna-to-exclude", "-X", help="Path to a fasta file with nucleotide sequences to exclude (Defaults '')", default="")
    parser.add_argument("--adapter", help="Adapter sequence (e.g. 'AGATCGGAAGAGCACACGTCT') (Defaults: '')", type=str, default="")
    parser.add_argument("--aligner", help="Choose your alignement tool : star or hisat2, (Defaults : 'hisat2')", type=str, default="hisat2")
    parser.add_argument("--min-read-length", help="Minimum read length for ribosome profiling (Defaults : 25)", type=int, default=25)
    parser.add_argument("--max-read-length", help="Maximum read length for ribosome profiling (Defaults: 35)", type=int, default=35)
    parser.add_argument("--gff-feature", help="Feature element to select during counting (Default : 'CDS')", default="CDS")
    parser.add_argument("--mean-threshold", help="Minimum mean of in-frame reads (Default : '70')", type=float, default=70)
    parser.add_argument("--median-threshold", help="Minimum median of in-frame reads (Default : '70')", type=float, default=70)
    parser.add_argument("--intergenic-features", help="List of features in the intergenic GFF (Default : 'nc_intergenic')", nargs="*", default=["nc_intergenic"])
    parser.add_argument("--multi_alignement", type=int, help="The maximum number of allowed multiple alignments for each read (Defaults : 10)", default=10)
    parser.add_argument("--introns_length", type=int, help="Intron length, applicable only when using STAR as the aligner (Defaults: 3000)", default=3000)
    parser.add_argument("--ram", help="Maximum allowed RAM to use (Mb). Defaults to 2000).", type=int, default=2000)
    parser.add_argument("--cores", help="Number of provided cores. Defaults to 1.", type=int, default=1)
    parser.add_argument("--threads", help="Maximum number of threads to use", type=int, default=3)
    parser.add_argument("-j", "--jobs", type=int, default=1, help="Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores. (default: 1)")
    parser.add_argument("-P", "--preview", action='store_true', default=False, help="Only dry-run the workflow (default False)")
    parser.add_argument('--dag', action='store_true', default=False, help='Generate a DAG image of the worfklow ("dag.svg")')
    parser.add_argument("-F", "--forceall", action='store_true', default=False, help="Force all output files to be re-created (default False)")
    parser.add_argument("--debug", action="store_true", default=False, help="Allow to debug rules with e.g. PDB. This flag allows to set breakpoints in run blocks.")

    trim_group = parser.add_mutually_exclusive_group(required=False)
    trim_group.add_argument("--trimmed", action='store_true', help="Flag indicating that the sequence adapters are already removed.")
    trim_group.add_argument("--not-trimmed", action='store_true', help="Flag indicating that the sequence adapters are not removed.")

    parser = add_container_args(parser=parser)

    return parser


def get_provided_args(parser: argparse.ArgumentParser, args: argparse.Namespace, ignore_args: List=[]):
    provided_args = {}
    for action in parser._actions:
        if any(arg in sys.argv for arg in action.option_strings if arg not in ignore_args):
            provided_args[action.dest] = getattr(args, action.dest)
            
    return provided_args


def check_provided_args(provided_args: dict, required_args: List[str], mutually_exclusive: List[tuple]):
    if not provided_args:
        get_parser().print_usage()
        print(f"\nMissing required arguments: {', '.join(required_args)}")
        exit(1)

    missing_args = []
    for r_arg in required_args:            
        if r_arg.lstrip("-").replace("-", "_") not in provided_args:
            missing_args.append(r_arg)

    inconsistencies = []
    if mutually_exclusive:
        for m_excl in mutually_exclusive:
            a, b = m_excl
            if a.lstrip("-").replace("-", "_") not in provided_args and b.lstrip("-").replace("-", "_") not in provided_args:
                missing_args.append(f"{a}/{b}")
            elif a.lstrip("-") in provided_args and b.lstrip("-") in provided_args:
                inconsistencies.append(f"{a}/{b}")


    if missing_args:
        get_parser().print_usage()
        print(f"\nMissing required arguments: {', '.join(missing_args)}")
    if inconsistencies:
        get_parser().print_usage()
        print(f"\Mutually exclusive arguments: {', '.join(inconsistencies)}")

    if missing_args or inconsistencies:
        exit(1)


def validate_required_args(config, required_args):
    missing_args = []
    for arg in required_args:        
        if arg.lstrip("-").replace("-", "_") not in config:
            missing_args.append(arg)

    if missing_args:
        get_parser().print_usage()
        print(f"\nMissing required arguments: {', '.join(missing_args)}")
        exit(1)


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


def update_config_from_args(provided_args: dict, config: dict):

    # set --trimmed args value. "--not-trimmed" flag is used as a proxy to set --trimmed value to False if --not-trimmed is in provided_args, True otherwise
    if "not_trimmed" in provided_args:
        is_trimmed = not provided_args.pop("not_trimmed")
        provided_args["trimmed"] = is_trimmed     

    # update default config with given command line args
    for key, value in provided_args.items():
        config[key] = value

    return config


def load_config(args: argparse.Namespace):
    required_args = ["--fna", "--gff", "--gff-intergenic", "--fastq"]
    mutually_exclusive_args = [("--trimmed", "--not-trimmed")]
    # get provided arguments into a dictionary
    provided_args = get_provided_args(parser=get_parser(), args=args)

    # check that provided arguments contains mandatories one
    check_provided_args(provided_args=provided_args, required_args=required_args, mutually_exclusive=mutually_exclusive_args)
    

    # load default yaml config file
    config = get_default_config()

    # update default config with config file, if given
    if args.config:
        update_config_from_file(default_config=config, configfile=args.config)

    # update default config with given command line args
    config = update_config_from_args(provided_args=provided_args, config=config)

    # check that mandatory arguments are given and valid
    validate_required_args(config, required_args+["--trimmed"])

    return config
