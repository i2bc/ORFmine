"""

# If the user does not provide a minimum number of cpu available :
# a quarter of computer resources is assigned for the analysis (doubled on specific steps)
used_memory=""
if [ $# -eq 0 ];
then
    cpu=`nproc --all`;
    echo "Total CPU available : "${cpu};
    cpu_use=$((${cpu}/4));
    available_memory_mb=""
else
    cpu_use=$(($1/2));
    available_memory_mb=$(($2*1000))
    echo "Maximum RAM used for the analysis : ${available_memory_mb}Mb";
fi;
echo "Maximum CPU used for the analysis : "$((${cpu_use}*2));

used_memory="--resources mem_mb=${available_memory_mb}"

mkdir /workdir/orfribo/;
# conda list;
snakemake -s /ORFmine/orfribo/RiboDoc_BAM2Reads/Snakefile -j --dag -np --directory /workdir/orfribo/ --nolock | dot -Tsvg > /workdir/orfribo/dag_last-run.svg;
snakemake -s /ORFmine/orfribo/RiboDoc_BAM2Reads/Snakefile -j --dag -np --forceall --directory /workdir/orfribo/ --nolock | dot -Tsvg > /workdir/orfribo/dag_all.svg;
snakemake -s /ORFmine/orfribo/RiboDoc_BAM2Reads/Snakefile -j ${cpu_use} ${used_memory} --directory /workdir/orfribo/ -k --nolock;

"""
import argparse
import pkg_resources
import multiprocessing

import snakemake


def get_args() -> argparse.Namespace:
    """Return command line parameters

    Returns:
        argparse.Namespace: command line parameters
    """
    
    parser = argparse.ArgumentParser(description='Run ORFribo pipeline.')
    parser._action_groups.pop()  # Remove original optional argument group

    mandatory_arguments = parser.add_argument_group('Mandatory arguments')
    optional_arguments = parser.add_argument_group('Optional arguments')

    mandatory_arguments.add_argument(
        "-c", "--configfile",
        required=True,
        nargs="?",
        help="Configuration file defining mandatory inputs, output paths and other optional parameters."
    )

    optional_arguments.add_argument(
        "-C", "--config",
        required=False,
        nargs="*",
        metavar="KEY=VALUE",
        default=dict(),
        help="Configuration file defining mandatory inputs, output paths and other optional parameters."
    )

    optional_arguments.add_argument(
        "-j", "--jobs",
        required=False,
        type=int,
        nargs="?",
        default=1,
        help="Use at most N CPU cluster/cloud jobs in parallel. For local execution this is an alias for --cores. (default: 1)"
    )

    optional_arguments.add_argument(
        "-m", "--mem_mb",
        required=False,
        type=int,
        nargs="?",
        default=50,
        help="Maximum allowed RAM to be used in Mb. (default: 50Mb)"
    )

    optional_arguments.add_argument(
        "-n", "--dry-run",
        required=False,
        action='store_true',
        default=False,
        help="Only dry-run the workflow (default False)"
    )

    optional_arguments.add_argument(
        "--dag",
        required=False,
        action='store_true',
        default=False,
        help="Print the dag in the graphviz dot language (default False)"
    )

    optional_arguments.add_argument(
        "-F", "--forceall",
        required=False,
        action='store_true',
        default=False,
        help="Force all output files to be re-created (default False)"
    )

    args = parser.parse_args()

    if args.config:
        args.config = { x.split("=")[0]:x.split("=")[1] for x in args.config }

    args.mem_mb = {"mem_mb": args.mem_mb * 1000} 
    
    return args


def main():
    args = get_args()

    snakefile = pkg_resources.resource_filename("orfribo", 'Snakefile')

    snakemake.snakemake(
        snakefile,
        # configfiles=args.configfile,
        dryrun=args.dry_run,
        nodes=args.jobs,
        resources=args.mem_mb,
        forceall=args.forceall,
        printdag=args.dag, 
        config=args.config
    )


if __name__ == "__main__":
    main()