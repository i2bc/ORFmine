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
from argparse import Namespace
import json
from pathlib import Path
import pkg_resources
import snakemake
import time

from orfmine import DOCKER_IMAGE
from orfmine.utilities.container import ContainerCLI
from orfmine.orfribo.lib import argparser 

def generate_dag_svg(snakefile, output_svg_path):
    import subprocess

    cmd = (f"snakemake -s {snakefile} -j --dag -np --forceall --nolock | dot -Tsvg > {output_svg_path}")
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Command '{cmd}' failed with error:\n{e.stderr}")
    except Exception as e:
        raise RuntimeError(f"Error executing the command: {cmd}.\nError: {str(e)}")


def set_outdir(config: dict, args: Namespace):
    """Set root directory of orfribo results; it must be created here for container usage

    Args:
        config (dict): preset config dictionary
        args (Namespace): argparse.Namespace instance
    """
    # if outdir already defined, do nothing
    if Path(config["out"]).stem:
        return 

    # else generate an oudtir name
    suffix_date = time.strftime("%Y%m%d-%H%M%S") 
    outdir = f"orfribo_{suffix_date}"
    config["out"] = outdir

    # set args.out value as the generated outdir name
    args.out = outdir
    Path(outdir).mkdir(parents=True, exist_ok=True)


def start_orfribo(args: Namespace, config: dict):
    # get the orfribo snakefile 
    snakefile = pkg_resources.resource_filename("orfmine.orfribo", 'Snakefile')

    if args.dag:
        generate_dag_svg(snakefile=snakefile, output_svg_path="orfribo_dag.svg")
        exit()

    resources = {"mem_mb": args.ram}

    snakemake.snakemake(
        snakefile=snakefile,
        dryrun=args.preview,
        nodes=args.jobs,
        resources=resources,
        forceall=args.forceall,
        printshellcmds=True,
        config=config,
        force_incomplete=True,
        cores=args.cores,
        debug=args.debug,
        # omit_from="find_adapter_sequence"
    )


def run_orfribo_containerized(args: Namespace):
    # load config file
    #config = argparser.load_config(args=args)

    # list of flags related to input files
    input_args = ["--fna", "--gff", "--gff-intergenic", "--fastq"]
    if args.config:
        input_args += ["--config"]


    # update default config with from config file, if given; add input file if present
    if args.config:
        input_args += ["config"]

    # add input file if present
    if args.rna_to_exclude:
        input_args += ["rna_to_exclude"]

    # set root directory of orfribo results
    set_outdir(config=config, args=args)

    # instantiate containerCLI handler
    cli = ContainerCLI(
            input_args=input_args,
            output_arg="--out",
            args=args,
            workdir="/output",
            image_base=DOCKER_IMAGE,
            prog="orfribo",
            container_type="docker" if args.docker else "singularity",
            dev_mode=args.dev,
            package_binding={"orfmine": "/home/orfuser/orfmine/orfmine"}
        )

    cli.show()
    if not args.dry_run:
        cli.run()


def run_orfribo_locally(args: Namespace):
    # load config file. Config sequence setting: default config.yaml -> optional given config file -> provided arguments
    config = argparser.load_config(args=args)
    # exit()
    #args = get_args()

    # set root directory of orfribo results
    set_outdir(config=config, args=args)

    # if not exist, create empty file of ribosomic RNAs to exclude
    #if not Path(config["rna_to_exclude"]).exists():
    #    with open(Path(config["rna_to_exclude"]), "x") as _f:
    #        pass

    # print config
    print(json.dumps(config, indent=2))

    # start orfribo
    start_orfribo(args=args, config=config)



def main():
    args = argparser.get_args()
    if args.docker or args.singularity:
        run_orfribo_containerized(args=args)
    else:
        run_orfribo_locally(args=args)


if __name__ == "__main__":
    # orfribo --fna database/Scer.fna --gff database/Scer.gff --gff_intergenic database/mapping_orf_Scer.gff --fastq fastq/ --already_trimmed yes --rna_to_exclude database/Scer_rRNA.fa -j 4
    main()
