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
from pathlib import Path
import pkg_resources
import snakemake
import time
from yaml import safe_load as yaml_safe_load

from orfmine.orfribo.lib import argparser


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


def update_config(default_config: dict, configfile: str):
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


def generate_dag_svg(snakefile, output_svg_path):
    import subprocess

    cmd = (f"snakemake -s {snakefile} -j --dag -np --forceall --nolock | dot -Tsvg > {output_svg_path}")
    subprocess.run(cmd, shell=True, check=True)



def main():
    args = argparser.get_args()

    config = get_default_config()
    if args.config:
         update_config(default_config=config, configfile=args.config)

    provided_args = argparser.get_provided_args(parser=argparser.get_parser(), args=args)
    for key, value in provided_args.items():
        config[key] = value
         
    required_args = ["--fna", "--gff", "--gff_intergenic", "--path_to_fastq", "--already_trimmed"]
    argparser.validate_required_args(config, required_args)

    # get the orfribo snakefile 
    snakefile = pkg_resources.resource_filename("packages.orfribo", 'Snakefile')

    import json
    print(json.dumps(config, indent=2))
    # exit()

    if args.dag:
        generate_dag_svg(snakefile=snakefile, output_svg_path="orfribo_dag.svg")
        exit()

    # set root directory of orfribo results
    if not config["out_base"]:
        suffix_date = time.strftime("%Y%m%d-%H%M%S") 
        config["out_base"] = Path('orfribo_' + suffix_date)

    if not Path(config["fasta_outRNA"]).exists():
        with open(Path(config["fasta_outRNA"]), "x") as _f:
            pass

    snakemake.snakemake(
        snakefile=snakefile,
        dryrun=args.dry_run,
        nodes=args.jobs,
        resources=args.mem_mb,
        forceall=args.forceall,
        printshellcmds=True,
        config=config,
        force_incomplete=True,
        # omit_from="select_read_lengths"
    )


if __name__ == "__main__":
    # orfribo --fna database/Scer.fna --gff database/Scer.gff --gff_intergenic database/mapping_orf_Scer.gff --path_to_fastq fastq/ --already_trimmed yes --fasta_outRNA database/Scer_rRNA.fa -j 4
    main()
