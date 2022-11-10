################################################################################
#######################           RiboDoc           ############################
################################################################################

################################################################################
### Tool for Ribosome Profiling
### GST team
################################################################################

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
