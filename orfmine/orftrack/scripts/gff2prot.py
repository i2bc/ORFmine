import argparse
from pathlib import Path
import re
import sys
from typing import List
import time
from os import stat as ostat
from orfmine import DOCKER_IMAGE

from orfmine.orftrack.lib import gff_parser, fasta_parser
from orfmine.orftrack.lib.utils import CDSQueue, get_chunks, index_genomes, validate_chromosomes, merge_outfiles

import multiprocessing

from orfmine.utilities.container import ContainerCLI
from orfmine.utilities.container import add_container_args


def get_args():
    """
    Returns:
        Parameters
    """

    parser = argparse.ArgumentParser(description='Returns amino acid and nucleic fasta files from genomic data')
    group = parser.add_mutually_exclusive_group()

    parser.add_argument(
        "--fna",
        required=True,
        nargs="?",
        help="Genomic fasta file (.fna)"
    )
    
    parser.add_argument(
        "--gff",
        required=True,
        nargs="?",
        help="GFF annotation file (.gff)"
    )

    parser.add_argument(
        "--chromosomes",
        required=False,
        nargs="*",
        type=str,
        default=[],
        help="Chromosome names to process. By default, all chromosomes are processed (default: [])."
    )

    parser.add_argument(
        "-O", "--outdir",
        required=False,
        nargs="?",
        default='./',
        type=str,
        help="Output directory"
    )

    parser.add_argument(
        "-B", "--out-basename",
        required=False,
        nargs="?",
        default="",
        type=str,
        help="Basename for output files (default: '_proteins' suffix added to fasta file name)"
    )

    parser.add_argument(
        "--codon-table",
        required=False,
        type=str,
        default="1",
        help="Codon table ID to be used for considered chromosomes. \
            The ID must be consistent with the NCBI codon table IDs \
            (see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). \
            By default, the standard codon table (id = 1) is used. Defaults to '1'"
    )
    
    parser.add_argument(
            "--features",
            type=str,
            required=True, 
            nargs="*",
            help="Annotation features to be considered (By definition is all)"
    )

    group.add_argument(
        "-P", "--proteic",
        required=False,
        action='store_true',
        default=False,
        help="For only amino acids output format (.pfasta)."
    )

    group.add_argument(
        "-N", "--nucleic",
        required=False,
        action='store_true',
        default=False,
        help="For nucleic output format (.nfasta)."
    )

    parser.add_argument(
        "--elongate",
        type=int,
        required=False, 
        nargs="?",
        default=None,
        help="Number of nucleotides the sequence of matching GFF elements must be elongated from both extremities. Defaults to None"
    )

    parser.add_argument(
        "--chr-exclude",
        type=str,
        required=False, 
        nargs="*",
        default=[],
        help="List of chromosomes to be excluded. Defaults to None"
    )

    parser.add_argument(
        "--cpus",
        type=int,
        required=False, 
        nargs="?",
        default=multiprocessing.cpu_count(),
        help="Maximum number of CPUs to use. By default, all available CPUs are used."
    )

    parser.add_argument(
        "-E", "--stop-end",
        required=False,
        action='store_true',
        default=False,
        help="Quality flag used to process only proteins with a STOP codon at the end of their sequence. Defauts to False."
    )

    parser.add_argument(
        "-I", "--no-stop-intra",
        required=False,
        action='store_true',
        default=False,
        help="Quality flag used to process only proteins that does not possess a STOP codon inside their sequence. Defauts to False."
    )

    parser = add_container_args(parser=parser)


    args = parser.parse_args()

    return args


def process_chromosome(
        gff_filename: str=None,
        fasta_chr: fasta_parser.Fasta=None,
        chr_id: str="",
        features: list=[],
        basename_out: str="",
        out_formats: list=[".pfasta", ".nfasta"],
        elongation: int=None,
        gff_indexes: dict={},
        check_stop_end: bool=False,
        check_stop_intra: bool=False,
        ) -> dict:
    """
    Return a dictionary with chromosome name as key, and gff_parser.Chromosome instance as value
    if asked GFF elements that are CDS only
    and / or 
    write on the fly fasta / gff file for asked GFF elements that are not CDS.

    Args:
        gff_filename (str, optional): GFF filename. Defaults to None.
        fasta_chr (str, optional): fasta_parser.Fasta instance. Defaults to None.
        asked_chrs (list, optional): Chromosome names to be processed. By default, all chromosomes are processed (default: []).
        features (list, optional): List of GFF features to process. Empty list means all features will be treated. Defaults to [].
        basename_out (str, optional): Output file basename. Defaults to "".
        out_formats (list, optional): List of output formats to be generated. Outputs can be .fna, .pfasta or both.
        asked_chrs (list, optional): Chromosome names not to be processed. By default, all chromosomes are processed (default: []).

    Returns:
        dict: Dictionary with chromosome name as key, and gff_parser.Chromosome instance as value
              with CDS GFF elements only.  
    """
    start_time = time.time()

    process_name = multiprocessing.current_process().name
    print(f"{process_name} -  Starting to process chromosome {chr_id}")

    features_to_keep = "|".join(features)

    # queue used to process CDS particular (CDS can't be written 'on the fly')
    queue = CDSQueue(check_stop_end=check_stop_end, check_stop_intra=check_stop_intra)
    
    try:
        out_files = { _ext:open(basename_out+_ext, "w") for _ext in out_formats }
        has_elements = False

        with open(gff_filename, 'r') as gff_file:
            # retrieve end-of-file value
            eof = gff_file.seek(0, 2)

            # place cursor to the position of the 1st given chromosome line 
            gff_file.seek(gff_indexes[chr_id], 0)

            # read line as long as it is not starting with '#' and chr name = chr_id
            while True:
                line = gff_file.readline()
                # skip lines starting with a comment mark
                if line.startswith("#"):
                    continue

                # quit the loop if chromosome name is not the chromosome of interest
                chr_name = line.split('\t')[0]
                if chr_name != chr_id:
                    break

                element_type = line.split('\t')[2]
                
                # skip all gff elements not matching features_to_keep:
                if not re.findall(features_to_keep, element_type):
                    continue

                has_elements = True

                try:
                    gff_element = gff_parser.GffElement(gff_line=line, fasta_chr=fasta_chr)
                except KeyError as e:
                    print(f"Warning: {chr_name} seems to be absent from the fasta file. Continue...")
                

                # process non CDS elements
                if gff_element.type != "CDS":
                    if elongation:
                        # warning: element coords in gff are not elongated coords but they are in the fasta files
                        if ".gff" in out_files:
                            out_files[".gff"].write(gff_element.get_gffline())
                        gff_element.start = gff_element.start - elongation
                        gff_element.end = gff_element.end + elongation
                    
                    # writing non CDS elements
                    if ".pfasta" in out_files:
                        out_files[".pfasta"].write(gff_element.get_fastaline())
                    if ".nfasta" in out_files:
                        out_files[".nfasta"].write(gff_element.get_fastanuc_line())

                # CDS elements must be treated differently (CDS must be merged for multi-exonic proteins)
                else:
                    queue.add(cds=gff_element)
                    if queue.cds_completed:
                        if elongation:
                            if ".gff" in out_files:
                                out_files[".gff"].write(queue.build_fake_gff(elongate=elongation))
                            queue.elongate(value=elongation)

                        # writing CDS elements
                        if ".pfasta" in out_files:
                            out_files[".pfasta"].write(queue.get_fasta(_type="proteic"))
                        if ".nfasta" in out_files:
                            out_files[".nfasta"].write(queue.get_fasta(_type="nucleic"))

                if gff_file.tell() == eof:
                    break
        
        # The last CDS is still in cds_list
        if queue.cds_list:            
            # must copy cds_list in cds_completed since queue.get_fasta() applies only on cds_completed
            queue.cds_completed = queue.cds_list
            if elongation:
                if ".gff" in out_files:
                    out_files[".gff"].write(queue.build_fake_gff(elongate=elongation))                
                queue.elongate(value=elongation)

            # writing CDS elements
            if ".pfasta" in out_files:
                out_files[".pfasta"].write(queue.get_fasta(_type="proteic"))
            if ".nfasta" in out_files:
                out_files[".nfasta"].write(queue.get_fasta(_type="nucleic"))
    finally:
        # close opened files
        for _, _file in out_files.items():
            _file.close()

        # remove empty files
        if not has_elements:
            print(f"No asked features found for chromosome {chr_id}")
            for _ext in out_formats:
                if ostat(basename_out + _ext).st_size == 0:
                    Path(basename_out + _ext).unlink()            
        

    print(f"{process_name} -  Chromosome {chr_id} successfully processed in {round(time.time()-start_time, 2)} seconds")
                                

def run_gff2prot(args):

    genomic_fna = args.fna
    genomic_gff = args.gff
    chromosomes_to_process = args.chromosomes
    chromosomes_to_exclude = args.chr_exclude
    features = args.features
    elongation = args.elongate
    cpus = args.cpus
    codon_table_id = args.codon_table
    check_stop_end = args.stop_end
    check_stop_intra = args.no_stop_intra

    out_formats = [".pfasta", ".nfasta"]
    if args.proteic:
        out_formats.remove(".nfasta")
    elif args.nucleic:
        out_formats.remove(".pfasta")


    outpath = Path(args.outdir)
    outpath.mkdir(parents=True, exist_ok=True)

    features_in_name = "_".join(features)
    if not args.out_basename:
        basename_out = str(outpath / f"{Path(genomic_gff).stem}_{features_in_name}")
    else:
        basename_out = str(outpath / f"{Path(args.out_basename).stem}")

    if elongation:
        out_formats.append(".gff")

    # get useful indexes of fasta & gff files 
    fasta_hash, gff_indexes = index_genomes(fasta_file=genomic_fna, gff_file=genomic_gff, codon_table_id=codon_table_id)

    # chromosomes must be consistent between the fasta file and gff file, and the asked chromosomes
    valid_chromosomes = validate_chromosomes(to_include=chromosomes_to_process, to_exclude=chromosomes_to_exclude, in_gff=gff_indexes.keys(), in_fasta=fasta_hash.keys())
    print("Chromosomes that will be processed:")
    for chr in valid_chromosomes:
        print(f" - {chr}")
    print()


    chunks = get_chunks(l=valid_chromosomes, chunks_size=cpus)

    i = 0
    merged_files = []
    for n, chunk in enumerate(chunks):
        processes = []
        out_filenames = [] 

        for chr_id in chunk:
            basename = f"{basename_out}_{chr_id}_{i}"
            out_filenames.append(basename)

            p = multiprocessing.Process(
                target=process_chromosome,
                args=(genomic_gff, fasta_hash[chr_id], chr_id, features, basename, out_formats, elongation, gff_indexes, check_stop_end, check_stop_intra)
            )
            p.chr_id = chr_id
            processes.append(p)
            i += 1

        for p in processes:
            p.start()

        for p in processes:
            p.join()

        merged_files.append(outpath/f"out_merged_chunk-{n}")
        merge_outfiles(inpath=outpath, infiles=out_filenames, extensions=out_formats, outbasename=f"out_merged_chunk-{n}")
    merge_outfiles(inpath=outpath, infiles=merged_files, extensions=out_formats, outbasename=str(Path(basename_out).stem))
    

def run_gff2prot_containerized(args):
    # instantiate containerCLI handler
    cli = ContainerCLI(
            input_args=["--fna", "--gff"],
            output_arg="--outdir",
            args=args,
            image_base=DOCKER_IMAGE,
            prog="gff2prot",
            container_type="docker" if args.docker else "singularity",
            dev_mode=args.dev,
            package_binding={"orfmine": "/home/orfuser/orfmine/orfmine"}
        )
    
    cli.show()
    if not args.dry_run:
        cli.run()


def main():
    args = get_args()

    if args.docker or args.singularity:
        run_gff2prot_containerized(args=args)
    else:
        run_gff2prot(args=args)

if __name__ == '__main__':
    main()