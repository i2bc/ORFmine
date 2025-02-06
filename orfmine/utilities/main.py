"""
Utility script used to handle the download, update and generation of codon_tables JSON files.

"""
import json
from pathlib import Path
from typing import Union
import pkg_resources

from orfmine.utilities.lib import table_handlers


NCBI_GC_FILENAME = "gc.prt"
NCBI_GC_URL = "https://ftp.ncbi.nih.gov/entrez/misc/data/" + NCBI_GC_FILENAME
OUTPATH = pkg_resources.resource_filename(__name__, 'data')


def download_ncbi_tables_file():
    out_filename = str(Path(OUTPATH) / Path(NCBI_GC_URL).name)
    print(f"Downloading {str(NCBI_GC_URL)} into {str(OUTPATH)} ...")
    table_handlers.download_file(url=NCBI_GC_URL, outname=out_filename)
    print("Download process done.\n")


def dump_tables(gc_filename: Union[str, Path]):
    if not Path(gc_filename).exists():
        print(f"{gc_filename} not found. Process aborted.\n")
        exit(1)
    
    table_basename = str(Path(gc_filename).parent / "codon_table_{}.json")    
    codon_tables = table_handlers.CodonTables(gc_file=gc_filename)

    for element in codon_tables:
        out_table = table_basename.format(element.id)
        print(f"Writing codon table {element.id} in {out_table} ...")
        with open(out_table, 'w') as outjson:
            outjson.write(json.dumps(element.table, indent=2))
    print("Writing codon table JSON files done.\n")


def update_ncbi_tables_file(gc_filename: Union[str, Path]):
    download_ncbi_tables_file()
    dump_tables(gc_filename=gc_filename)


if __name__ == "__main__":
    update_ncbi_tables_file(gc_filename=str(Path(OUTPATH) / NCBI_GC_FILENAME))