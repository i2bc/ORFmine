# https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt

from pathlib import Path
import re
from typing import List, Union
import requests
import pkg_resources

list

class CodonTables():
    """
    """
    def __init__(self, gc_file: Union[str, Path]):
        self.tables: List[CodonTable] = []
        self._parse_ncbi_gc_table(filename=gc_file)

    def _parse_ncbi_gc_table(self, filename: Union[str, Path]):
        START_KEYWORD = "Genetic-code-table"
        START_TABLE_BLOCK = "{"
        END_TABLE_BLOCK = "}"
        ALLOWED_KEY = ["name", "id", "ncbieaa", "sncbieaa", "--"]

        last_read_key = ""

        start_flag = False
        with open(filename, "r", encoding="utf-8") as gc_file:
            for line in gc_file:
                # skip starting comment lines
                if line.startswith("--"):
                    continue
                # mark start flag as ready if the keyword is met
                if line.startswith(START_KEYWORD):
                    start_flag = True
                    continue
                elif not start_flag:
                    continue

                infos = line.split()
                key = infos[0]

                # instantiate new table dict at new table block
                if key == START_TABLE_BLOCK:
                    gc_table = {}
                    name_idx = 0
                # instantiate a CodonTable at end table block and add it to self.tables
                elif re.search(re.escape(END_TABLE_BLOCK), line):
                    self.tables.append(CodonTable(**gc_table))

                elif key in ALLOWED_KEY:
                    if key == "name":
                        key = f"{key}_{name_idx}"
                        name_idx += 1
                    elif key == "--":
                        key = f"{infos[1].lower()}"
                        
                    if key not in gc_table:
                        content = self._parse_info(infos[1:])
                        gc_table[key] = content
                        last_read_key = key
                else:
                    if last_read_key == "name":
                        gc_table[last_read_key] += infos[:]

    def _parse_info(self, info: List[str]):
        joined_info = ' '.join(info)
        match = re.search(r"Base[0-9]", joined_info)
        if match:
            joined_info = joined_info.replace(match.group(), '')

        return joined_info.replace('"', '').replace("'", "").strip(",").strip()

    def __str__(self) -> str:
        return f"CodonTables: {len(self.tables)} available."

    def __repr__(self) -> str:
        return f"CodonTables: {len(self.tables)} available."

    def get_ids(self):
        return [ x.id for x in self.tables]

    def get_names(self):
        return [ x.name for x in self.tables]

    def get_alt_names(self):
        return [ x.alt_name for x in self.tables]

    def get_infos(self, id: int=1):
        return [ x for x in self.tables ]

    def get_table(self, id: int=1):
        return [ x.table for x in self.tables if x.id == id ][0]

    def codon_table(self, id: int=1, **kwargs):
        return [ x for x in self.tables if x.id == id ][0]

    # Creating an iterator object
    def __iter__(self):
        self.index = 0
        return self

    # Move to next element using __next__
    def __next__(self):
        # Storing the current i value
        index = self.index

        # Stop iteration if limit is reached
        if index > len(self.tables) - 1:
            raise StopIteration

        # double & return old value
        self.index = index + 1 
        return self.tables[index]

class CodonTable():
    """
    Object used to handle a codon table.

    The instantiation requires the unpacking of key:value pairs
    such as the one generated from parse_ncbi_gc_table():

    gc_table = {
        'name_0': 'Standard',
        'name_1': 'SGC0', 
        'id': '1',
        'ncbieaa': 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG',
        'sncbieaa': '---M------**--*----M---------------M----------------------------',
        'base1': 'TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG',
        'base2': 'TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG',
        'base3': 'TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG'
    }

    codonTable = CodonTable(**gc_table)
    """
    def __init__(self, *args, **kwargs) -> None:
        self.name: str = kwargs["name_0"]
        self.id: int = int(kwargs["id"])
        self.alt_name: str = f"SGC{self.id-1}"
        self.table, self.alt_table = self.set_tables(kwargs)

    def set_tables(self, kwargs: dict) -> tuple[dict, dict]:
        table = {}
        alt_table = {}

        for n, aa in enumerate(kwargs["ncbieaa"]):
            codon = kwargs["base1"][n] + kwargs["base2"][n] + kwargs["base3"][n]
            table[codon] = aa

            if kwargs["sncbieaa"][n] == "-":
                alt_table[codon] = aa
            else:
                alt_table[codon] = kwargs["sncbieaa"][n]

        return table, alt_table

    def __str__(self) -> str:
        return f"CodonTable: name: {self.name}, id: {self.id}, name_id: {self.alt_name}"

    def __repr__(self) -> str:
        return f"CodonTable: name: {self.name}, id: {self.id}, name_id: {self.alt_name}"

    def start_codons(self):
        return [ (k,v) for k,v in self.table.items() if v == "M"]

    def stop_codons(self):
        return [ (k,v) for k,v in self.table.items() if v == "*"]


def download_file(url: Union[str, Path], outname: Union[str, Path]):
    with requests.get(url=str(url)) as response:
        response.raise_for_status()

        with open(outname, 'wb') as _f:
            for chunk in response.iter_content(50*1024*1024):
                _f.write(chunk)
