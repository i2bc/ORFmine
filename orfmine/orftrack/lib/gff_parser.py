# -*- coding: utf-8 -*-
"""
Created on Sun Jul 12 16:59:28 2020

@author: nicolas
"""
import sys

from orfmine.orftrack.lib import logHandler, inspect
from orfmine.orftrack.lib.parameters import Parameter
from orfmine.utilities.lib.logging import get_logger


logger = get_logger(name=__name__)


class GffElement:

    def __init__(self, gff_line=None, fasta_chr=None):
        self.gff_line = gff_line.split('\t') if gff_line else None
        self.fasta_chr = fasta_chr if fasta_chr else None
        self.len_chr = fasta_chr.nucid_max if fasta_chr else None

        self.seqid = self.gff_line[0] if gff_line else None
        self.source = self.gff_line[1] if gff_line else None
        self.type = self.gff_line[2] if gff_line else None
        self.start = int(self.gff_line[3]) if gff_line else None
        self.end = int(self.gff_line[4]) if gff_line else None
        self.score = self.gff_line[5] if gff_line else '.'
        self.strand = self.gff_line[6] if gff_line else None

        self.idx = None
        self.nb = None

        if gff_line:
            if self.gff_line[7] != '.':
                self.phase = int(self.gff_line[7])
            else:
                self.phase = self.gff_line[7]

            self.frame = self._get_frame()
            self._get_attributes()
        else:
            self.phase = '.'
            self.frame = None
            self.id_ = None
            self.name = None
            self.parent = None
            self.status = None
            self.color = None

        self.ovp_phased = []
        self.ovp_unphased = []
        self.suborfs = []
        self.orf_ovp = None

    def _get_attributes(self):
        attributes = self.gff_line[8:][0]

        if 'ID' in attributes:
            self.id_ = self._parse_attributes(key='ID')
        else:
            self.id_ = '_'.join([self.type, str(self.start), str(self.end), self.strand])
        if 'Name' in attributes:
            self.name = self._parse_attributes(key='Name')
        else:
            self.name = self.id_

        if 'Parent' in attributes:
            self.parent = self._parse_attributes(key='Parent')
        else:
            self.parent = None
        self.status = None
        self.color = None

    def _parse_attributes(self, key):
        attributes_col = self.gff_line[8:][0].split(';')
        attribute = [x for x in attributes_col if key in x][0]

        return attribute.split('=')[-1]

    def _get_frame(self):
        if self.strand == '+':
            return (self.get_coors()[0] - 1) % 3
        else:
            return (self.len_chr - self.get_coors()[1]) % 3

    def get_coors(self):
        return self._coors_adjusted()

    def _coors_adjusted(self):
        start = self.start
        end = self.end
        if isinstance(self.phase, int):
            offset = (3 - ((self.end-self.start+1-self.phase) % 3)) % 3

            if self.strand == '+':
                start = self.start + self.phase
                end = self.end+offset if self.end+offset <= self.len_chr else self.len_chr
            elif self.strand == '-':
                start = self.start-offset if self.start-offset > 0 else self.start
                end = self.end-self.phase

        return start, end

    def get_len(self):
        return self.get_coors()[1] - self.get_coors()[0] + 1

    def sequence(self):
        phase = 0 if not isinstance(self.phase, int) else self.phase

        return self.fasta_chr.sequence(start=self.start, end=self.end, strand=self.strand, phase=phase)

    def translate(self):
        return self.fasta_chr.translate(start=self.start, end=self.end, strand=self.strand, phase=self.phase)

    def get_fastaline(self):
        fastaline = '>'+self.id_+'\n'+self.translate()+'\n'

        return fastaline

    def get_fastanuc_line(self):
        fastaline = '>'+self.id_+'\n'+self.sequence()+'\n'

        return fastaline

    def get_gffline(self):
        # if self.gff_line and 'frag' not in self.type:
        #     return '\t'.join(self.gff_line)
        if not self.status:
            return '\t'.join(self.gff_line)
        else:
            gff_line = self.seqid
            gff_line += '\t' + self.source
            gff_line += '\t' + self.type
            gff_line += '\t' + str(self.start)
            gff_line += '\t' + str(self.end)
            gff_line += '\t' + self.score
            gff_line += '\t' + self.strand
            gff_line += '\t' + self.phase
            gff_line += '\tID=' + self.format_id()
            gff_line += ';Parent=' + self.parent
            gff_line += ';Status=' + self.status
            gff_line += ';color=' + self.color

            if self.ovp_phased:
                gff_line += ';Ovp_with=' + '|'.join([x.name.strip() for x in self.ovp_phased])
            elif self.ovp_unphased:
                gff_line += ';Ovp_with=' + '|'.join([x.name.strip() for x in self.ovp_unphased])

            return gff_line.strip() + '\n'  #  @BUG added strip() to prevent extra '\n', but dont know yet where it comes from - to identify

    def run_assignment(self, elements, param, is_fragment=False):
        self.set_ovp_elements(elements=elements, co_ovp=param.co_ovp)
        self.set_attributes(is_fragment=is_fragment)

    def set_attributes(self, is_fragment=False):
        self.set_type(is_fragment=is_fragment)
        self.id_ = self.format_id()
        self.set_parent()
        self.set_color()
        self.set_status()

    def set_ovp_elements(self, elements=None, co_ovp=0.7):
        """

        Adds ORF overlapping elements to either self.ovp_phased or self..ovp_unphased lists.

        Args:
            elements (list[GffElement]): elements of the GFF input files
            co_ovp (float): cutoff used to define an overlapping element with an ORF

        Returns:
            None

        """
        if elements:
            # orf_ovp_max = -1
            for element in elements:
                orf_ovp, element_ovp, is_overlap = get_overlap(orf_coors=self.get_coors(), other_coors=element.get_coors())
                if element_ovp == 1.0 or orf_ovp >= co_ovp:
                    if isinstance(element.phase, int) and element.frame == self.frame and element.strand == self.strand:
                        if element not in self.ovp_phased:
                            self.ovp_phased.append(element)
                    else:
                        if element not in self.ovp_unphased:
                            element.orf_ovp = orf_ovp
                            self.ovp_unphased.append(element)

    def set_type(self, is_fragment=False):
        if self.ovp_phased:
            self.type = 'c_CDS'

        elif self.ovp_unphased:
            if 'CDS' in [x.type for x in self.ovp_unphased]:
                ovp_unphased_elements = [x for x in self.ovp_unphased if x.type == 'CDS']
            else:
                ovp_unphased_elements = self.ovp_unphased

            ovp_elements_same_strand = [x for x in ovp_unphased_elements if x.strand == self.strand]
            if ovp_elements_same_strand:
                highest_overlapping_element = sorted([x for x in ovp_elements_same_strand], key=lambda x: x.orf_ovp)[-1]
            else:
                highest_overlapping_element = sorted([x for x in ovp_unphased_elements], key=lambda x: x.orf_ovp)[-1]

            if not is_fragment:
                if highest_overlapping_element.strand == self.strand:
                    self.type = 'nc_ovp_same-' + highest_overlapping_element.type
                else:
                    self.type = 'nc_ovp_opp-' + highest_overlapping_element.type
            else:
                if highest_overlapping_element.strand == self.strand:
                    self.type += '_ovp_same-' + highest_overlapping_element.type
                else:
                    self.type += '_ovp_opp-' + highest_overlapping_element.type
        else:
            if not is_fragment:
                self.type = 'nc_intergenic'

    def format_id(self):
        return '_'.join([self.seqid, self.strand,
                         str(self.start) + '-' + str(self.end),
                         str(self.frame), self.type])

    def set_parent(self):
        self.parent = self.seqid+'_'+'1-'+str(self.len_chr)

    def set_color(self):
        if self.type == 'c_CDS':
            self.color = '#DF7D2E'  # ff4d4d
        else:
            if self.type == 'nc_intergenic':
                self.color = '#5E5E5E'
            else:
                if 'opp' in self.type:
                    self.color = '#009051'
                else:
                    self.color = '#2657AF'

    def set_status(self):
        if self.ovp_phased:
            self.status = 'coding'
        else:
            self.status = 'non-coding'

    def fragment_phased_cds(self, orf_len=60):
        """

        Function that fragments an ORF sequence at the borders of its overlapping CDS in the same frame.

        Given this illustrative configuration:

        *-----------|CDS1|-------|CDS2|------------*


        The fragmentation process will give three shorter ORFs called suborfs:

        nc_5-CDSfrag
        *-----------|
                      nc_intra-CDSfrag
                         |-------|
                                       nc_3-CDSfrag
                                      |------------*

        Note that a fragment is considered a suborf only if its length is at least equal to orf_len.

        """
        cds_ele = None
        cds_elements = sorted(self.ovp_phased, key=lambda x: x.start, reverse=False)

        for i, cds_ele in enumerate(cds_elements):
            if cds_ele.start == self.start:
                continue
            start = self.start if i == 0 else cds_elements[i-1].end + 1
            end = cds_elements[i].start - 1

            if end - start + 1 >= orf_len:
                fragment = GffElement(gff_line=self.get_gffline(), fasta_chr=self.fasta_chr)
                fragment.start = start
                fragment.end = end

                if self.strand == '+':
                    if cds_elements[i].idx == 1:
                        fragment.type = 'nc_5-CDSfrag'
                    else:
                        fragment.type = 'nc_intra-CDSfrag'
                else:
                    if cds_ele.idx == cds_ele.nb:
                        fragment.type = 'nc_3-CDSfrag'
                    else:
                        fragment.type = 'nc_intra-CDSfrag'

                self.suborfs.append(fragment)

        if cds_ele.end != self.end:
            start = cds_ele.end + 1
            end = self.end
            if end - start + 1 >= orf_len:
                fragment = GffElement(gff_line=self.get_gffline(), fasta_chr=self.fasta_chr)
                fragment.start = start
                fragment.end = end

                if cds_ele.idx == cds_ele.nb:
                    fragment.type = 'nc_3-CDSfrag' if self.strand == '+' else 'nc_5-CDSfrag'
                else:
                    fragment.type = 'nc_intra-CDSfrag'

                if self.strand == '+':
                    if cds_ele.idx == cds_ele.nb:
                        fragment.type = 'nc_3-CDSfrag'
                    else:
                        fragment.type = 'nc_intra-CDSfrag'
                else:
                    if cds_ele.idx == 1:
                        fragment.type = 'nc_5-CDSfrag'
                    else:
                        fragment.type = 'nc_intra-CDSfrag'

                self.suborfs.append(fragment)


class Chromosome:

    def __init__(self, id_, fasta_chr):
        self.id_ = id_
        self.fasta_chr = fasta_chr
        self.start = 1
        self.end = fasta_chr.nucid_max
        self.source = None
        self.coors_intervals = self._set_intervals()
        self.gff_elements = []
        
    def _set_intervals(self, value=50000):
        """

        Defines a dictionary containing a set of intervals as keys and an empty list as values.

        For instance, for @param value = 50000:

        {(1, 49999): [],
         (50000, 99999): [],
         ...
         }

        """
        return {(x, x + value - 1): [] for x in range(1, self.end, value)}
        
    def _get_intervals(self, coors):
        """

        Returns a tuple/generator of interval coordinates that overlap with the coordinates given in coors.

        If self.coors_intervals is defined as describe in _set_intervals() and coors=(47500, 52000),
        then the function will return ((1, 49999), (50000, 99999))

        """
        return (x for x in self.coors_intervals if get_overlap(x, coors)[2])
        
    def _add_to_intervals(self):
        """

        Adds the index of the self.gff_elements last element in the correct intervals list of self.coors_intervals.

        For example, if self.coors_intervals is defined as describe in _set_intervals() and
        the coordinates of the last element are (47500, 52000), then the index of the last element will be added as follow:

        {(1, 49999): [..., last_element_idx],
         (50000, 99999): [..., last_element_idx],
         (100000, 149999): [...],
         ...
         }

        """
        last_element = self.gff_elements[-1]
        for interval in self._get_intervals(coors=last_element.get_coors()):
            self.coors_intervals[interval].append(self.gff_elements.index(last_element))
        
    def add(self, gff_element):
        """
        Adds a Gff_element instance in self.gff_elements.
        
        Arguments:
            - gff_element: instance of Gff_element()
        """
        self.gff_elements.append(gff_element)
        self._add_to_intervals()
        
    def get_elements_in_intervals(self, coors):
        """

        Returns all gffElement instances overlapping with coors.

        intervals_mx is a tuple/generator of lists of self.coors_intervals values overlapping with coors. For example,
        if coors=(47500, 52000), then intervals_mx will be: ([..., element_i, element_j, ...], [element_j, element_k, ...])
        where element_i (and others) are indexes of elements in self.gff_elements.

        intervals_flat flattens the matrix, remove redundant elements and sort them. For example, according to the case
        above, intervals_flat will be: (..., element_i, element_j, element_k, ...).


        """
        intervals_mx = (self.coors_intervals[x] for x in self._get_intervals(coors=coors))
        intervals_flat = sorted(set([val for sublist in intervals_mx for val in sublist]))

        return (self.gff_elements[x] for x in intervals_flat)
        
    def get_elements(self, coors=None, frame=None, strand=None, types=None):
        """
        Returns a list of Gff_element instances of a given type. If the frame is
        given, only CDS in this frame will be returned, all CDS otherwise.
        
        Arguments:
            - frame: None or int in 0, 1 or 2
            - strand: None or str in '+' or '-'
            - coors: None or list of coors in the form [104, 395]
            - types: None or list of str (e.g. ['CDS', 'tRNA']
            
        Returns:
            - list of Gff_element instances        
        """

        if types:
            if coors:
                if strand:
                    elements = (x for x in self.get_elements_in_intervals(coors) if x.type in types and x.strand == strand)
                else:
                    elements = (x for x in self.get_elements_in_intervals(coors) if x.type in types)
            else:
                if strand:
                    elements = (x for x in self.gff_elements if x.type in types and x.strand == strand)
                else:
                    elements = (x for x in self.gff_elements if x.type in types)
        else:
            if coors:
                if strand:
                    elements = (x for x in self.get_elements_in_intervals(coors) if x.strand == strand)
                else:
                    elements = (x for x in self.get_elements_in_intervals(coors))
            else:
                if strand:
                    elements = (x for x in self.gff_elements if x.strand == strand)
                else:
                    elements = (x for x in self.gff_elements)

        if not frame:
            return elements
        else:
            return (x for x in elements if x.frame == frame)
        
    def get_types(self):
        return set([x.type for x in self.gff_elements])
        
    def sequence(self, start: int, end: int, strand='+', phase=0):
        start = start if start else self.start
        end = end if end else self.end
        return self.fasta_chr.sequence(start=start, end=end, strand=strand, phase=phase)
                                      
    def rev_comp(self):
        return self.fasta_chr.reverse_complement(self.sequence(start=1, end=10))

    def get_cds(self):
        return (x for x in self.gff_elements if x.type == 'CDS')

    def group_cds(self):
        proteins_dict = {}

        for cds in self.get_cds():
            if cds.name not in proteins_dict:
                proteins_dict[cds.name] = []

            if cds.strand == "-":
                proteins_dict[cds.name].insert(0, cds)
            else:
                proteins_dict[cds.name].append(cds)
    
        return proteins_dict

    def proteins_fasta(self):
        proteins = self.group_cds()
        for protein in sorted(proteins):            
            fasta = '>' + protein + ':' + self.id_ + '\n'
            fasta += ''.join([cds.translate() for cds in proteins[protein]]) + '\n'
            yield fasta

    def proteins_fastanuc(self):
        proteins = self.group_cds()
        for protein in sorted(proteins):
            fastanuc = '>' + protein + ':' + self.id_ + '\n'
            fastanuc += ''.join([cds.sequence() for cds in proteins[protein]]) + '\n'
            yield fastanuc

    def index_cds(self):
        proteins = self.group_cds()
        for protein in sorted(proteins):
            for i, cds in enumerate(proteins[protein], start=1):
                cds.idx = i
                cds.nb = len(proteins[protein])


def get_overlap(orf_coors=(), other_coors=None):
    """
    Function defining if ORF coordinates overlap with another genomic element
    coordinates.
    
    Arguments:
        - orf_coors: start and end coordinates of an ncORF (list)
        - other_coors: start and end coordinates of a genomic element (list)
        
    Returns:
        - a tuple of the overlapping fraction between the ORF and the other element
        (float if overlap, None otherwise)
    """
    is_overlap = False
    orf_ovp, other_ovp = 0, 0
    x_max = max(orf_coors[0], other_coors[0])
    y_min = min(orf_coors[1], other_coors[1])
    if x_max < y_min:
        is_overlap = True
        len_ovp = y_min - x_max + 1
        len_orf = orf_coors[1] - orf_coors[0] + 1
        len_other = other_coors[1] - other_coors[0] + 1
        orf_ovp = len_ovp / float(len_orf)
        other_ovp = len_ovp / float(len_other)
        
    return orf_ovp, other_ovp, is_overlap


GFF_DESCR = {}

def set_gff_descr(gff_fname):
    global GFF_DESCR
    GFF_DESCR = {}
        
    with open(gff_fname, 'rb') as gff_file:
        line = gff_file.readline().decode(encoding='utf-8')
        while line:
            if not line.startswith('#'):
                name = line.split('\t')[0]
                pos_chr = gff_file.tell()-len(line)
                if name not in GFF_DESCR:
                    GFF_DESCR[name] = pos_chr
                
            line = gff_file.readline().decode(encoding='utf-8')


def parse(param: Parameter, fasta_hash: dict):
    """ chr_asked=param.chr, chr_exclude=param.chr_exclude
    @param fasta_hash:
    @param param:
    @type chr_exclude: list
    @type chr_asked: list
    """
    
    gff_fname = param.gff_fname
    fasta_hash = fasta_hash
    chr_asked = param.chr if param.chr else []
    chr_exclude = param.chr_exclude if param.chr_exclude else []

    if not GFF_DESCR:
        set_gff_descr(gff_fname)

    logger.debug('Checking chromosome IDs consistency between GFF and fasta file...')
    # Get chromosomes present both in the GFF file and the fasta file
    chrs_common = inspect.check_chrids(chrs_gff=sorted(GFF_DESCR), chrs_fasta=sorted(fasta_hash))

    # Get chromosomes to be treated
    if not chr_asked:
        chr_ids = sorted([x for x in chrs_common if x not in chr_exclude])
    else:
        chr_ids = sorted([x for x in chr_asked if x not in chr_exclude])

    # If chr_asked, check that asked chromosomes are valid (i.e. present in gff and fasta files)
    if chr_ids == chr_asked:
        invalid_chrs = []
        for _chr in chr_ids:
            if _chr not in chrs_common:
                invalid_chrs.append(_chr)

        if invalid_chrs:
            logger.error('Error: wrong chromosome ID have been asked:')
            for invalid in invalid_chrs:
                logger.error(' - {}'.format(invalid))
                logger.error('')
            sys.exit(1)

    gff_data = {}
    with open(gff_fname, 'r') as gff_file:
        eof = gff_file.seek(0, 2)
        for chr_id in chr_ids:
            gff_file.seek(GFF_DESCR[chr_id], 0)
            line = gff_file.readline()
            chr_name = line.split('\t')[0]
            while chr_name == chr_id:
                if chr_name not in gff_data:
                    logger.debug('  - Reading chromosome: ' + chr_name)
                    gff_data[chr_name] = Chromosome(id_=chr_name, fasta_chr=fasta_hash[chr_id])
                    chromosome = gff_data[chr_name]
                    chromosome.source = line.split('\t')[1]

                element_type = line.split('\t')[2]
                if element_type not in ['chromosome', 'region', 'match']:
                    if not param.types_except and not param.types_only:
                        chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))
                    else:
                        if param.types_except:
                            if element_type not in param.types_except:
                                chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))
                        elif param.types_only:
                            if element_type in param.types_only:
                                chromosome.add(gff_element=GffElement(gff_line=line, fasta_chr=fasta_hash[chr_id]))
                
                # read next line as long as it is not starting with '#'
                while True:
                    line = gff_file.readline()
                    if not line.startswith("#"):
                        break

                if gff_file.tell() == eof:
                    break
                else:
                    chr_name = line.split('\t')[0]

    # Assign an index to CDS to facilitate their fusion into protein
    for chr_name in sorted(gff_data):
        gff_data[chr_name].index_cds()

    return gff_data


def get_parent(gff_line):
    attributes = gff_line.split('\t')[8:][0].split(';')
    parent = [x for x in attributes if 'Parent' in x][0].split('=')[-1]

    return parent
