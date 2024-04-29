# from orfmine.orftrack.lib import logHandler


def get_infos(_input: str, option: str):
    """
    @param _input: gff filename
    @param option: 'chrs' or 'types'
    """
    if option == 'chrs':
        chr_names = []
        with open(_input, 'r') as gff_file:
            for line in gff_file:
                if not line.startswith('#'):
                    chr_name = line.split('\t')[0]
                    if chr_name not in chr_names:
                        chr_names.append(chr_name)
                        print(chr_name)

    elif option == 'types':
        types_by_chr = {}
        with open(_input, 'r') as gff_file:
            for line in gff_file:
                if not line.startswith('#'):
                    chr_name = line.split('\t')[0]
                    _type = line.split('\t')[2]
                    if chr_name not in types_by_chr:
                        types_by_chr[chr_name] = []
                        print(chr_name)
                    if _type != 'region' and _type not in types_by_chr[chr_name]:
                        types_by_chr[chr_name].append(_type)
                        print(' - ' + _type)
