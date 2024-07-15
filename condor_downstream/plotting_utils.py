import re
from Bio.SeqUtils import seq1

def rgb_string_to_hex(rgb):
    """
    e.g. input: 'rgb(141,211,199)'
    """
    rgb = tuple(map(int, rgb[4:-1].split(',')))
    return '#%02x%02x%02x' % rgb

def __hgvsp_3to1(hgvsp_string):
    """
    1. identify all the 3-letter amino acid codes (3 letters after every capital letter)
    2. replace those 3-letter amino acid codes with 1-letter amino acid codes using Bio.SeqUtils.se1(). Keep other characters intact
    """
    pattern = re.compile(r"[A-Z][a-z][a-z]")
    return re.sub(pattern, lambda x: seq1(x.group()), hgvsp_string)