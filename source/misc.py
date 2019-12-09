import re

SPEC = {'het_gla': 'Heterocephalus_glaber',
       'cav_por': 'Cavia_porcellus',
       'hom_sap': 'Homo_sapiens', }

COLOR = {'b': (221, 132, 82),
        'c': (85, 168, 104),
        'd': (196, 78, 82),
        'e': (129, 114, 179),
        'f': (147, 120, 96),
        'g': (218, 139, 195),
        'h': (140, 140, 140),
        'i': (204, 185, 116),
        'j': (100, 181, 205),
        'p': (228, 228, 239),
        'r': (157, 157, 206),
        'black': (0, 0, 0),
        'white': (255, 255, 255),
        'hm1': (69, 117, 180),
        'hm2': (124, 170, 208),
        'hm3': (195, 224, 237),
        'hm4': (233, 246, 230),
        'hm5': (250, 253, 199),
        'hm6': (254, 235, 161),
        'hm7': (252, 178, 113),
        'hm8': (241, 115, 75),
        'hm9': (215, 48, 39),} 

FILTER = {'evalue': re.compile(r'(?<=<parameters_expect>)(.*?)(?=</parameters_expect>)'),
         'chrom': re.compile(r'GCF(.*?)/'),
         'genome_file': re.compile(r'GCF(.*?)\d_genomic\.fna\.gz'),
         'pfam-result': re.compile(r'(?<=<result_url>)(.*?)(?=</result_url>)'),
         'domain': re.compile(r'(?<=alihmmname=")(.*?)(?=")'),
         'dom_start': re.compile(r'(?<=ienv=")(.*?)(?=")'),
         'dom_end': re.compile(r'(?<=jenv=")(.*?)(?=")'),
         'hit': r'<hit_num>',
         'hsp': r'<hsp_num>',
         'query_len': re.compile(r'(?<=<blastoutput_query-len>)(.*?)(?=</blastoutput_query-len>)'),
         'hsp_evalue': re.compile(r'(?<=<hsp_evalue>)(.*?)(?=</hsp_evalue>)'),
         'hsp_score': re.compile(r'(?<=hsp_bit-score>)(.*?)(?=</hsp_bit-score)'),
         'hsp_seq': re.compile(r'(?<=hsp_hseq>)(.*?)(?=</hsp_hseq)'),
         'hsp_frame': re.compile(r'(?<=hsp_hit-frame>)(.*?)(?=</hsp_hit-frame)'),
         'hsp_query-from': re.compile(r'(?<=hsp_query-from>)(.*?)(?=</hsp_query-from)'),
         'hsp_query-to': re.compile(r'(?<=hsp_query-to>)(.*?)(?=</hsp_query-to)'),
         'hsp_hit-from': re.compile(r'(?<=hsp_hit-from>)(.*?)(?=</hsp_hit-from)'),
         'hsp_hit-to': re.compile(r'(?<=hsp_hit-to>)(.*?)(?=</hsp_hit-to)'), }

ADDRESS = {'uniprot': 'https://www.uniprot.org/uniprot/',
          'hmmer': 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan',
          'ncbi': 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/',}

def groupattr(insts, attr):
    ret = []
    for i in insts:
        ret.append(eval(f'i.{attr}'))
    return ret