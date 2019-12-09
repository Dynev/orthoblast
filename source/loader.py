import os
import re
import gzip
import shutil
import requests

from . import blaster
from . import misc

def verifyDirs(outname):
    if not os.path.isdir(os.path.join(os.getcwd(), 'proteins')):
        os.mkdir(os.path.join(os.getcwd(), 'proteins'))
    if not os.path.isdir(os.path.join(os.getcwd(), 'genomes')):
        os.mkdir(os.path.join(os.getcwd(), 'genomes'))
    if not os.path.isdir(os.path.join(os.getcwd(), 'results')):
        os.mkdir(os.path.join(os.getcwd(), 'results'))
    if outname:
        if not os.path.isdir(os.path.join(os.getcwd(), 'results', outname)):
            os.mkdir(os.path.join(os.getcwd(), 'results', outname))
    return True

def checkProtein(prot):
    adr = os.path.join(os.getcwd(), 'proteins', prot + '.fasta')
    if os.path.exists(adr):
        return True
    else:
        return False

def checkGenome(spec):
    adr = os.path.join(os.getcwd(), 'genomes', spec)
    if os.path.isfile(os.path.join(adr, f'{spec}.nsq')):
        return True
    else:
        return False

def loadProtein(prot):
    full = prot + '_human.fasta'
    req = requests.get(misc.ADDRESS['uniprot'] + full)
    ret = req.text
    if '|' in ret[:40]:
        ret = '>' + ret[4:]
        prot_path = os.path.join(os.getcwd(), 'proteins')
        with open(os.path.join(prot_path, f'{prot}.fasta'), 'w+') as f:
            f.write(ret)
        return True
    else:
        return False

def installGenome(spec):
    if spec not in misc.SPEC.keys():
        return False
    
    adr = misc.ADDRESS['ncbi'] + misc.SPEC[spec] + '/latest_assembly_versions/'
    for i in (misc.FILTER['chrom'], misc.FILTER['genome_file']):
        req = requests.get(adr)
        ret = req.text
        res = i.search(ret).group(0)
        adr = adr + res
    
    with requests.get(adr, stream = True) as req:
        with open(os.path.join(os.getcwd(), 'genomes', f'{spec}.fna.gz'), 'wb+') as gz:
            shutil.copyfileobj(req.raw, gz)
    with gzip.open(os.path.join(os.getcwd(), 'genomes', f'{spec}.fna.gz')) as gz:
        with open(os.path.join(os.getcwd(), 'genomes', f'{spec}.fna'), 'wb+') as fna:
            shutil.copyfileobj(gz, fna)
    os.remove(os.path.join(os.getcwd(), 'genomes', f'{spec}.fna.gz'))
    blaster.runMakeDB(os.path.join(os.getcwd(), 'genomes', f'{spec}.fna'), f'{spec}')
    os.remove(os.path.join(os.getcwd(), 'genomes', f'{spec}.fna'))
    return True

def loadSeq(prot):
    adr = os.path.join(os.getcwd(), 'proteins', f'{prot}.fasta')
    with open(adr, 'r') as f:
        f = f.read()
        seq_start = misc.FILTER['seq'].search(f).start()
        seq = f[seq_start:-1]
        seq = re.sub('\n', '', seq)
    return seq