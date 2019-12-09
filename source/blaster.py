import os
import re
from bs4 import BeautifulSoup as bs
from Bio.Blast.Applications import NcbitblastnCommandline as tblastn
from Bio.Blast.Applications import NcbimakeblastdbCommandline as makeblastdb

from . import misc

def runBlast(spec, prot, ev, outname):
    output = os.path.join(os.getcwd(), 'results', outname if outname else '', 'output', f'{spec}_{prot}.xml')
    if os.path.isfile(output):
        with open(output) as res:
            res = str(bs(res, 'lxml'))
            tev = misc.FILTER['evalue'].search(res).group(0)
        if ev == float(tev):
            return
    
    cline = tblastn(query = os.path.join(os.getcwd(), 'proteins', f'{prot}.fasta'), 
                    db = os.path.join(os.getcwd(), 'genomes', spec, spec), 
                    out = output,
                    evalue = ev,
                    num_threads = 4,
                    outfmt = 5)
    cline()

def runMakeDB(inname, outname):
    cline = makeblastdb(input_file = inname, dbtype = 'nucl', out = outname)
    cline()