import os
import re
import sys
import copy
import gzip
import shutil
import operator
import subprocess

os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import numpy
import pandas
import pygame
import requests
from bs4 import BeautifulSoup as bs
from Bio.Blast.Applications import NcbitblastnCommandline as tblastn

from source import misc
from source import drawer
from source import loader
from source import parser
from source import blaster
from source import interpreter

def main():
    mode = parser.getMode(sys.argv[1:])
    if isinstance(mode, str): 
        print(mode)
        return

    check = loader.verifyDirs(mode.outname)
    if isinstance(check, str):
        print(check)
        return

    print(f'Running set: {len(mode.proteins)} proteins across {len(mode.species)} species')
    
    task = len(mode.proteins)
    toremove = []
    result = ['.fasta already exists', '.fasta downloaded', 
              '.fasta download failed']
    status = 0
    pfin = False
    for i in range(len(mode.proteins)):    
        print(f'Acquiring proteins: {mode.proteins[i]} ({i+1}/{task})', end = '\r')
        if not loader.checkProtein(mode.proteins[i]):
            if loader.loadProtein(mode.proteins[i]):
                status = 1
            else:
                status = 2
                toremove.append(i)
        else: status = 0
    if toremove:
        if len(toremove) == len(mode.proteins):
            print('None of specified proteins were acquired. Exiting     ')
            return
        else:
            msg = f'{mode.proteins[toremove[0]]}'
            if len(toremove) > 1:
                for i in toremove[1:]:
                    msg += f', {mode.proteins[i]}'
            msg += (f'{" was" if len(toremove) == 1 else " were"} not downloaded, '
                    f'continuing with {len(mode.proteins) - len(toremove)}')
            print(msg)
            for i in sorted(toremove, reverse = True):
                if mode.pdf:
                    mode.pdf = mode.pdf[mode.pdf['Proteins'] != mode.proteins[i]]
                mode.proteins.pop(i)
            pfin = True
    else:
        print('All proteins were successfully acquired          ')
        
    task = len(mode.species)
    toremove = []
    result = [' database already exists', ' database installed', 
              ' database installation failed']
    sfin = False
    for i in range(len(mode.species)):
        print(f'Acquiring genomes: {mode.species[i]} ({i+1}/{task})                     ', end = '\r')
        if not loader.checkGenome(mode.species[i]):
            if loader.installGenome(mode.species[i]):
                status = 1
            else:
                status = 2
                toremove.append(i)
        else: status = 0
    if toremove:
        if len(toremove) == len(mode.species):
            print('None of specified genomes were acquired. Exiting             ')
            return
        else:
            msg = f'{mode.species[toremove[0]]}'
            if len(toremove) > 1:
                for i in toremove[1:]:
                    msg += f', {mode.species[i]}'
            msg += (f'{" was" if len(toremove) == 1 else " were"} not downloaded, '
                    f'continuing with {len(mode.species) - len(toremove)}          ')
            print(msg)
            for i in sorted(toremove, reverse = True):
                mode.species.pop(i)
            sfin = True
    else:
        print('All genomes were successfully acquired               ')

    if not pfin and not sfin:
        print('All proteins and genomes were successfully acquired. Continuing')
    else:
        answer = ''
        while True:
            answer = input(f'Some {"proteins" if pfin else ""} '
                           f'{"and" if pfin and sfin else ""} {"genomes" if sfin else ""} '
                           f'were not acquired. Continue with reduced set? [y/n] ')
            if answer not in ['y', 'n']:
                continue
            else:
                if answer == 'n':
                    return
                if answer == 'y':
                    break

    task = len(mode.proteins) * len(mode.species)
    count = 1
    for s in mode.species:
        for p in mode.proteins:
            print(f'Running TBLAST: {p} vs {s} ({count}/{task})       ', end = '\r')
            blaster.runBlast(s, p, mode.evalue, mode.outname)
            count += 1
    print('TBLAST run completed                    ')

    rdf = pandas.DataFrame()
    count = 1
    for s in mode.species:
        results = []
        for p in mode.proteins:
            print(f'Analyzing results: {p} vs {s} ({count}/{task})           ', end = '\r')
            result = interpreter.interpretBlast(s, p, mode.coverage, mode.outname)
            results.append(result)
            count += 1
        rdf[s] = pandas.Series(results)
    print(f'Analysis of results completed')

    count = 1
    if mode.draw == True:
        for s in mode.species:
            for p in mode.proteins:
                print(f'Creating graphs ({count}/{task})', end ='\r')
                drawer.drawResult(s, p, rdf[s], mode.outname)
        print(f'Creation of graphs completed')

    for s in mode.species:
        rdf[s] = misc.groupattr(rdf[s], 'nhomo')
    if mode.pdf:
        rdf = pandas.concat([mode.pdf, rdf], axis = 1)
    else:
        rdf['Proteins'] = mode.proteins
    rdf.to_excel(os.path.join(os.getcwd(), 'results', mode.outname if mode.outname else '', 'result.xlsx'), index = False)
    print(f'Run finished. Results are in results/{mode.outname if mode.outname else ""}')

if __name__ == "__main__":
    main()