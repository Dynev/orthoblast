import os
import copy
import pygame as pg
from pygame.locals import *

from . import misc

def drawResult(spec, prot, ev, cov, result, outname):
    def drawRect(x1, x2, y, color):
        pg.draw.rect(screen, misc.COLOR['black'], Rect(x1, y, x2-x1, 14))
        pg.draw.rect(screen, misc.COLOR[color], Rect(x1+2, y+2, x2-x1-4, 10))
    def prepHit(hit):
        chit = copy.deepcopy(hit)
        for h in chit.hsps:
            h.qstart /= chit.query.qlen / 500
            h.qend /= chit.query.qlen / 500
            h.qend = round(h.qend)
        if chit.domains != []:
            for d in chit.domains:
                d.start /= chit.query.qlen / 500
                d.end /= chit.query.qlen / 500
        for i in range(1, len(chit.hsps)):
            if hit.hsps[i].qstart - hit.hsps[i-1].qend == 1:
                chit.hsps[i].qstart = chit.hsps[i-1].qend - 2
        return chit
    def prepQuery(query):
        q = copy.deepcopy(query)
        for d in q.domains:
            d.start /= q.qlen / 500
            d.end /= q.qlen / 500
        return q
    def drawTitle():
        title = font.render(f'HUMAN {prot.upper()} VS {spec.upper()} GENOME (E={ev if ev else 10}, COV={cov if cov else "CUSTOM"})', True, misc.COLOR['black'])
        screen.blit(title, (10, 7))
    def drawHeatMap():
        for i in range(1, 10):
            pg.draw.rect(screen, misc.COLOR[f'hm{i}'], Rect(346+i*24, 6, 24, 10))
            score = font.render(str(25*i), True, misc.COLOR['black'])
            screen.blit(score, (358+i*24-score.get_width()/2, 18))
    def drawDomain(x1, x2, y, name):
        drawRect(x1, x2, y, 'e')
        dname = font.render(name, True, misc.COLOR['black'])
        coord = (x1 + (x2-x1)/2) - dname.get_width()/2
        screen.blit(dname, (coord, y+2))
    def drawQuery(query):
        query = prepQuery(query)
        drawRect(50, 550, 32, 'c')
        qname = font.render('QUERY', True, misc.COLOR['black'])
        qsize = font.render(str(query.qlen)+'bp', True, misc.COLOR['black'])
        screen.blit(qname, (25-qname.get_width()/2, 34))
        screen.blit(qsize, (575-qsize.get_width()/2, 34))
        for d in query.domains:
            drawDomain(50+d.start, 50+d.end, 32, d.name)
    def drawHit(y, hit, flag):
        hit = prepHit(hit)
        for hsp in hit.hsps:
            if hsp.score > 200:
                color = 'hm9'
            elif hsp.score > 175:
                color = 'hm8'
            elif hsp.score > 150:
                color = 'hm7'
            elif hsp.score > 125:
                color = 'hm6'
            elif hsp.score > 100:
                color = 'hm5'
            elif hsp.score > 75:
                color = 'hm4'
            elif hsp.score > 50:
                color = 'hm3'
            elif hsp.score > 25:
                color = 'hm2'
            else:
                color = 'hm1'
            drawRect(50+hsp.qstart, 50+hsp.qend, y, color)
        nhit = font.render('HIT', True, misc.COLOR['black'])
        nseq = font.render('SEQ', True, misc.COLOR['black'])
        screen.blit(nhit, (25-nhit.get_width()/2, y+2))
        screen.blit(nseq, (575-nseq.get_width()/2, y+2))
        if flag:
            y += 20
            for d in hit.domains:
                drawDomain(50+d.start, 50+d.end, y, d.name)
            ndom = font.render('DOM', True, misc.COLOR['black'])
            screen.blit(ndom, (575-ndom.get_width()/2, y+2))
    pg.init()
    length = 40*len([hit for hit in result.hits if hit.seq != '']) + 20*len([hit for hit in result.hits if hit.seq == ''])
    screen = pg.display.set_mode((600, 60+length))
    screen.fill((255,255,255))
    font = pg.font.SysFont('Arial Bold', 15)
    drawTitle()
    drawHeatMap()
    drawQuery(result.query)
    y = 60
    for i in range(len(result.hits)):
        if result.hits[i].seq != '':
            drawHit(y, result.hits[i], True)
            y += 40
        else:
            drawHit(y, result.hits[i], False)
            y += 20

    if not os.path.isdir(os.path.join(os.getcwd(), 'results', f'{outname if outname else ""}', 'graphs')):
        os.mkdir(os.path.join(os.getcwd(), 'results', f'{outname if outname else ""}', 'graphs'))
    pg.image.save(screen, os.path.join(os.getcwd(), 'results', f'{outname if outname else ""}', 'graphs', f'{spec}_{prot}.png'))
    pg.quit()
