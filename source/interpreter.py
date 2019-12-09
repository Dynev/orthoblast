import os
import re
import operator as op
import requests as rq
from bs4 import BeautifulSoup as bs

from . import loader
from . import misc

class Domain:
    def __init__(self, name, start, end):
        self.name = name
        self.start = int(start)
        self.end = int(end)

def findDomains(seq):
    req = rq.post(misc.ADDRESS['hmmer'], params = {'hmmdb':'pfam','seq':seq},
                  allow_redirects = False)
    adr = req.headers['location']
    req = rq.get(adr, params = {'output':'xml'})
    ret = str(bs(req.text, 'lxml'))
    names = misc.FILTER['domain'].findall(ret)
    starts = misc.FILTER['dom_start'].findall(ret)
    ends = misc.FILTER['dom_end'].findall(ret)
    isincluded = re.compile(r'(?<=is_included=")(\d)').findall(ret)
    isoutcompeted = re.compile(r'(?<=outcompeted=")(\d)').findall(ret)
    domains = zip(names, starts, ends)
    domains = [dom for i,dom in enumerate(domains) if isincluded[i] == '1']
    domains = [dom for i,dom in enumerate(domains) if isoutcompeted[i] == '0']
    domains = [Domain(*dom) for dom in domains]
    return domains

def loadSeq(prot):
    with open(os.path.join(os.getcwd(), 'proteins', f'{prot}.fasta')) as fasta_file:
        fasta_file = fasta_file.read()
        seq_filter = re.compile('[A-Z]{10}')
        seq_begin = seq_filter.search(fasta_file).start()
        seq = fasta_file[seq_begin:-1]
    return re.sub('\n', '', seq)

class Hsp:
    def __init__(self, qstart, qend, seq, score, frame):
        self.qstart = qstart
        self.qend = qend
        self.seq = re.sub('-', '', seq)
        self.score = score
        self.frame = frame
        self.truncate()
    def __repr__(self):
        return f'q{self.qstart}-{self.qend} ({self.alen()}); score={self.score}; frame={self.frame}\n{self.seq} ({len(self.seq)})'
    def alen(self):
        return self.qend - self.qstart + 1
    def truncate(self):
        l = self.alen()
        if l > len(self.seq):
            self.seq += 'X'*(l-len(self.seq))
        elif l < len(self.seq):
            self.seq = self.seq[:l]

class Segview:
    def __init__(self, hsp1, qlen):
        self.qlen = qlen
        self.view = [[-1, -1], [hsp1.qstart, hsp1.qend], [qlen+2, qlen+2]]
    def reclaim(self, new):
        self.view = sorted(self.view, key = lambda x: x[0])
        for i in range(new+2):
            for j in range(len(self.view) - 1):
                if self.view[j][1] == self.view[j+1][0] - 1:
                    self.view[j+1][0] = self.view[j][0]
                    self.view.pop(j)
                    break
    def insert(self, hsp):
        for i in range(len(self.view) - 1):
            if hsp.qstart >= self.view[i][0]:
                spos = i
                if hsp.qstart <= self.view[i][1]:
                    bsin = True
                else:
                    bsin = False
            if hsp.qend > self.view[i][1]:
                epos = i
                if hsp.qend >= self.view[i+1][0]:
                    bein = True
                else:
                    bein = False
        starts = []
        ends = []
        for i in range(spos, epos+1):
            if i == spos:
                if bsin == True:
                    starts.append(self.view[i][1] + 1)
                else:
                    starts.append(hsp.qstart)
            else:
                starts.append(self.view[i][1] + 1)
            if i == epos:
                if bein == True:
                    ends.append(self.view[i+1][0] - 1)
                else:
                    ends.append(hsp.qend)
            else:
                ends.append(self.view[i+1][0] - 1)
        result = [list(x) for x in zip(starts, ends)]
        for r in result:
            self.view.append(r)
        self.reclaim(len(result))
        result = [list(x) for x in zip(starts, ends)]
        return result

class Hit:
    def __init__(self, tcov, query, hsps):
        self.tcov = tcov
        self.query = query
        self.hsps = hsps
        self.qcov = 0
        self.seq = ''
        self.domains = []
        self.homo = False
        self.converge()
        self.relap()
        self.detHomo()
    def __repr__(self):
        rep = ''
        for i, hsp in enumerate(self.hsps):
            rep += f'Hsp {i+1}:\n{hsp}\n'
        if self.qcov != 0:
            rep += f'Query Coverage: {self.qcov}%\n'
        if self.seq != '':
            rep += f'Full Seq:\n{self.seq} ({len(self.seq)})\n'
        if self.domains != []:
            rep += f'Domains:\n'
            for d in self.domains:
                rep += f'{d.name} ({d.start}-{d.end})\n'
        return rep
    def converge(self):
        fp = 0
        fn = 0
        for hsp in self.hsps:
            if hsp.frame > 0:
                fp += hsp.alen()*hsp.score
            else:
                fn += hsp.alen()*hsp.score
        if fp > fn:
            self.hsps = [hsp for hsp in self.hsps if hsp.frame > 0]
        else:
            self.hsps = [hsp for hsp in self.hsps if hsp.frame < 0]
    def relap(self):
        self.hsps = sorted(self.hsps, key = op.attrgetter('score'), reverse = True)
        s = Segview(self.hsps[0], self.query.qlen)
        todelete = []
        if len(self.hsps) != 1:
            for i in range(1, len(self.hsps)):
                ncoords = s.insert(self.hsps[i])
                if ncoords == []:
                    todelete.append(i)
                elif len(ncoords) == 1:
                    sshift = ncoords[0][0] - self.hsps[i].qstart
                    eshift = self.hsps[i].qend - ncoords[0][1]
                    self.hsps[i].qstart = ncoords[0][0]
                    self.hsps[i].qend = ncoords[0][1]
                    nslice = slice(sshift, (-eshift if eshift != 0 else None))
                    self.hsps[i].seq = self.hsps[i].seq[nslice]
                else:
                    for n in ncoords:
                        sshift = n[0] - self.hsps[i].qstart
                        eshift = self.hsps[i].qend - n[1]
                        nslice = slice(sshift, (-eshift if eshift != 0 else None))
                        self.hsps.append(Hsp(n[0], n[1], self.hsps[i].seq[nslice], self.hsps[i].score, self.hsps[i].frame))
                    todelete.append(i)
        if todelete != []:
            for j in sorted(todelete, reverse = True):
                self.hsps.pop(j)
    def getCov(self):
        for hsp in self.hsps:
            self.qcov += hsp.alen()
        self.qcov /= self.query.qlen
        self.qcov = round(self.qcov, 2)
    def getSeq(self):
        self.hsps = sorted(self.hsps, key = op.attrgetter('qstart'))
        pos = 1
        for hsp in self.hsps:
            if hsp.qstart == pos:
                self.seq += hsp.seq
            else:
                self.seq += 'X'*(hsp.qstart - pos)
                self.seq += hsp.seq
            pos = hsp.qend + 1
    def getDomains(self):
        self.getSeq()
        self.domains = findDomains(self.seq)
    def detHomo(self):
        self.getCov()
        if self.tcov:
            if self.qcov >= self.tcov:
                self.getDomains()
                hdn = misc.groupattr(self.domains, 'name')
                qdn = misc.groupattr(self.query.domains, 'name')
                if len(self.domains) == len(self.query.domains):
                    if set(hdn) == set(qdn):
                        self.homo = True
                    else:
                        self.homo = False
                else:
                    self.homo = False
            else:
                self.homo = False
        else:
            if True: #self.qcov >= 0.5:
                if (self.qcov >= self.query.dlen / self.query.qlen - 0.1):
                    self.getDomains()
                    hdn = misc.groupattr(self.domains, 'name')
                    qdn = misc.groupattr(self.query.domains, 'name')
                    if len(self.domains) == len(self.query.domains):
                        if set(hdn) == set(qdn):
                            self.homo = True
                        else:
                            self.homo = False
                    else:
                        self.homo = False
                else:
                    self.homo = False

class Query:
    def __init__(self, qlen, name):
        self.qlen = qlen
        self.name = name
        self.seq = loadSeq(name)
        self.domains = findDomains(self.seq)
        self.dlen = self.getDlen()
    def getDlen(self):
        dlen = 0
        for d in self.domains:
            dlen += d.end - d.start
        return dlen

class Result: 
    def __init__(self, query, hits):
        self.query = query 
        self.hits = hits
        self.nhomo = 0
        self.countHomo()
    def countHomo(self):
        self.nhomo = misc.groupattr(self.hits, 'homo').count(True)

def interpretBlast(spec, prot, coverage, outname):
    with open(os.path.join(os.getcwd(), 'results', outname if outname else '', 'output', f'{spec}_{prot}.xml')) as res:
        res = str(bs(res, 'lxml'))
        qlen = int(misc.FILTER['query_len'].search(res).group(0))
        query = Query(qlen, prot)
        hitsl = re.finditer(misc.FILTER['hit'], res)
        hitsl = [h for h in hitsl]
        hits = []
        for i in range(len(hitsl)):
            hitsp = res[hitsl[i].start():(hitsl[i+1].start() if (i < len(hitsl) - 1) else len(res))]
            hspsl = re.finditer(misc.FILTER['hsp'], hitsp)
            hspsl = [h for h in hspsl]
            hsps = []
            for j in range(len(hspsl)):
                hspsp = hitsp[hspsl[j].start():(hspsl[j+1].start() if (j < len(hspsl) - 1) else len(hitsp))]
                qf = int(misc.FILTER['hsp_query-from'].search(hspsp).group(0))
                qt = int(misc.FILTER['hsp_query-to'].search(hspsp).group(0))
                sq = str(misc.FILTER['hsp_seq'].search(hspsp).group(0))
                sc = float(misc.FILTER['hsp_score'].search(hspsp).group(0))
                fr = int(misc.FILTER['hsp_frame'].search(hspsp).group(0))
                hsps.append(Hsp(qf, qt, sq, sc, fr))
            hits.append(Hit(coverage, query, hsps))
    return Result(query, hits)