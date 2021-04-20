#!/usr/bin/env python3
import sys, os.path, os, re, glob, encodings
from operator import itemgetter
from optparse import OptionParser
from collections import defaultdict

from psmparse import getParser, getParsers, AddMZMLFields
from psmformat import SequenceDatabase, PSMCache

psmparsers = getParsers()

from version import VERSION
parser = OptionParser(version=VERSION)
parser.add_option("--format",type="choice",choices=psmparsers,default=psmparsers[0],dest="psmfmt",help="PSM input file format")
parser.add_option("--seqdb",type="string",default="RefSeq:Human,UniProt:Human",
                    dest="seqdb",help="Comma separated list of sequence databases to remap peptides against. Default: RefSeq:Human,UniProt:Human.")
parser.add_option("--specdir",type="string",default=None,
                    dest="specdir",help="Directory containing the spectra. Default: Same as .psm file(s).")
parser.add_option("--psmfilter",type="string",default=None,
                  dest="psmfilter",help="Criteria for accepting PSMs (Python expression applied original data). Default: None.")
parser.add_option("--seqdir",type="string",default="sequence",
                  dest="seqdir",help="Directory containing the sequence files. Default: \"sequence\".")

opts,args = parser.parse_args()
if opts.seqdb:
    opts.seqdb = [_f for _f in map(str.strip,opts.seqdb.split(',')) if _f]
else:
    opts.seqdb = []

PSMParser = getParser(opts.psmfmt)
psms = PSMParser(args,filter=opts.psmfilter)

if len(set(psms.seqdbs())|set(opts.seqdb)) > 0:
    opts.seqdb.extend(psms.seqdbs())
seqdbs = SequenceDatabase.seqdbs(opts.seqdir,opts.seqdb)
if "Decoy" in seqdbs:
    seqdbs["Decoy"].setDecoyPrefix(psms.decoyPrefix())
# print seqdbs

scores = psms.scores
params = psms.params
specparams = AddMZMLFields.specparams

md = psms.metadata()
md['Software'].append("psmextract v%s"%VERSION)

psmcache = PSMCache(AddMZMLFields(psms,opts.specdir),scores,params,specparams)

for sdb in list(seqdbs.values()):
    sdb.peptidemap(psmcache.peptides())

# Perhaps we will want a more generic one for non-CPTAC uses?
from psmformat import PSMFormater

fmt = PSMFormater(scores,params,specparams)
fmt.write_metadata(md)

from StringIndex import StringIndex
protindex = StringIndex()
peptindex = StringIndex()
prot2pept = defaultdict(set)

sdbseenids = set()
for psm in psmcache.psms():

    peptid = peptindex.add(psm['Peptide'])

    origProteins = psm['Protein']
    psm['_Protein'] = list(origProteins)
    newProteins = []

    # print origProteins

    sdbneeded = set()
    for sdb in list(seqdbs.values()):
        # print sdb.id(),sdb.remapper()
        if not sdb.remapper():
            continue
        for pr in sdb.proteins(psm['Peptide']):
            sdbneeded.add(sdb.id())
            origProteins = psms.removeProtein(sdb,pr,origProteins)
            newProteins.append(pr)

    # print origProteins
    # print newProteins

    if len(origProteins) > 0:
        sdbneeded.update(psms.seqdbs(origProteins))

    for sdbid in sdbneeded:
        if sdbid not in sdbseenids:
            sdbseenids.add(sdbid)
            fmt.write_seqdb(seqdbs[sdbid])

    psm['Protein'] = newProteins
    psm['Protein'].extend(origProteins)

    for pr in psm['Protein']:
        protid = protindex.add(":".join(pr[:2]))
        prot2pept[protid].add(peptid)

    fmt.write_psm(psm)

grpmd = dict()
grpmd['AnalysisSoftware'] = "psmextract"
fmt.write_prgrp_metadata(grpmd)

from parsimony import Dominator, Components

allprot = []
for prid in protindex:
    sdbid,pracc = protindex.string(prid).rsplit(':',1)
    allprot.append([prid,sdbid,pracc])

def cmp(a,b):
    if a < b:
        return -1
    elif a > b:
        return 1
    return 0

def prcmp(a,b):
    # compare sequence databases and choose the one with smallest
    # priority
    c = cmp(seqdbs[a[1]].priority(),seqdbs[b[1]].priority())
    if c != 0:
        return c
    # Sequence databases have equal priority, they must be
    # the same sequence database.
    assert(a[1] == b[1])
    if seqdbs[a[1]]._acc:
        dla = seqdbs[a[1]].protein_defline(a[2])
        dlb = seqdbs[b[1]].protein_defline(b[2])
        return seqdbs[a[1]].prefer((a[0],dla),(b[0],dlb))
    return cmp(a[0],b[0])

from functools import cmp_to_key

allprot.sort(key=cmp_to_key(prcmp))

pr2key = dict()
for i,prid in enumerate(map(itemgetter(0),allprot)):
    pr2key[prid] = i

dom = Dominator(peptides=set(peptindex),edges=prot2pept,proteins=set(protindex))
dom.dominate(prsortkey=pr2key.get)

selected=set(dom.equivalentto)

comp = Components(peptides=set(peptindex),edges=prot2pept,proteins=set(protindex))
for i,(prids,pepids) in enumerate(sorted(comp,key=lambda t: -max([len(prot2pept[prid]) for prid in t[0]]))):
    prdata = []
    prseen = set()
    first = True
    for prid in sorted(prids&selected,key=lambda prid: -len(prot2pept[prid])):
        if prid in prseen:
            continue
        prstr = protindex.string(prid)
        sdbid,pracc = prstr.rsplit(':',1)
        coverage = seqdbs[sdbid].protein_coverage(pracc)
        if coverage != None:
            coverage *= 100.0
        prdata.append(dict(id=prstr,anchor=first,
                           peptides=len(prot2pept[prid]),
                           coverage=coverage))
        first=False
        prseen.add(prid)
        for prid1 in sorted(dom.equivalentto[prid],key=pr2key.get):
            if prid1 == prid or prid1 in prseen:
                continue
            prstr1 = protindex.string(prid1)
            sdbid1,pracc1 = prstr1.rsplit(':',1)
            coverage = seqdbs[sdbid1].protein_coverage(pracc1)
            if coverage != None:
                coverage *= 100.0
            prdata.append(dict(id=prstr1,anchor=False,sameset=prstr,
                               peptides=len(prot2pept[prid1]),
                               coverage=coverage))
            prseen.add(prid1)
        for prid1 in set(dom.containedby[prid] - prseen):
            prstr1 = protindex.string(prid1)
            sdbid1,pracc1 = prstr1.rsplit(':',1)
            coverage = seqdbs[sdbid1].protein_coverage(pracc1)
            if coverage != None:
                coverage *= 100.0
            prdata.append(dict(id=prstr1,anchor=False,subset=prstr,
                               peptides=len(prot2pept[prid1]),
                               coverage=coverage))
            prseen.add(prid1)

    fmt.write_prgrp(i+1,len(pepids),prdata)
