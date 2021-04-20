#!/usr/bin/env python27
import sys, os.path, os, re, glob
from operator import itemgetter
from optparse import OptionParser
from version import VERSION
from collections import defaultdict
import mzml

VERSION='1.2.1'
parser = OptionParser(version=VERSION)
parser.add_option("--taxa",type="string",default="Human",
                    dest="taxa",help="Comma separated taxa from Human, Mouse. Default: Human.")
parser.add_option("--spectradir",type="string",default=None,
                    dest="tsvdir",help="Directory containing the spectra. Default: Same as .psm file.")
parser.add_option("--source",type="string",default="UniProt,RefSeq",
                    dest="source",help="Comma separated sources from UniProt, RefSeq. Default: UniProt,RefSeq.")
parser.add_option("--prgrp",action="store_true",default=False,
		  dest="prgrp",help="Compute protein groups")
opts,args = parser.parse_args()
if opts.taxa:
    opts.taxa = filter(None,map(str.strip,opts.taxa.split(',')))
else:
   opts.taxa = []
if opts.source:
    opts.source = filter(None,map(str.strip,opts.source.split(',')))
else:
   opts.source = []

from csv import DictReader

# FileName ScanNum QueryPrecursorMz OriginalPrecursorMz
# PrecursorError(ppm) QueryCharge OriginalCharge PrecursorScanNum
# PrecursorArea PrecursorRelAb RTAtPrecursorHalfElution PeptideSequence
# AmbiguousMatch Protein DeNovoScore MSGFScore Evalue Qvalue PepQvalue
# PrecursorPurity FractionDecomposition HCDEnergy iTRAQ114 iTRAQ115
# iTRAQ116 iTRAQ117 iTRAQFlags iTRAQTotalAb iTRAQFractionOfTotalAb
# PhosphoRSPeptide nPhospho FullyLocalized

ignore = """
ScanNum FileName Protein QueryCharge PeptideSequence OriginalCharge
OriginalPrecursorMz PrecursorScanNum PhospoRSPeptide
""".split()
scores = """
    DeNovoScore MSGFScore SpecEvalue Evalue Qvalue PepQvalue 
""".split()
params = """
    QueryPrecursorMz 
    PrecursorError(ppm)
    PrecursorArea PrecursorRelAb RTAtPrecursorHalfElution
    AmbiguousMatch PrecursorPurity FractionDecomposition HCDEnergy
    iTRAQ114 iTRAQ115 iTRAQ116 iTRAQ117
    iTRAQFlags iTRAQTotalAb iTRAQFractionOfTotalAb
    PhosphoRSPeptide nPhospho FullyLocalized
""".split()
specparams=['intensity of precursor ion','retention time']
# specparams=[]
headers = [ 'SpectrumFile', 'Location', 'SpectrumID', 'SpectrumIDFormat', 'FileFormat', 'Scan', 
	    'ExperimentalMassToCharge', 'Rank', 'Peptide', 
	    'retention time', 'intensity of precursor ion',
	    'ChargeState', 'Modification', 'Protein' ] + scores + params + specparams
expected_keys = set(scores + params + ignore)

ptms = """
[	+144.102	UNIMOD:214	144.102063	name
K	+144.102	UNIMOD:214	144.102063	name
M	+15.995		Oxidation       15.994915	name
C	+57.021		Carbamidomethyl 57.021464	name
Q	+0.984		Deamidated	0.984016	name
N	+0.984		Deamidated	0.984016	name
[Q	-17.027		Gln->pyro-Glu	-17.026549 	name
[Q	-17.026		Gln->pyro-Glu	-17.026549 	name
[E	-18.011		Glu->pyro-Glu 	-18.010565	name
S	+79.966		Phospho		79.966331 	name
T	+79.966		Phospho		79.966331 	name
Y	+79.966		Phospho		79.966331 	name
[       +42.011         Acetyl          42.010565       name
"""
# [Q	+144.102-17.027 UNIMOD:214,Gln->pyro-Glu 144.102063,-17.026549 name,name  
# [E	+144.102-18.011 UNIMOD:214,Glu->pyro-Glu 144.102063,-18.010565 name,name  
cvmods = {}
for l in ptms.splitlines():
    if l.strip() == "":
	continue
    aa,delstr,name,delta,spec = l.split()
    delta = float(delta)
    cvmods[(aa,delstr)] = (name,delta,spec)

root = os.path.split(sys.argv[0])[0]
seqdb = defaultdict(list)
for t in opts.taxa:
  for s in opts.source:
    matches = glob.glob('%s/sequence/%s:%s:*.fasta'%(root,s,t))
    if len(matches) == 1:
        seqdb[s].extend(matches)
    if len(matches) > 1:
        raise RuntimeError("Sequence database problem, too many matching databases for %s"%('%s/sequence/%s:%s:*.fasta'%(root,s,t),))
	
md = dict()
md['SpectrumIDFormat'] = 'Thermo nativeID format'
# md['FileFormat'] = 'Thermo RAW file'
md['FileFormat'] = 'mzML file'
md['Threshold'] = 'MS-GF:QValue 0.01'
# md['Enzyme'] = 'unspecific cleavage'
# md['EnzymeSemi'] = 0
md['AnalysisSoftware'] = 'MS-GF+'
md['Software'] =  filter(None,("""
MS-GF+ Beta (v9176)
CPTAC-CDAP v1.0
CPTAC-DCC:nistcdap2psm v%s (r2281)
CPTAC-DCC:textpsm2mzid (md5:74f4441bb7507235465d0f03144c86b5)
ProteoWizard r5701
"""%VERSION).splitlines())
md['OutputFormat'] = "CPTAC-DCC:mzIdentML v1.2.1"

def generatepsms(args):
  for f in args:
   m = re.search(r'^(.*).raw.cd?ap.psm',f,re.I)
   assert m
   spectra = m.group(1)
   mzmlfile = "%s.mzML.gz"%spectra
   if opts.tsvdir:
	mzmlfile = os.path.join(opts.tsvdir,os.path.split(mzmlfile)[1])
   specmd = {}
   for nativeID,rtsec,pmz,z,pit in mzml.iterprecursor(mzmlfile):
       specmd[nativeID] = (pmz,z,rtsec,pit)
   h = open(f,'r')	
   reader = DictReader(h,dialect='excel-tab')
   for r in reader:
	if not r.get('PeptideSequence'):
	    continue
	psm = {}
	modpepstr = r['PeptideSequence']
	modpep = re.split(r'([A-Z])',modpepstr)
	pepseq = ""
	mods = []
	if modpep[0] != "":
	    pos = 0; res = "-"; delstr = modpep[0]
	    key0 = ('[',delstr)
	    key1 = ('['+modpep[1],delstr)
	    assert (key0 in cvmods or key1 in cvmods), "Keys not in cvmods: %s,%s"%(key0,key1,) + '\n' + str(r)
	    if key1 in cvmods:	
	        name,delta,spec = cvmods[key1]
		pos = 1; res = modpep[1]
	    elif key0 in cvmods:
	        name,delta,spec = cvmods[key0]
	    if spec == 'name':
	        mods.append((pos,res,delta,name))
	    else:
	        mods.append((pos,res,delta))
	j = 0
	for i in range(1,len(modpep),2):
	    j += 1
	    pepseq += modpep[i]
	    if modpep[i+1] != "":
		aa = modpep[i]
		delstr = modpep[i+1]
		name,delta,spec = cvmods[(aa,delstr)]
		if spec == 'name':
		    mods.append((j,aa,delta,name))
		else:
		    mods.append((j,aa,delta))
        mods.sort(key=lambda t: (t[0],t[2]))
	psm['Peptide'] = pepseq
	psm['Rank'] = 1
	psm['Modification'] = mods
        psm['SpectrumFile'] = os.path.split(spectra)[1]
        psm['Location'] = "%s.mzML"%(psm['SpectrumFile'],)
	psm['SpectrumID'] = "controllerType=0 controllerNumber=1 scan=%s"%(r['ScanNum'],)
        nativeID = psm['SpectrumID']
	psm['Scan'] = int(r['ScanNum'])
	psm['ExperimentalMassToCharge'] = specmd[nativeID][0]
	psm['retention time'] = (specmd[nativeID][2],"second")
        if specmd[nativeID][3]:
            psm['intensity of precursor ion'] = (specmd[nativeID][3],"number of counts")
	psm['ChargeState'] = int(r['QueryCharge'])
	for k in scores + params:
	    if k in r:
	        psm[k] = r[k]
	if 'PhospoRSPeptide' in r:
	    psm['PhosphoRSPeptide'] = r['PhospoRSPeptide']
	for k,v in r.items():
	    if k and v:
	        assert k in expected_keys, "Bad key %s in row %r"%(k,r)
	
	praccs = r['Protein'].split(';')
	psm['Protein'] = []
	ndecoy = 0; ntarget = 0;
	for pracc in praccs:
	    if pracc.startswith('XXX_'):
		ndecoy += 1
	    else:
		m = re.search(r'ref\|([^|]*)\|\(pre=(.),post=(.)\)',pracc)
	        if m:
		    psm['Protein'].append(("RefSeq:"+m.group(1),m.group(2),m.group(3)))
		ntarget += 1
	if ntarget == 0:
	    continue
	yield psm

rowlimit=1e+20
def limit(rows):
    for i,r in enumerate(rows):
	if i >= rowlimit:
	    break
	yield r

from peptidescan.OutOfCoreTable import OutOfCoreSortedTable
from peptidescan.PeptideRemapper import PeptideRemapper, UniProtIsoformAcc, RefSeqAcc, UCSCKGAcc, UniProtIsoformRefSeq, RefSeqUniProtIsoform
t = OutOfCoreSortedTable(rows=limit(generatepsms(args)),headers=headers,
			 key=lambda r: (r['SpectrumFile'],r['Scan']))
pepmaps1 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, UniProtIsoformAcc(), preprocess=True), seqdb['UniProt'])
pepmaps2 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, RefSeqAcc(), preprocess=True), seqdb['RefSeq'])
pepmaps3 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, UCSCKGAcc(), preprocess=True, translation='F'), seqdb['UCSC'])

def getdbid(f):
    sdbid = os.path.split(f)[1].rsplit('.',1)[0].split(':')
    if len(sdbid) == 1:
        return sdbid[0],sdbid[0],None,None,os.path.split(f)[1]
    elif len(sdbid) == 2:
        return ':'.join(sdbid[:2]),sdbid[0],sdbid[1],None
    return ':'.join(sdbid[:2]),sdbid[0],sdbid[1],sdbid[2],os.path.split(f)[1]

refseqorg=dict(human='H_sapiens',mouse='M_musculus')

def getdburl(name,organism,version):
    if organism == 'Human':
	taxa = 9606
    elif organism == 'Mouse':
	taxa = 10090
    else:
	raise RuntimeError("Bad organism!")
    if name == 'UniProt':
	return "http://www.uniprot.org/uniprot/?query=taxonomy%%3a%d+AND+keyword%%3a1185&force=yes&format=fasta&include=yes"%taxa
    elif name == 'RefSeq':
	return "ftp://ftp.ncbi.nlm.nih.gov/refseq/%s/mRNA_Prot/%s.protein.faa.gz"%(refseqorg[organism.lower()],
										   organism.lower())
    raise RuntimeError("Bad name!")

def getdbsrc(name,organism,version):
    if name == 'UniProt':
	return "DB source UniProt"
    elif name == 'RefSeq':
	return "DB source NCBI"
    elif name in ('Unknown','Decoy'):
	return None
    raise RuntimeError("Bad name!")

dborder = {'UniProt:Human':2,
	   'RefSeq:Human':1,
	   'UniProt:Mouse':4,
	   'RefSeq:Mouse':3,
	   'UniProt':6,
	   'RefSeq':5,
	   'Unknown':7,
	   'Decoy':8}
accprefer = RefSeqUniProtIsoform().prefer

print "MDBEGIN"
for k,v in md.items():
    if k != "Software":
        print k,v
    else:
        for v1 in md[k]:
            print k,v1
print "MDEND"


from StringIndex import StringIndex
protindex = StringIndex()
peptindex = StringIndex()
prot2pept = defaultdict(set)
dbidoutput = set()
dbidmap = dict()
for psm in t:

    dbidneeded = set()
    dbaccneeded = set(map(lambda acc: acc.rsplit('.',1)[0],map(itemgetter(0),psm['Protein'])))
    for pm in pepmaps1 + pepmaps2 + pepmaps3:
        dbid = getdbid(pm.seqdb)
	dbidmap[dbid[0]] = pm
        for pr in pm.proteins(psm['Peptide']):
            dbidneeded.add(dbid)
            try:
                dbaccneeded.remove("RefSeq:"+pr[0].rsplit('.',1)[0])
            except KeyError:
                pass
                
    if len(dbaccneeded) > 0:
        dbidneeded.add(("RefSeq","RefSeq",None,None,None))

    for dbid,name,organism,version,filename in (dbidneeded - dbidoutput):
        dbidoutput.add((dbid,name,organism,version,filename))
        print "SEQDBBEGIN"
        print "ID",dbid
        print "Name",name
        if organism:
            print "Organism",organism
        if version:
            print "Release",version
	if organism and name:
	    print "URI",getdburl(name,organism,version)
	if name and getdbsrc(name,organism,version):
            print "DBSource",getdbsrc(name,organism,version)
	if filename:
	    print "Location",filename
        print "SEQDBEND"

    print "PSMBEGIN"
    peptid = peptindex.add(psm['Peptide'])
    for key in "SpectrumFile Location FileFormat SpectrumID SpectrumIDFormat Scan Rank ChargeState ExperimentalMassToCharge Peptide".split():
	if psm.get(key) not in (None,""):
	    print key,psm[key]
    for m in psm['Modification']:
	print "Modification"," ".join(map(str,m))
    for key in scores:
	if psm.get(key) not in (None,""):
	    value = psm[key]
	    prefix = 'CPTAC-CDAP:'
	    if key == 'DeNovoScore':
		prefix = 'MS-GF:'
	    elif key == 'MSGFScore':
		prefix = 'MS-GF:'
		key = 'RawScore'
	    elif key == 'SpecEvalue':
		prefix = 'MS-GF:'
		key = 'SpecEValue'
	    elif key == 'Evalue':
		prefix = 'MS-GF:'
		key = 'EValue'
	    elif key == 'Qvalue':
		prefix = 'MS-GF:'
		key = 'QValue'
	    elif key == 'PepQvalue':
		prefix = 'MS-GF:'
		key = 'PepQValue'
	    key1 = re.sub(' ','\\ ',key)
	    print "Score",prefix+key1,value
    for key in params:
	if psm.get(key) not in (None,""):
	    key1 = re.sub(' ','\\ ',key)
	    print "Param","CPTAC-CDAP:"+key1,psm[key]
    for key in specparams:
	if psm.get(key) not in (None,""):
	    key1 = re.sub(' ','\\ ',key)
	    if isinstance(psm[key],basestring):
	        print "SpecParam",key1,psm[key]
	    else:
                units = re.sub(' ','\\ ',psm[key][1])
                print "SpecParam",key1,psm[key][0],units
                
    seen = set()
    found = False
    for pm in pepmaps1:
	dbid = getdbid(pm.seqdb)[0]
	for pracc,laa,start,end,raa,prdefline,prlen in map(itemgetter(0,2,3,4,5,6,9),pm.proteins(psm['Peptide'])):
	    # print pracc,laa,start,end,raa
	    # if "UniProt:"+pracc == psm['Protein'][0]:
	    # 	seen = True
	    # seen.add("UniProt:"+pracc)
	    print " ".join(map(str,["Protein",dbid+":"+pracc,laa[-1],raa[0],start+1,end,prlen,prdefline]))
            protid = protindex.add(dbid+":"+pracc)
	    prot2pept[protid].add(peptid)
	    found = True
    for pm in pepmaps2:
	dbid = getdbid(pm.seqdb)[0]
	for pracc,laa,start,end,raa,prdefline,prlen in map(itemgetter(0,2,3,4,5,6,9),pm.proteins(psm['Peptide'])):
	    # print pracc,laa,start,end,raa
	    # if "RefSeq:"+pracc == psm['Protein'][0]:
	    #	seen = True
	    seen.add("RefSeq:"+pracc.rsplit('.',1)[0])
	    print " ".join(map(str,["Protein",dbid+":"+pracc,laa[-1],raa[0],start+1,end,prlen,prdefline]))
            protid = protindex.add(dbid+":"+pracc)
	    prot2pept[protid].add(peptid)
	    found = True
    for pm in pepmaps3:
	dbid = getdbid(pm.seqdb)[0]
	for pracc,laa,start,end,raa,prdefline,prlen in map(itemgetter(0,2,3,4,5,6,9),pm.proteins(psm['Peptide'])):
	    print " ".join(map(str,["Protein",dbid+":"+pracc,laa[-1],raa[0],start+1,end,prlen,prdefline]))
            protid = protindex.add(dbid+":"+pracc)
	    prot2pept[protid].add(peptid)
	    found = True
    for pr in psm['Protein']:
	if pr[0].rsplit('.',1)[0] not in seen:
            print "Protein"," ".join(map(str,pr))
            protid = protindex.add(pr[0])
	    prot2pept[protid].add(peptid)
	    found = True
    assert found, "Cannot map peptide "+psm['Peptide']+" to any protein"
    print "PSMEND"

if not opts.prgrp:
    sys.exit(0)

def prprefer(prid1,prid2):
    prstr1 = protindex.string(prid1)
    prstr2 = protindex.string(prid2)
    dbid1,acc1 = prstr1.rsplit(':',1)
    dbid2,acc2 = prstr2.rsplit(':',1)
    c1 = cmp(dborder[dbid1],dborder[dbid2])
    if c1 != 0:
        return c1
    if dbid1 in ('RefSeq','Unknown','Decoy'):
	return accprefer((prid1,acc1),(prid2,acc2))
    dl1 = dbidmap[dbid1].protdb.get(acc1).defline
    dl2 = dbidmap[dbid2].protdb.get(acc2).defline
    return accprefer((prid1,dl1),(prid2,dl2))

grpmd = dict()
grpmd['AnalysisSoftware'] = 'CPTAC-DCC:nistcdap2psm'
grpmd['Threshold'] = 'no\\ threshold'
grpmd['AnalysisParams']=[]
# grpmd['AnalysisParams'] = filter(None,"""
# sequence\\ same-set\\ protein true
# sequence\\ sub-set\\ protein true
# """.splitlines())
preference = ",".join(sorted(map(itemgetter(0),dbidoutput),key=dborder.get))
grpmd['AnalysisParams'].append("same-set\\ protein\\ preference\\ order "+preference)

print "GRPMDBEGIN"
for k,v in grpmd.items():
    if k != "AnalysisParams":
        print k,v
    else:
        for v1 in grpmd[k]:
            print k,v1
print "GRPMDEND"

from parsimony import Dominator, Components
dom = Dominator(peptides=peptindex,edges=prot2pept,proteins=protindex)
dom.dominate()
dominant=set(dom.equivalentto)
comp = Components(peptides=set(peptindex),edges=prot2pept,proteins=set(protindex))
for i,(prids,pepids) in enumerate(sorted(comp,key=lambda t: -max(map(lambda prid: len(prot2pept[prid]), t[0])))):
    print "PRGRPBEGIN"
    print "Name Group%d"%(i+1,)
    print "Param - distinct\\ peptide\\ sequences",len(pepids)
    seen = set()
    for prid in sorted(prids,key=lambda prid: (1*(prid not in dominant),-len(prot2pept[prid]))):
      if prid not in seen:
        assert prid in dom.equivalentto[prid]
	for j,prid1 in enumerate(sorted(dom.equivalentto[prid],cmp=prprefer)):
	  if j == 0:
	    prid2 = prid1
	  prstr1 = protindex.string(prid1)
	  dbid1,pracc1 = prstr1.rsplit(':',1)
	  print "Protein",prstr1
          # print "Param",prstr1,"search\\ engine\\ specific\\ score",len(prot2pept[prid1])
	  if len(seen) == 0:
	    print "Param",prstr1,"anchor\\ protein"
	  elif j > 0:
	    print "Param",prstr1,"sequence\\ same-set\\ protein",protindex.string(prid2)
          print "Param",prstr1,"distinct\\ peptide\\ sequences",len(prot2pept[prid1])
	  if dbid1 not in ('RefSeq','Unknown','Decoy'):
	    print "Param",prstr1,"sequence\\ coverage",round(100.0*dbidmap[dbid1].protdb.get(pracc1).coverage(),2),"percent"
	  seen.add(prid1)
	for prid1 in set(dom.containedby[prid] - seen):
	  prstr1 = protindex.string(prid1)
	  dbid1,pracc1 = prstr1.rsplit(':',1)
	  print "Protein",prstr1
          # print "Param",prstr1,"search\\ engine\\ specific\\ score",len(prot2pept[prid1])
	  print "Param",prstr1,"sequence\\ sub-set\\ protein",protindex.string(prid2)
          print "Param",prstr1,"distinct\\ peptide\\ sequences",len(prot2pept[prid1])
	  if dbid1 not in ('RefSeq','Unknown','Decoy'):
	    print "Param",prstr1,"sequence\\ coverage",round(100.0*dbidmap[dbid1].protdb.get(pracc1).coverage(),2),"percent"
	  # print "Param",protindex.string(prid1),"non-leading\\ protein"
	  seen.add(prid1)
    print "PRGRPEND"
