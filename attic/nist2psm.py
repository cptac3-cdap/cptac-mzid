#!/usr/bin/env python
import sys, os.path, os, re, glob
from operator import itemgetter
from optparse import OptionParser
from version import VERSION
from collections import defaultdict
parser = OptionParser(version=VERSION)
opts,args = parser.parse_args()

from csv import DictReader

scores = [ 'Score', 'Dot Product', 'Prob(%)', 'Rev-Dot' ]
params = [ 'Pep', 'Nreps', 'Tfratio' ]
specparams=['retention time', 'MSPepSearch:Precursor m/z', "MSPepSearch:PrecursorMZ", "MSPepSearch:RT"]
headers = [ 'SpectrumFile', 'SpectrumID', 'SpectrumIDFormat', 'Scan', 
	    'ExperimentalMassToCharge', 'Rank', 'Peptide', 'retention time', 
	    'ChargeState', 'Modification', 'Protein' ] + scores + params + specparams


ptms = """
iTRAQ4plex           N-term 144.102063 
Acetyl               N-term   0.0
Oxidation            -        0.0
Carbamidomethyl      -        0.0
Gln->pyro-Glu        -        0.0
Glu->pyro-Glu        -        0.0
Pyro-carbamidomethyl -        0.0
"""
ptmactions = defaultdict(lambda :("-",0.0))
for l in ptms.splitlines():
    if l.strip() == "":
	continue
    name,cat,delta = l.split()
    delta = float(delta)
    ptmactions[name] = (cat,delta)

seqdb1 = glob.glob('sequence/UniProt:*.fasta')
seqdb2 = glob.glob('sequence/RefSeq:*.fasta')

def generatepsms(args):
  for f in args:
   m = re.search(r'^(.*).raw.mgf.tsv',os.path.split(f)[1],re.I)
   assert m
   spectra = m.group(1)
   specmd = {}
   h = open("%s.mzML.gz.spectrum_table.txt"%spectra,'r')
   h.next(); h.next();
   for l in h:
	sl = l.split()
	if sl[4] != 'ms2':
	    continue
	nativeID = sl[1]
	rtsec = float(sl[5])
	pmz = float(sl[12])
	specmd[nativeID] = (pmz,rtsec)
   h.close()	
   h = open(f,'r')	
   # There are three lines to ignore...
   h.next();h.next();h.next()
   reader = DictReader(h,dialect='excel-tab')
   for r in reader:
	if not r.get('Peptide'):
	    continue
	if int(r.get('Rank',1e+20)) > 1:
	    continue
	# print r
	psm = {}
	psm['Peptide'] = r['Peptide']
	kvs = re.split(r' ([\w()]+):'," "+r['Unknown'])
        kv = {}
	for i in range(1,len(kvs),2):
	  if '(' in kvs[i]:
	    kvs[i] = kvs[i].split('(',1)[0]
          kv[kvs[i]] = kvs[i+1]
        psm['SpectrumFile'] = spectra
	psm['SpectrumID'] = "controllerType=0 controllerNumber=1 scan=%s"%(kv['Scan'],)
        nativeID = "0.1.%s"%(kv['Scan'],)
        psm['SpectrumIDFormat'] = 'Thermo nativeID format'
	psm['Scan'] = int(kv['Scan'])
	psm['ExperimentalMassToCharge'] = specmd[nativeID][0]
	psm['retention time'] = (specmd[nativeID][1],"second")
        psm['MSPepSearch:Precursor m/z'] = r['Precursor m/z']
	psm['MSPepSearch:PrecursorMZ'] = kv['PrecursorMZ']
	psm['MSPepSearch:RT'] = kv['RT']
	psm['Rank'] = int(r['Rank'])
	psm['ChargeState'] = r['Charge']
	psm['Modification'] = []
	for k in scores + params:
	    psm[k] = r[k]
	mods = r['Mods'].split('/')
	for m in mods[1:]:
	    sm = m.split(',')
	    name = sm[2]
	    pos = int(sm[0])+1
	    assert name in ptmactions, "%s missing from PTM list..."%(name,)
	    ptmtype = ptmactions[name][0]
	    if ptmtype == 'N-term' and pos == 1:
		pos = 0
		res = "-"
	    elif ptmtype == 'C-term' and pos == len(r['Peptide']):
		pos = len(r['Peptide'])+1
		res = "-"
	    else:
	        res = sm[1]
	    delta = ptmactions[name][1]
	    if delta == 0.0:
	        psm['Modification'].append((pos,res,delta,name))
	    else:
	        psm['Modification'].append((pos,res,delta))
        psm['Modification'].sort(key=lambda t: (t[0],t[2]))
	flres = r['FlankRes']
	accs = r['Protein'].split()[0]
	m = re.search(r'ref\|([^|]+)',accs)
	if m:
	    psm['Protein'] = ("RefSeq:"+m.group(1),flres[0],flres[1],">"+r['Protein'])
	m = re.search(r'^(sp|tr)\|([^|]+)',accs)
	if m:
	    psm['Protein'] = ("UniProt:"+m.group(2),flres[0],flres[1],">"+r['Protein'])
	if 'Protein' not in psm:
	    print accs
	    sys.exit(1)
	yield psm

from peptidescan.OutOfCoreTable import OutOfCoreSortedTable
from peptidescan.PeptideRemapper import PeptideRemapper, UniProtIsoformAcc, RefSeqAcc
t = OutOfCoreSortedTable(rows=generatepsms(args),headers=headers,
			 key=lambda r: (r['SpectrumFile'],r['Scan']),
			 nodupkeys=True)
pepmaps1 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, UniProtIsoformAcc()), seqdb1)
pepmaps2 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, RefSeqAcc()), seqdb2)

for psm in t:
    print "PSMBEGIN"
    for key in "SpectrumFile SpectrumID SpectrumIDFormat Scan Rank ChargeState ExperimentalMassToCharge Peptide".split():
	print key,psm[key]
    for m in psm['Modification']:
	print "Modification"," ".join(map(str,m))
    for key in scores:
	if psm.get(key) != None:
	    key1 = re.sub(' ','\\ ',key)
	    print "Score","MSPepSearch:"+key1,psm[key]
    for key in params:
	if psm.get(key) != None:
	    key1 = re.sub(' ','\\ ',key)
	    print "Param","MSPepSearch:"+key1,psm[key]
    for key in specparams:
	if psm.get(key) != None:
	    key1 = re.sub(' ','\\ ',key)
	    if isinstance(psm[key],basestring):
	        print "SpecParam",key1,psm[key]
	    else:
	        print "SpecParam",key1," ".join(map(str,psm[key]))
    seen = False
    found = False
    for pm in pepmaps1:
	for pracc,laa,start,end,raa,prdefline in map(itemgetter(0,2,3,4,5,6),pm.proteins(psm['Peptide'])):
	    # print pracc,laa,start,end,raa
	    if "UniProt:"+pracc == psm['Protein'][0]:
		seen = True
	    print " ".join(map(str,["Protein","UniProt:"+pracc,laa,raa,start+1,end,prdefline]))
	    found = True
    for pm in pepmaps2:
	for pracc,laa,start,end,raa,prdefline in map(itemgetter(0,2,3,4,5,6),pm.proteins(psm['Peptide'])):
	    # print pracc,laa,start,end,raa
	    if "RefSeq:"+pracc == psm['Protein'][0]:
		seen = True
	    print " ".join(map(str,["Protein","RefSeq:"+pracc,laa,raa,start+1,end,prdefline]))
	    found = True
    # print "Protein"," ".join(map(str,psm['Protein']))
    assert found, "Cannot map peptide "+psm['Peptide']+" to any protein - claimed: " + " ".join(map(str,psm['Protein']))
    print "PSMEND"
