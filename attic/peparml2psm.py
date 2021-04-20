#!/usr/bin/env python
import sys, os.path, os, re, glob
from operator import itemgetter
from optparse import OptionParser
from version import VERSION
from collections import defaultdict
parser = OptionParser(version=VERSION)
parser.add_option('--spectraset',type='str',dest='spectraset',default=None)
opts,args = parser.parse_args()

from csv import DictReader

scores = [ 'estfdr', 'pepfdr', 'nagree', 'pred', 'Mascot:expectation value', 'OMSSA:evalue', 'X!Tandem:expect', 'K-Score:expect']
params = [ ]
specparams=['retention time']
headers = [ 'SpectrumFile', 'SpectrumID', 'SpectrumIDFormat', 'FileFormat', 'Location', 'Scan', 
	    'ExperimentalMassToCharge', 'Rank', 'Peptide', 'retention time', 'passThreshold',
	    'ChargeState', 'Modification', 'Protein' ] + scores + params + specparams


ptms = """
M	+15.994	Oxidation            -        0.0
M	+15.995	Oxidation            -        0.0
C	+57.021	Carbamidomethyl      -        0.0
Q	-17.027 Gln->pyro-Glu        -        0.0
E	-18.011 Glu->pyro-Glu        -        0.0
C	-17.027 Pyro-carbamidomethyl -        0.0
S	+79.966 Phospho		     -	      0.0
T	+79.966 Phospho		     -	      0.0
Y	+79.966 Phospho		     -	      0.0
"""
ptmactions = defaultdict(lambda :(None,"-",0.0))
for l in ptms.splitlines():
    if l.strip() == "":
	continue
    res,delstr,name,cat,delta = l.split()
    delta = float(delta)
    ptmactions[(res,delstr)] = (name,cat,delta)

# seqdb1 = []
seqdb2 = []
seqdb1 = glob.glob('sequence/UniProt:Rat.fasta')
# seqdb2 = glob.glob('sequence/RefSeq:Rat.fasta')

def generatepsms(args):
  for f in args:
     h = open(f,'r')	
     reader = DictReader(h)
     for r in reader:
	if opts.spectraset and r['spectra_set'] != opts.spectraset:
	    continue
	psm = {}
	psm['Peptide'] = r['peptide']
        psm['SpectrumFile'] = r['spectra_set']
        assert(r['start_scan'] == r['end_scan'])
	psm['SpectrumID'] = "mzMLid=controllerType=0 controllerNumber=1 scan=%s"%(r['start_scan'],)
        # nativeID = "%s"%(r['start_scan'],)
        psm['SpectrumIDFormat'] = 'mzML unique identifier'
        psm['FileFormat'] = 'mzML file'
        psm['Location'] = r['spectra_set']+'.mzML'
	psm['Scan'] = int(r['start_scan'])
	psm['ExperimentalMassToCharge'] = r['precursor_mz']
	if 'rentention_time' in r:
	     psm['retention time'] = (r['retention_time'],"second")
	psm['Rank'] = 1
	psm['ChargeState'] = r['assumed_charge']
        for k in scores + params:
	    if k in r:
                psm[k] = r[k]
        psm['passThreshold'] = 'false'
        if float(r['estfdr']) < 0.1:
            psm['passThreshold'] = 'true'
        for k in r:
            m = re.search('^eval-(.)\d+$',k)
            if m and r[k].strip():
                engine = m.group(1)
                if engine == 'm':
                    psm['Mascot:expectation value'] = float(r[k])
                elif engine == 'o':
                    psm['OMSSA:evalue'] = float(r[k])
                elif engine == 't':
                    psm['X!Tandem:expect'] = float(r[k])
                elif engine ==  'k':
                    psm['K-Score:expect'] = float(r[k])
	psm['Modification'] = []
	mods = r['mods'].split(',')
	for m in mods:
	    if m == "-":
		break
	    aapos,delta = m.split(':')
	    delta = float(delta)
	    aa = aapos[0]
	    pos = int(aapos[1:])
	    key = (aa,"%+.3f"%delta)
	    assert key in ptmactions, "%s missing from PTM list..."%(key,)
	    name = ptmactions[key][0]
	    ptmtype = ptmactions[key][1]
	    if ptmtype == 'N-term':
		pos = 0
		res = "-"
	    elif ptmtype == 'C-term':
		pos = len(r['Peptide'])+1
		res = "-"
	    else:
	        res = aa
	    delta = ptmactions[key][2]
	    if delta == 0.0:
	        psm['Modification'].append((pos,res,delta,name))
	    else:
	        psm['Modification'].append((pos,res,delta))
        psm['Modification'].sort(key=lambda t: (t[0],t[2]))
        deli = -1
        for i,m in enumerate(psm['Modification']):
            if len(m) == 4 and m[3] == 'Pyro-carbamidomethyl':
                deli = i+1
                break
        if deli >= 0 and len(psm['Modification'][deli]) == 4 and psm['Modification'][deli][3] == 'Carbamidomethyl':
            del psm['Modification'][deli]
        psm['Protein'] = r['protein'].split(';')
	yield psm

from peptidescan.OutOfCoreTable import OutOfCoreSortedTable
from peptidescan.PeptideRemapper import PeptideRemapper, UniProtIsoformAcc, RefSeqAcc
t = OutOfCoreSortedTable(rows=generatepsms(args),headers=headers,
			 key=lambda r: (r['SpectrumFile'],r['Scan']),
			 nodupkeys=True)
pepmaps1 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, UniProtIsoformAcc()), seqdb1)
pepmaps2 = map(lambda sdb: PeptideRemapper(map(lambda r: r['Peptide'], t), sdb, RefSeqAcc()), seqdb2)

for psm in t:
    found = False
    if not found:
      for pm in pepmaps1:
	for pracc,laa,start,end,raa,prdefline in map(itemgetter(0,2,3,4,5,6),pm.proteins(psm['Peptide'])):
	    found = True
	    break
    if not found:
      for pm in pepmaps2:
	for pracc,laa,start,end,raa,prdefline in map(itemgetter(0,2,3,4,5,6),pm.proteins(psm['Peptide'])):
	    found = True
	    break
    if not found:
	continue
    print "PSMBEGIN"
    for key in "SpectrumFile SpectrumID SpectrumIDFormat FileFormat Location Scan Rank ChargeState passThreshold ExperimentalMassToCharge Peptide".split():
	print key,psm[key]
    for m in psm['Modification']:
	print "Modification"," ".join(map(str,m))
    for key in scores:
	if psm.get(key) != None:
	    key1 = re.sub(' ','\\ ',key)
            if ':' in key1:
                print "Score",key1,psm[key]
            else:
                print "Score","PepArML:"+key1,psm[key]
    for key in params:
	if psm.get(key) != None:
	    key1 = re.sub(' ','\\ ',key)
            if ':' in key1:
                print "Param",key1,psm[key]
            else:
                print "Param","PepArML:"+key1,psm[key]
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
    if not found:
        print >>sys.stderr, "Problem peptide:",psm['Peptide']
        for pr in psm['Protein']:
            print "Protein","UniProt:"+pr
    # print "Protein"," ".join(map(str,psm['Protein']))
    # assert found, "Cannot map peptide "+psm['Peptide']+" to any protein - claimed: " + " ".join(map(str,psm['Protein']))
    print "PSMEND"
