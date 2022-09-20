
import mzml
import psmformat
from csv import DictReader
import re, gzip, os, sys, os.path
from collections import defaultdict
from operator import itemgetter

def fopen(filename):
    if filename.lower().endswith('.gz'):
        h = gzip.open(filename)
    else:
        h = open(filename,'r')
    return h

def typecast(d):
    d1 = {}
    for k,v in list(d.items()):
        try:
            v1 = v
            v1 = float(v)
            v1 = int(v)
        except ValueError:
            pass
        d1[k.replace(' ','')] = v1
    return d1

parsers = [ 'CDAP_NISTPSM_CPTAC2', 'CDAP_NISTPSM_CPTAC2_WITH_DECOYS', 'PepArML', 'UMich_CPTAC3_PGDAC', 'UMich_CPTAC3_PGDAC_WITH_DECOYS']

class CDAP_NISTPSM_CPTAC2(object):

    # Note that this extracter as written assumes Thermo raw files as the
    # original data-source...

    # Expected headers, tsv file
    #
    # FileName ScanNum QueryPrecursorMz OriginalPrecursorMz
    # PrecursorError(ppm) QueryCharge OriginalCharge PrecursorScanNum
    # PrecursorArea PrecursorRelAb RTAtPrecursorHalfElution PeptideSequence
    # AmbiguousMatch Protein DeNovoScore MSGFScore Evalue Qvalue PepQvalue
    # PrecursorPurity FractionDecomposition HCDEnergy iTRAQ114 iTRAQ115
    # iTRAQ116 iTRAQ117 iTRAQFlags iTRAQTotalAb iTRAQFractionOfTotalAb
    # PhosphoRSPeptide nPhospho FullyLocalized
    keepdecoys = False
    ignore = """
    ScanNum FileName Protein QueryCharge PeptideSequence OriginalCharge
    OriginalPrecursorMz PrecursorScanNum PhospoRSPeptide TMTFlags
    """.split()
    _scores = """
    DeNovoScore MSGFScore SpecEvalue Evalue Qvalue PepQvalue
    """.split()
    _params = """
    QueryPrecursorMz
    PrecursorError(ppm)
    PrecursorArea PrecursorRelAb RTAtPrecursorHalfElution
    AmbiguousMatch PrecursorPurity FractionDecomposition HCDEnergy
    iTRAQ114 iTRAQ115 iTRAQ116 iTRAQ117
    iTRAQFlags iTRAQTotalAb iTRAQFractionOfTotalAb
    TMT10-126 TMT10-127N TMT10-127C TMT10-128N TMT10-128C
    TMT10-129N TMT10-129C TMT10-130N TMT10-130C TMT10-131
    TMT10-Flags TMT10-FractionOfTotalAb TMT10-TotalAb
    TMT6-126 TMT6-127 TMT6-128 TMT6-129 TMT6-130 TMT6-131
    TMT6-Flags TMT6-FractionOfTotalAb TMT6-TotalAb
    TMT11-126C TMT11-127N TMT11-127C TMT11-128N TMT11-128C
    TMT11-129N TMT11-129C TMT11-130N TMT11-130C TMT11-131N TMT11-131C
    TMT11-Flags TMT11-FractionOfTotalAb TMT11-TotalAb
    TMT16-126C TMT16-127N TMT16-127C TMT16-128N TMT16-128C
    TMT16-129N TMT16-129C TMT16-130N TMT16-130C TMT16-131N
    TMT16-131C TMT16-132N TMT16-132C TMT16-133N TMT16-133C TMT16-134N
    TMT16-Flags TMT16-FractionOfTotalAb TMT16-TotalAb
    TMT18-126C TMT18-127N TMT18-127C TMT18-128N TMT18-128C
    TMT18-129N TMT18-129C TMT18-130N TMT18-130C TMT18-131N
    TMT18-131C TMT18-132N TMT18-132C TMT18-133N TMT18-133C
    TMT18-134N TMT18-134C TMT18-135N
    TMT18-Flags TMT18-FractionOfTotalAb TMT18-TotalAb
    PhosphoRSPeptide nPhospho FullyLocalized
    """.split()
    expected_keys = ignore + _scores + _params
    ptms = """
    [   +304.207        UNIMOD:2016     304.207146
    K   +304.207        UNIMOD:2016     304.207146
    [   +229.163        UNIMOD:737      229.162932
    K   +229.163        UNIMOD:737      229.162932
    [   +144.102        UNIMOD:214      144.102063
    K   +144.102        UNIMOD:214      144.102063
    M   +15.995         Oxidation       15.994915
    C   +57.021         Carbamidomethyl 57.021464
    Q   +0.984          Deamidated      0.984016
    N   +0.984          Deamidated      0.984016
    [Q  -17.027         Gln->pyro-Glu   -17.026549
    [Q  -17.026         Gln->pyro-Glu   -17.026549
    [E  -18.011         Glu->pyro-Glu   -18.010565
    S   +79.966         Phospho         79.966331
    T   +79.966         Phospho         79.966331
    Y   +79.966         Phospho         79.966331
    [   +42.011         Acetyl          42.010565
    K   +42.011         Acetyl          42.010565
    K   +114.043        GG              114.042927
    """
    keymap = {
        'DeNovoScore': ('MS-GF','DeNovoScore'),
        'MSGFScore':   ('MS-GF','RawScore'),
        'SpecEvalue':  ('MS-GF','SpecEValue'),
        'Evalue':      ('MS-GF','EValue'),
        'Qvalue':      ('MS-GF','QValue'),
        'PepQvalue':   ('MS-GF','PepQValue')
    }
    def __init__(self,filenames,**kwargs):
        self.filter = kwargs.get('filter')
        self.cdapversion = kwargs.get('cdapversion','1.1')
        self.filenames = filenames
        self.cvmods = {}
        for l in self.ptms.splitlines():
            if l.strip() == "":
                continue
            aa,delstr,name,delta = l.split()
            delta = float(delta)
            self.cvmods[(aa,delstr)] = (name,delta)
        self.scores = []
        for s in self._scores:
            if s in self.keymap:
                self.scores.append(':'.join(self.keymap[s]))
            else:
                self.scores.append('CPTAC-CDAP:'+s)
        self.params = []
        for p in self._params:
            if p in self.keymap:
                self.params.append(':'.join(self.keymap[p]))
            else:
                self.params.append('CPTAC-CDAP:'+p)

    def metadata(self):
        md = dict()
        md['SpectrumIDFormat'] = 'Thermo nativeID format'
        md['Threshold'] = 'MS-GF:QValue 0.01'
        md['AnalysisSoftware'] = 'MS-GF+'
        md['OutputFormat'] = "CPTAC-DCC:mzIdentML v1.2.2"
        software = """
             MS-GF+ Release (v2017.01.27) (27 Jan 2017)
             CPTAC-CDAP v%s
        """%(self.cdapversion,)
        md['Software'] =  [_f for _f in map(str.strip,(software.splitlines())) if _f]
        return md

    def seqdbs(self,proteins=None):
        if proteins:
            return list(map(itemgetter(0),proteins))
        if self.keepdecoys:
            return ["RefSeq","UniProt","Decoy"]
        return ["RefSeq","UniProt"]

    def removeProtein(self,sdb,pr,proteins):

        # if sdb.name() != 'RefSeq':
        #     return proteins

        retval = []
        for pri in proteins:
            if pr[1] == pri[1]:
                continue
            pracc = pr[1].rsplit('.',1)[0]
            priacc = pri[1].rsplit('.',1)[0]
            if pracc == priacc:
                continue
            retval.append(pri)

        return retval

    def extract_mods(self,r):
        modpep = re.split(r'([A-Z])',r['PeptideSequence'])
        pepseq = ""
        mods = []
        if modpep[0] != "":
            pos = 0; res = "-"; delstr = modpep[0]
            key0 = ('[',delstr)
            key1 = ('['+modpep[1],delstr)
            # assert (key0 in self.cvmods or key1 in self.cvmods), \
            #        "Keys not in cvmods: %s,%s"%(key0,key1,) + '\n' + str(r)
            if key1 in self.cvmods:
                name,delta = self.cvmods[key1]
                pos = 1; res = modpep[1]
            elif key0 in self.cvmods:
                name,delta = self.cvmods[key0]
            else:
                delta = float(delstr)
                name = ""
            mods.append((pos,res,delta,name))
        j = 0
        for i in range(1,len(modpep),2):
            j += 1
            pepseq += modpep[i]
            if modpep[i+1] != "":
                aa = modpep[i]
                delstr = modpep[i+1]
                if (aa,delstr) in self.cvmods:
                    name,delta = self.cvmods[(aa,delstr)]
                else:
                    delta = float(delstr)
                    name = ""
                mods.append((j,aa,delta,name))
        mods.sort(key=lambda t: (t[0],t[2]))
        return pepseq,mods

    def extract_specfile(self,r):
        m = re.search(r'^(.*)\.(raw|mzml)',r['FileName'],re.I)
        assert m, r['FileName']
        return os.path.split(m.group(1))[1]

    def decoyPrefix(self):
        return "XXX_"

    def psms(self):
        for filename in self.filenames:
            reader = DictReader(fopen(filename),dialect='excel-tab')
            for r in reader:
                if not r.get('PeptideSequence'):
                    continue
                if self.filter != None and not eval(self.filter,None,typecast(r)):
                    continue
                psm = {}
                spectra = self.extract_specfile(r)
                pepseq,mods = self.extract_mods(r)
                psm['Peptide'] = pepseq
                psm['Rank'] = 1
                psm['Modification'] = mods
                psm['SpectrumFile'] = spectra
                psm['_fullpath'] = filename
                psm['Location'] = "%s.mzML"%(spectra,)
                psm['SpectrumID'] = "controllerType=0 controllerNumber=1 scan=%s"%(r['ScanNum'],)
                psm['Scan'] = int(r['ScanNum'])
                psm['ChargeState'] = int(r['QueryCharge'])
                if 'PhospoRSPeptide' in r:
                    psm['PhosphoRSPeptide'] = r['PhospoRSPeptide']
                if 'TMTFlags' in r and 'TMT10-126' in r:
                    r['TMT10-Flags'] = r['TMTFlags']
                for k,v in list(r.items()):
                    if k and v:
                        assert k in self.expected_keys, "Bad key %s in row %r"%(k,r)
                for k in self._scores + self._params:
                    if k in r:
                        if k in self.keymap:
                            psm[":".join(self.keymap[k])] = r[k]
                        else:
                            psm["CPTAC-CDAP:"+k] = r[k]

                praccs = r['Protein'].split(';')
                psm['Protein'] = []
                ndecoy = 0; ntarget = 0;
                decoy_prefix = 'XXX_'
                for pracc in praccs:
                    isdecoy = False
                    if pracc.startswith(decoy_prefix):
                        isdecoy = True
                        pracc = pracc[len(decoy_prefix):]
                    m = re.search(r'^(.*)\(pre=(.),post=(.)\)$',pracc)
                    pracc = m.group(1)
                    laa = m.group(2)
                    raa = m.group(3)

                    source = None;

                    m = re.search(r'^(gi\|\d+\|)?ref\|([^|]*)\|$',pracc)
                    if m:
                        source = "RefSeq"; pracc = m.group(2)

                    m = re.search(r'^([NXYA]P_[^|]*)$',pracc)
                    if m:
                        source = "RefSeq"; pracc = m.group(1)

                    m = re.search(r'^(sp|tr)\|([^|]*)\|([^|]*)$',pracc)
                    if m:
                        source = "UniProt"; pracc = m.group(2)

                    if not source:
                        raise RuntimeError

                    if not isdecoy:
                        psm['Protein'].append((source,pracc,laa,raa))
                        ntarget += 1
                    elif self.keepdecoys:
                        psm['Protein'].append(("Decoy",decoy_prefix+pracc,laa,raa))
                        ndecoy += 1

                if self.keepdecoys:
                    if ndecoy > 0:
                        psm['decoy'] = True
                    else:
                        psm['decoy'] = False
                else:
                    if ntarget == 0:
                        continue

                yield psm

    def __iter__(self):
        return self.psms()


class CDAP_NISTPSM_CPTAC2_WITH_DECOYS(CDAP_NISTPSM_CPTAC2):
    keepdecoys = True

class PepArML(object):

    # spectra_set,start_scan,end_scan,precursor_mz,assumed_charge,precursor_neutral_mass,calc_neutral_pep_mass,massdiff,leftaa,peptide,rightaa,mods,protein,decoy,nagree,estfdr,DeNovoScore-g1,EValue-g1,RawScore-g1,SpecEValue-g1,eval-g1,rank-g1,eval-o1,mzhits-o1,pval-o1,rank-o1,b_ions-t1,b_score-t1,eval-t1,fI-t1,hyperscore-t1,maxI-t1,pval-t1,rank-t1,sumI-t1,y_ions-t1,y_score-t1,peptide_len,precursor_intensity,activation_method,retention_time,basepeak_intensity,total_ion_current,c13massdiff,c13peak,c13relint,icscore,hydrophobicity-kd,gravy,hydrophilicity-hw,isoelectric-point,PPP0008,PPP0160,PPP0161,PPP0180,PPP0196,PPP0197,PPP0313,PPP0502,PPP0787,PPP0788,PPP0823,PPP0877,PPP1006,mapped,ntspec,ctspec,nspect,misscl,rtpred,rtdelta,siteprob,pred,pepfdr,decoypeps,2udpr,3udpr

    _scores = """
    nagree estfdr pred eval-g1 eval-o1 eval-t1 eval-s1 eval-m1 eval-y1 eval-c1 eval-k1
    """.split()

    _params = """
    """.split()

    ptms = """
    [   +144.102        UNIMOD:214      144.102063
    K   +144.102        UNIMOD:214      144.102063
    M   +15.995         Oxidation       15.994915
    C   +57.021         Carbamidomethyl 57.021464
    Q   +0.984          Deamidated      0.984016
    N   +0.984          Deamidated      0.984016
    [Q  -17.027         Gln->pyro-Glu   -17.026549
    [Q  -17.026         Gln->pyro-Glu   -17.026549
    [C  -17.027         Ammonia-loss    -17.026549
    [E  -18.011         Glu->pyro-Glu   -18.010565
    S   +79.966         Phospho         79.966331
    T   +79.966         Phospho         79.966331
    Y   +79.966         Phospho         79.966331
    [   +42.011         Acetyl          42.010565
    K   +42.011         Acetyl          42.010565
    K   +114.043        GG              114.042927
    """

    keymap = {
        'eval-g1':     ('MS-GF','SpecEValue'),
        'eval-o1':     ('OMSSA','evalue'),
        'eval-t1':     ('X!Tandem','expect'),
        'eval-m1':     ('Mascot','Mascot:expectation value'),
        'eval-k1':     ('X!Tandem','k-score plugin expect'),
        'eval-s1':     ('X!Tandem','s-score plugin expect'),
        'eval-c1':     ('Comet','expectation value'),
        'eval-y1':     ('MyriMatch','expect'),
        'nagree':      ('PepArML','Agree'),
        'estfdr':      ('PepArML','QValue'),
        'pred':        ('PepArML','Confidence'),
        'siteprob':    ('PepArML','SiteLocalization'),
    }

    def __init__(self,filenames,**kwargs):
        self.filter = kwargs.get('filter')
        self.filenames = filenames
        self.cvmods = {}
        for l in self.ptms.splitlines():
            if l.strip() == "":
                continue
            aa,delstr,name,delta = l.split()
            delta = float(delta)
            self.cvmods[(aa,delstr)] = (name,delta)
        self.scores = []
        for s in self._scores:
            if s in self.keymap:
                self.scores.append(':'.join(self.keymap[s]))
            else:
                self.scores.append('PepArML:'+s)
        self.params = []
        for p in self._params:
            if p in self.keymap:
                self.params.append(':'.join(self.keymap[p]))
            else:
                self.params.append('PepArML:'+p)

    def metadata(self):
        md = dict()
        md['SpectrumIDFormat'] = 'Thermo nativeID format'
        md['AnalysisSoftware'] = 'PepArML'
        md['Threshold'] = 'PepArML:QValue 0.0031'
        md['Software'] =  [_f for _f in map(str.strip,("""
        PepArML http://grg.tn/PepArML
        """.splitlines())) if _f]
        return md

    def seqdbs(self,proteins=None):
        return ["UniProt"]

    def removeProtein(self,sdb,pr,proteins):

        if sdb.name() != 'UniProt':
            return proteins

        retval = []
        for pri in proteins:
            if pr[1] == pri[1]:
                continue
            retval.append(pri)

        return retval

    def extract_mods(self,r):
        if r['mods'] == "-":
            return []
        seq=r['peptide']
        mods=[]
        for m in r['mods'].split(','):
            match = re.search(r'^(.)(\d+):(.*)$',m)
            assert match
            aa = match.group(1)
            pos = int(match.group(2))
            deltastr = match.group(3)
            key0 = (aa,deltastr)
            key1 = None
            if pos == 0:
                key1 = (aa+seq[0],deltastr)
            elif pos == 1:
                key1 = ('['+aa,deltastr)
            if key0 in self.cvmods:
                name,delta = self.cvmods[key0]
            elif key1 in self.cvmods:
                name,delta = self.cvmods[key1]
                pos = 1; aa = seq[0]
            else:
                raise RuntimeError("keys not in cvmods: %s,%s, mods = %r"%(key0,key1,m))
            mods.append((pos,aa,delta,name))
        mods.sort(key=lambda t: (t[0],t[2]))
        return mods

    def psms(self):
        for filename in self.filenames:
            reader = DictReader(fopen(filename))
            for r in reader:
                if not r.get('peptide'):
                    continue
                if r.get('decoy') == '1':
                    continue
                if self.filter != None and not eval(self.filter,None,typecast(r)):
                    continue
                psm = {}
                spectra = r['spectra_set']
                psm['Peptide'] = r['peptide']
                psm['Rank'] = 1
                psm['Modification'] = self.extract_mods(r)
                psm['SpectrumFile'] = spectra
                psm['_fullpath'] = filename
                psm['Location'] = "%s.mzML"%(spectra,)
                psm['SpectrumID'] = "controllerType=0 controllerNumber=1 scan=%s"%(r['start_scan'],)
                psm['Scan'] = int(r['start_scan'])
                psm['ChargeState'] = int(r['assumed_charge'])
                for k in self._scores + self._params:
                    if r.get(k,"").strip():
                        if k in self.keymap:
                            psm[":".join(self.keymap[k])] = r[k].strip()
                        else:
                            psm["PepArML:"+k] = r[k].strip()
                praccs = r['protein'].split(';')
                psm['Protein'] = []
                for pracc in praccs:
                    psm['Protein'].append(("UniProt",pracc))
                yield psm

    def __iter__(self):
        return self.psms()

class UMich_CPTAC3_PGDAC(object):

    # Spectrum Peptide Modified Peptide Charge Retention Calculated
    # M/Z Observed M/Z Original Delta Mass Adjusted Delta Mass
    # Experimental Mass Peptide Mass Expectation Hyperscore Nextscore
    # PeptideProphet Probability Intensity Is Unique Is Used Assigned
    # Modifications Observed Modifications Number of Phospho Sites
    # Phospho Site Localization Protein Protein ID Entry Name Gene
    # Protein Description Mapped Proteins Purity 7316-321 Abundance
    # 7316-341 Abundance 7316-3069 Abundance 7316-466 Abundance
    # 7316-475 Abundance 7316-487 Abundance 7316-206 Abundance
    # 7316-1851 Abundance 7316-3319 Abundance 7316-496 Abundance
    # Bridge-1 Abundance

    keepdecoys = False

    _scores = list(map(str.strip,[_f for _f in """
Expectation
Hyperscore
Nextscore
PeptideProphet Probability
probability
    """.splitlines() if _f]))
    _params = list(map(str.strip,[_f for _f in """
Number of Phospho Sites
Phospho Site Localization
TMT10-126
TMT10-127N
TMT10-127C
TMT10-128N
TMT10-128C
TMT10-129N
TMT10-129C
TMT10-130N
TMT10-130C
TMT10-131
TMT11-126C
TMT11-127N
TMT11-127C
TMT11-128N
TMT11-128C
TMT11-129N
TMT11-129C
TMT11-130N
TMT11-130C
TMT11-131N
TMT11-131C
    """.splitlines() if _f]))
    keymap = {
      'probability': ('probability',),
    }
    ptms = """
    [   +229.163        UNIMOD:737      229.162932
    K   +229.163        UNIMOD:737      229.162932
    [   +144.102        UNIMOD:214      144.102063
    K   +144.102        UNIMOD:214      144.102063
    M   +15.995         Oxidation       15.994915
    C   +57.021         Carbamidomethyl 57.021464
    C   +57.022         Carbamidomethyl 57.021464
    Q   +0.984          Deamidated      0.984016
    N   +0.984          Deamidated      0.984016
    [Q  -17.027         Gln->pyro-Glu   -17.026549
    [Q  -17.026         Gln->pyro-Glu   -17.026549
    [E  -18.011         Glu->pyro-Glu   -18.010565
    S   +79.966         Phospho         79.966331
    T   +79.966         Phospho         79.966331
    Y   +79.966         Phospho         79.966331
    [   +42.011         Acetyl          42.010565
    K   +42.011         Acetyl          42.010565
    K   +114.043        GG              114.042927
    """

    def __init__(self,filenames,**kwargs):
        self.filter = kwargs.get('filter')
        self.filenames = filenames
        self.cvmods = {}
        for l in self.ptms.splitlines():
            if l.strip() == "":
                continue
            aa,delstr,name,delta = l.split()
            delta = float(delta)
            self.cvmods[(aa,delstr)] = (name,delta)
        self.scores = []
        for s in self._scores:
            if s in self.keymap:
                self.scores.append(':'.join(self.keymap[s]))
            else:
                self.scores.append('UMich:'+s)
        self.params = []
        for p in self._params:
            if p in self.keymap:
                self.params.append(':'.join(self.keymap[p]))
            else:
                self.params.append('UMich:'+p)

    def metadata(self):
        md = dict()
        md['SpectrumIDFormat'] = 'Thermo nativeID format'
        md['Threshold'] = ""
        md['AnalysisSoftware'] = ""
        md['OutputFormat'] = "CPTAC-DCC:mzIdentML v1.2.2"
        md['Software'] =  [_f for _f in map(str.strip,("""
        """.splitlines())) if _f]
        return md

    def seqdbs(self,proteins=None):
        if proteins:
            return list(map(itemgetter(0),proteins))
        if self.keepdecoys:
            return ["RefSeq","UniProt","Decoy"]
        return ["RefSeq","UniProt"]

    def decoyPrefix(self):
        return "XXX_"

    def removeProtein(self,sdb,pr,proteins):

        # if sdb.name() != 'RefSeq':
        #     return proteins

        retval = []
        for pri in proteins:
            if pr[1] == pri[1]:
                continue
            pracc = pr[1].rsplit('.',1)[0]
            priacc = pri[1].rsplit('.',1)[0]
            if pracc == priacc:
                continue
            retval.append(pri)

        return retval

    def extract_specfile(self,r):
        m = re.search(r'^(.*)\.(\d+)\.(\d+)\.\d+$',r['Spectrum'])
        assert m, r['Spectrum']
        assert m.group(2) == m.group(3), r['Spectrum']
        return os.path.split(m.group(1))[1],int(m.group(2))

    def extract_mods(self,r):
        modstr = r['Assigned Modifications']
        mods = []
        if modstr:
            for mod in map(str.strip,modstr.split(',')):
                m = re.search(r'^(\d*)(.)\((.*)\)$',mod)
                assert m, r
                aa = m.group(2)
                deltastr = "%+.3f"%(float(m.group(3)),)
                if aa == "n":
                    aa = '-'
                    pos = 0
                    key = ('[',deltastr)
                else:
                    pos = int(m.group(1))
                    key = (aa,deltastr)
                if key in self.cvmods:
                    name,delta = self.cvmods[key]
                else:
                    raise RuntimeError("key not in cvmods: %s, mods = %r"%(key,mod))
                mods.append((pos,aa,delta,name))
            mods.sort(key=lambda t: (t[0],t[2]))
        return mods

    def psms(self):
        for filename in self.filenames:
            reader = DictReader(fopen(filename),dialect='excel-tab')
            for r in reader:
                if self.filter != None and not eval(self.filter,None,typecast(r)):
                    continue
                psm = {}
                spectra,scan = self.extract_specfile(r)
                pepseq = r['Peptide']
                mods = self.extract_mods(r)
                psm['Peptide'] = pepseq
                psm['Rank'] = 1
                psm['Modification'] = mods
                psm['SpectrumFile'] = spectra
                psm['_fullpath'] = filename
                psm['Location'] = "%s.mzML"%(spectra,)
                psm['SpectrumID'] = "controllerType=0 controllerNumber=1 scan=%s"%(scan,)
                psm['Scan'] = scan
                psm['ChargeState'] = int(r['Charge'])
                psm['probability'] = float(r['PeptideProphet Probability'])

                for k in self._scores + self._params:
                    if k in r:
                        if k in self.keymap:
                            psm[":".join(self.keymap[k])] = r[k]
                        else:
                            psm["UMich:"+k] = r[k]

                quant = []
                for h in reader.fieldnames:
                    if h.endswith(' Abundance'):
                        quant.append(h)
                labels = None
                if len(quant) == 10:
                    labels = """TMT10-126 TMT10-127N TMT10-127C TMT10-128N TMT10-128C
                                TMT10-129N TMT10-129C TMT10-130N TMT10-130C TMT10-131""".split()
                elif len(quant) == 11:
                    labels = """TMT11-126C TMT11-127N TMT11-127C TMT11-128N TMT11-128C
                                TMT11-129N TMT11-129C TMT11-130N TMT11-130C TMT11-131N TMT11-131C""".split()
                if labels:
                    assert len(labels) == len(quant)
                    for l,h in zip(labels,quant):
                        psm['UMich:'+l] = r[h]
                else:
                    if len(quant) > 0:
                        raise RuntimeError

                praccs = [r['Protein']]
                praccs.extend(list(map(str.strip,r['Mapped Proteins'].split(','))))
                ndecoy = 0; ntarget = 0;
                decoy_prefix = 'rev_'
                psm['Protein'] = []
                for pracc in praccs:

                    if not pracc:
                        continue

                    isdecoy = False
                    if pracc.startswith(decoy_prefix):
                        isdecoy = True
                        pracc = pracc[len(decoy_prefix):]

                    source = None;

                    m = re.search(r'^([NXYA]P_[^|]*)$',pracc)
                    if m:
                        source = "RefSeq"; pracc = m.group(1)

                    m = re.search(r'^(sp|tr)\|([^|]*)\|([^|]*)$',pracc)
                    if m:
                        source = "UniProt"; pracc = m.group(2)

                    if not source:
                        raise RuntimeError

                    if not isdecoy:
                        psm['Protein'].append((source,pracc,None,None))
                        ntarget += 1
                    elif self.keepdecoys:
                        psm['Protein'].append(("Decoy","XXX_"+pracc,None,None))
                        ndecoy += 1

                if self.keepdecoys:
                    if ndecoy > 0:
                        psm['decoy'] = True
                    else:
                        psm['decoy'] = False
                else:
                    if ntarget == 0:
                        continue

                yield psm

    def __iter__(self):
        return self.psms()

class UMich_CPTAC3_PGDAC_WITH_DECOYS(UMich_CPTAC3_PGDAC):
    keepdecoys = True

class AddMZMLFields(object):

    specparams = ['intensity of precursor ion','retention time']

    def __init__(self,input,specdir=None):
        self.input = input
        self.specdir = specdir

    def psms(self):
        currentmzmlfile = None
        for psm in self.input.psms():
            mzmlfilename = psm['Location']
            if mzmlfilename != currentmzmlfile:
                if not self.specdir:
                    filename = os.path.join(os.path.split(psm['_fullpath'])[0],mzmlfilename)
                else:
                    filename = os.path.join(self.specdir,mzmlfilename)
                if not os.path.exists(filename):
                    filename += '.gz'
                assert os.path.exists(filename)
                currentmzmlfile = mzmlfilename
                specmd = {}
                for nativeID,rtsec,pmz,z,pit in mzml.iterprecursor(filename):
                    specmd[nativeID] = (pmz,z,rtsec,pit)
            nativeID = psm['SpectrumID']
            psm['ExperimentalMassToCharge'] = specmd[nativeID][0]
            psm['retention time'] = (specmd[nativeID][2],"second")
            if specmd[nativeID][3]:
                psm['intensity of precursor ion'] = (specmd[nativeID][3],"number of counts")
            yield psm

    def __iter__(self):
        return self.psms()

class AddReporterIonFields(object):

    def __init__(self,input,labels,specdir=None):
        self.input = input
        self.labels = labels
        self.specdir = specdir

    def psms(self):
        currentmzmlfile = None
        for psm in self.input.psms():
            mzmlfilename = psm['Location']
            if mzmlfilename != currentmzmlfile:
                if not self.specdir:
                    filename = os.path.join(os.path.split(psm['_fullpath'])[0],mzmlfilename)
                else:
                    filename = os.path.join(self.specdir,mzmlfilename)
                if not os.path.exists(filename):
                    filename += '.gz'
                assert os.path.exists(filename)
                currentmzmlfile = mzmlfilename
                specmd = {}
                for nativeID,reporters in mzml.iterreporters(filename,self.labels):
                    specmd[nativeID] = reporters
            nativeID = psm['SpectrumID']
            for k,v in specmd[nativeID].items():
                if not k.startswith("_"):
                    psm["CPTAC-CDAP:"+self.labels+"-"+k] = "%.6g/%s"%(v[0],"%.2f"%v[2] if v[2] != None else '?')
                elif k == "_total":
                    psm["CPTAC-CDAP:"+self.labels+"-TotalAb"] = "%.6g"%v
                elif k == "_frac":
                    psm["CPTAC-CDAP:"+self.labels+"-FractionOfTotalAb"] = "%.6g"%v
            yield psm

    def __iter__(self):
        return self.psms()

def getParser(fmt):
    if fmt in parsers:
        return eval(fmt)
    raise RuntimeError("Bad PSM parser: %s"%fmt)

def getParsers():
    return parsers
