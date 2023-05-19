
import os, sys, os.path, re

from peptidescan.OutOfCoreTable import OutOfCoreSortedTable
from peptidescan.PeptideRemapper import PeptideRemapper, UniProtIsoformAcc, RefSeqAcc, UCSCKGAcc, UniProtIsoformRefSeq, RefSeqUniProtIsoform, GencodeENSP
from operator import itemgetter

class SequenceDatabase(object):
    _orgmap = {}
    _acc = None
    def __init__(self,id,filename=None,organism=None,version=None,priority=None,decoyPrefix=None):
        self._id = id
        self._organism = organism
        self._version = version
        self._filename = filename
        self._decoyPrefix = decoyPrefix
        self._priority = priority
        self._pm = None
        assert self._name
        assert not organism or organism in self._orgmap

    def id(self):
        return self._id

    def name(self):
        return self._name

    def filename(self):
        return self._filename

    def priority(self):
        return self._priority

    def location(self):
        if self._filename:
            return os.path.split(self._filename)[1]
        return None

    def url(self):
        return None

    def organism(self):
        return self._organism

    def version(self):
        return self._version

    def dbsource(self):
        return None

    def orgmap(self):
        if not self._organism in self._orgmap:
            return None
        return self._orgmap[self._organism]

    def remapper(self):
        return bool(self._pm)

    def peptidemap(self,peptides):
        if not self._acc:
            return
        if not self._filename:
            return
        # print(self._filename, self._acc)
        self._pm = PeptideRemapper(peptides, self._filename, self._acc, preprocess=True)

    def proteins(self,peptide):
        assert(self._pm)
        dbid = self.id()
        for pracc,laa,start,end,raa,prdefline,prlen in map(itemgetter(0,2,3,4,5,6,9),self._pm.proteins(peptide)):
            yield dbid,pracc,laa[-1],raa[0],start+1,end,prlen,prdefline

    def protein_defline(self,acc):
        if not self._pm:
            return None
        return self._pm.protdb.get(acc).defline

    def protein_coverage(self,acc):
        if not self._pm:
            return None
        return self._pm.protdb.get(acc).coverage()

    def prefer(self,a,b):
        assert(self._acc)
        return self._acc.prefer(a,b)

    def decoyPrefix(self):
        return self._decoyPrefix

    def setDecoyPrefix(self,prefix):
        self._decoyPrefix = prefix

    @staticmethod
    def seqdbs(seqdir,seqdbs):
        from configparser import ConfigParser
        seqdbconfig = ConfigParser()
        seqdbconfig.read(os.path.join(seqdir,'seqdb.ini'))

        for sec in seqdbconfig.sections():
            if sec.startswith('Organism:'):
                SequenceDatabase._orgmap[sec.split(':',1)[1]] = dict(seqdbconfig.items(sec))
        sdbs = {}
        seen = set()
        for i,sec in enumerate(seqdbs):
            if sec in seen:
                continue
            seen.add(sec)
            if seqdbconfig.has_section(sec):
                kwargs = dict(seqdbconfig.items(sec))
                kwargs['priority'] = i+1
                kwargs['id'] = sec
                sdb = SequenceDatabase.new(seqdir,**kwargs)
                sdbs[sdb.id()] = sdb
        return sdbs

    @staticmethod
    def new(seqdir,**kwargs):
        filename = kwargs.get('filename')
        if filename and not os.path.exists(filename):
            filename = os.path.join(seqdir,filename)
        assert (not filename) or os.path.exists(filename)
        id = kwargs.get('id')
        assert id
        version = kwargs.get('version')
        organism = kwargs.get('organism')
        source = kwargs.get('source')
        priority = kwargs.get('priority')
        decoyPrefix = kwargs.get('decoyPrefix')
        try:
            return eval(source)(id,filename,organism,version,priority,decoyPrefix)
        except NameError:
            pass
        raise RuntimeError("Bad sequence database source %s"%source)

class RefSeq(SequenceDatabase):
    _name = 'RefSeq'
    _acc = RefSeqAcc()
    def __init__(self,*args,**kwargs):
        super(RefSeq,self).__init__(*args,**kwargs)

    def url(self):
        if not self.orgmap():
            return None
        return "ftp://ftp.ncbi.nlm.nih.gov/refseq/%(refseqsci)s/mRNA_Prot/%(refseqorg)s.*.protein.faa.gz"%self.orgmap()

    def dbsource(self):
        return "DB source RefSeq"

class UniProt(SequenceDatabase):
    _name = 'UniProt'
    _acc = UniProtIsoformAcc()
    def __init__(self,*args,**kwargs):
        super(UniProt,self).__init__(*args,**kwargs)

    def url(self):
        if not self.orgmap():
            return None
        return "https://rest.uniprot.org/uniprotkb/stream?query=keyword%3AKW-1185+AND+organism_id%3A%(taxid)s&format=fasta&force=true&compressed=true&includeIsoform=true"%self.orgmap()

    def dbsource(self):
        return "DB source UniProt"

class Gencode(SequenceDatabase):
    _name = 'Gencode'
    _acc = GencodeENSP()
    def __init__(self,*args,**kwargs):
        super(Gencode,self).__init__(*args,**kwargs)

    def url(self):
        if not self.orgmap():
            return None
        assert self.organism() in ("Human","Mouse")
        metadata = {}
        metadata['release'] = self.version().lower().strip('v')
        metadata['organism'] = self.organism().lower()
        return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_%(organism)s/release_%(release)s/gencode.v%(release)s.pc_translations.fa.gz"%metadata

    def dbsource(self):
        return "DB source Gencode"

class Decoy(SequenceDatabase):
    _name = 'Decoy'

class Unknown(SequenceDatabase):
    _name = 'Unknown'

class Contaminant(SequenceDatabase):
    _name = 'Contaminant'

class PSMCache(object):

    def __init__(self,input,*extraheaders):
        headers = [ 'SpectrumFile', 'Location', 'SpectrumID',
                    'SpectrumIDFormat', 'FileFormat', 'Scan',
                    'ExperimentalMassToCharge', 'Rank', 'Peptide',
                    'retention time', 'intensity of precursor ion',
                    'ChargeState', 'Modification', 'Protein' ]
        for eh in extraheaders:
            headers.extend(eh)
        self.t = OutOfCoreSortedTable(rows=input,headers=headers,
                                      key=lambda r: (r['SpectrumFile'],r['Scan']))

    def peptides(self):
        return map(lambda r: r['Peptide'], self.t)

    def psms(self):
        return self.t

class PSMFormater(object):

    def __init__(self,scores,params,specparams):
        self.scores = scores
        self.params = params
        self.specparams = specparams

    def write_metadata(self,md):
        md['FileFormat'] = 'MS:1000584' #mzML file
        md['Software'].append('textpsm2mzid (md5:b8bee0ddf08678274d7f5e55822fae32)')
        md['Software'].append('ProteoWizard r22286')

        print("MDBEGIN")
        for k,v in list(md.items()):
            if k != "Software":
                print(k,v)
            else:
                for v1 in md[k]:
                    print(k,v1)
        print("MDEND")

    def write_prgrp_metadata(self,grpmd):
        grpmd['Threshold'] = 'no\\ threshold'
        print("GRPMDBEGIN")
        for k,v in list(grpmd.items()):
            if k != "AnalysisParams":
                print(k,v)
            else:
                for v1 in grpmd[k]:
                    print(k,v1)
        print("GRPMDEND")

    def write_seqdb(self,sdb):
        print("SEQDBBEGIN")
        print("ID",sdb.id())
        print("Name",sdb.name())
        if sdb.organism():
            print("Organism",sdb.organism())
        if sdb.version():
            print("Release",sdb.version())
        if sdb.url():
            print("URI",sdb.url())
        if sdb.dbsource():
            print("DBSource",sdb.dbsource())
        if sdb.location():
            print("Location",sdb.location())
        if sdb.decoyPrefix():
            print("DecoyPrefix",sdb.decoyPrefix())
        print("SEQDBEND")

    def write_psm(self,psm):
        print("PSMBEGIN")
        for key in "SpectrumFile Location FileFormat SpectrumID SpectrumIDFormat Scan Rank ChargeState ExperimentalMassToCharge Peptide".split():
            if psm.get(key) not in (None,""):
                print(key,psm[key])
        for m in psm['Modification']:
            print("Modification"," ".join(map(str,m)))
        for key in self.scores:
            if psm.get(key) not in (None,""):
                value = psm[key]
                key1 = re.sub(' ','\\ ',key)
                print("Score",key1,value)
        for key in self.params:
            if psm.get(key) not in (None,""):
                key1 = re.sub(' ','\\ ',key)
                print("Param",key1,psm[key])
        for key in self.specparams:
            if psm.get(key) not in (None,""):
                key1 = re.sub(' ','\\ ',key)
                if isinstance(psm[key],str):
                    print("SpecParam",key1,psm[key])
                else:
                    units = re.sub(' ','\\ ',psm[key][1])
                    print("SpecParam",key1,psm[key][0],units)

        for pr in psm['Protein']:
            print("Protein",pr[0]+":"+pr[1]," ".join(map(str,pr[2:])))

        print("PSMEND")

    def write_prgrp(self,index,npep,proteins):
        print("PRGRPBEGIN")
        print("Name Group%d"%(index,))
        print("Param - distinct\\ peptide\\ sequences",npep)
        for pr in proteins:
            prstr = pr['id']
            print("Protein",prstr)
            if pr.get('anchor'):
                print("Param",prstr,"anchor\\ protein")
            if pr.get('sameset'):
                print("Param",prstr,"sequence\\ same-set\\ protein",pr['sameset'])
            if pr.get('subset'):
                print("Param",prstr,"sequence\\ sub-set\\ protein",pr['subset'])
            print("Param",prstr,"distinct\\ peptide\\ sequences",pr['peptides'])
            if pr.get('coverage'):
                print("Param",prstr,"sequence\\ coverage",round(pr['coverage'],2),"percent")
        print("PRGRPEND")
