#!/usr/bin/env python3

import sys, os, os.path, gzip, zlib, math, re, csv
from base64 import b64decode
from collections import defaultdict
from array import array
from find_elementtree import ET

def present(s):
    return (s != None)

specparams = [
    ('ms level', 'MS:1000511', int),
    ('MS1 spectrum', 'MS:1000579', present),
    ('MSn spectrum', 'MS:1000580', present),
    ('master scan number',None,int)
    ]

selionparams = [
    ('selected ion m/z', 'MS:1000744', str),
    ('charge state', 'MS:1000041', str),
    ('peak intensity', 'MS:1000042', str),
    ]
activationparams = [ 
    ('beam-type collision-induced dissociation','MS:1000422',present),
    ('collision energy','MS:1000045',float),
]
scanparams = [
    # convert to seconds...
    ('scan start time', 'MS:1000016', lambda s: 60.0*float(s)),
    ]

bdaparams = [
    ('32-bit float','MS:1000521', present),
    ('64-bit float','MS:1000523', present),
    ('zlib compression','MS:1000574', present),
    ('m/z array','MS:1000514', present),
    ('intensity array','MS:1000515', present),
]

def getcvparams(ele,paramlist):
    values = {}
    for chele in list(ele.findall(CVP)) + list(ele.findall(USP)):
        for p in paramlist:
            if (p[0] != None and chele.attrib['name'] == p[0]) or \
               (p[1] != None and chele.attrib.get('accession') == p[1]):
                values[chele.attrib['name']] = p[2](chele.attrib.get('value'))
    return values

ns = '{http://psi.hupo.org/ms/mzml}'
SPEC = ns+'spectrum'
CVP  = ns+'cvParam'
USP  = ns+'userParam'
PRE  = ns+'precursorList/' + ns + 'precursor'
SEL  = ns+'selectedIonList/' + ns + 'selectedIon'
ACT  = ns+'activation'
SCAN = ns+'scanList/' + ns + 'scan'
BDA  = ns+'binaryDataArrayList/' + ns + 'binaryDataArray'
BINARY = ns+'binary'

def iterprecursor(infile):
    if infile.lower().endswith('.gz'):
        h = gzip.open(infile)
    else:
        h = open(infile,'rb')
        if h.read(2) == b'\x1f\x8b':
            h = gzip.open(infile)
        else:
            h.seek(0)
    for event, ele in ET.iterparse(h):
        if ele.tag == SPEC:
            specid = ele.attrib['id']
            cvparams = getcvparams(ele,specparams)
            msn = cvparams.get('MSn spectrum',False)
            mslevel = cvparams['ms level']
            scan = ele.find(SCAN)
            cvparams = getcvparams(scan,scanparams)
            rt = cvparams['scan start time']
            if msn:
                prec = ele.find(PRE)
                selected = prec.find(SEL)
                cvparams = getcvparams(selected,selionparams)
                precursormz = cvparams['selected ion m/z']
                precursorz = cvparams.get('charge state')
                precursorit = cvparams.get('peak intensity')
                yield specid,rt,precursormz,precursorz,precursorit
    h.close()

reporters = None

import sys, os.path, configparser

def labelfilepath():
    if getattr(sys, 'frozen', False):
        prog = sys.executable
    else:
        prog = __file__
    prog = os.path.abspath(os.path.realpath(prog))
    path = os.path.split(prog)[0]
    return os.path.join(path,"labels.ini")

def get_reporter_ions(labelname):
    global reporters
    if not reporters:
        reporters = dict()
        config = configparser.ConfigParser()
        config.optionxform = str
        config.read(labelfilepath())
        alias = []
        for sec in config.sections():
            reporters[sec] = (dict(),dict())
            for k,v in config.items(sec):
                if k not in ("tolerance","resolution","type","plex","tags","extra","alias"):
                    reporters[sec][0][k] = float(v)
                elif k in ("tolerance","resolution"):
                    reporters[sec][1][k] = float(v)
                elif k in ("plex",):
                    reporters[sec][1][k] = int(v)
                elif k in ("type",):
                    reporters[sec][1][k] = v
                elif k in ("tags","extra"):
                    reporters[sec][1][k] = [ vi.strip() for vi in v.split(',') ]
                elif k == "alias":
                    alias.append((v,sec))
            assert(len(reporters[sec][0]) == (reporters[sec][1]['plex']+len(reporters[sec][1].get('extra',[]))))
            assert(len(reporters[sec][1]['tags']) == reporters[sec][1]['plex'])
        for a,sec in alias:
            reporters[a] = (dict(reporters[sec][0].items()),dict(reporters[sec][1].items()))
    return reporters[labelname]

def iterreporters(infile,labels,**kw):
    ions,ionmd = get_reporter_ions(labels)
    tolerance = kw.get('tolerance',ionmd['tolerance'])
    resolution = kw.get('resolution',ionmd['resolution'])
    if labels.startswith('MS3-'):
        mslevel = 3
    else:
        mslevel = 2
    lowmz = min(ions.values())-2*tolerance
    highmz = max(ions.values())+2*tolerance

    ionmap = defaultdict(list)
    for lab,labmz in ions.items():
        labmzlow = int(math.floor((labmz-tolerance)*1000))
        labmzhigh = int(math.ceil((labmz+tolerance)*1000))
        for i in range(labmzlow,labmzhigh+1):
            ionmap[i].append(lab)

    if infile.lower().endswith('.gz'):
        h = gzip.open(infile)
    else:
        h = open(infile,'rb')
        if h.read(2) == b'\x1f\x8b':
            h = gzip.open(infile)
        else:
            h.seek(0)
    for event, ele in ET.iterparse(h):
        if ele.tag == SPEC:
            specid = ele.attrib['id']
            cvparams = getcvparams(ele,specparams)
            msn = cvparams.get('MSn spectrum',False) and cvparams.get('ms level',2) == mslevel
            if msn:
                if mslevel == 3:
                    specidelts = {}
                    keys = []
                    for elt in specid.split():
                        k,v = elt.split('=')
                        keys.append(k)
                        specidelts[k] = v
                    specidelts['scan'] = str(cvparams['master scan number'])
                    specid = " ".join("%s=%s"%(k,specidelts[k]) for k in keys)
                mz = []; it = []
                for bda in ele.findall(BDA):
                    cvparams1 = getcvparams(bda,bdaparams)
                    assert cvparams1.get('64-bit float',False) or cvparams1.get('32-bit float',False)
                    ftype = ('f' if cvparams1.get('32-bit float',False) else 'd')
                    uncomp = ((lambda x: x) if not cvparams1.get('zlib compression',False) else zlib.decompress)
                    bintxt = bda.findtext(BINARY)
                    if bintxt:
                        values = array(ftype,uncomp(b64decode(bintxt)))
                        if sys.byteorder == 'big':
                            values.byteswap()
                        if cvparams1.get('m/z array',False):
                            mz  = values
                        if cvparams1.get('intensity array',False):
                            it = values
                assert(len(mz) == len(it))
                totalab = sum(it)
                peaks = list(sorted([ (mz,it) for (mz,it) in zip(mz,it) if lowmz <= mz <= highmz ]))
                labit = dict((k,0.0) for k in ions.keys())
                labdel = dict((k,None) for k in ions.keys())
                labqual = dict((k,None) for k in ions.keys())
                for mz,it in peaks:
                    for lab in ionmap[int(round(mz*1000))]:
                        labmz = ions[lab]
                        bucketdel = lambda d: int(math.floor(abs(d*1000)))
                        if mz >= labmz-tolerance and mz <= labmz+tolerance and (labdel[lab] == None or bucketdel(mz-labmz) < bucketdel(labdel[lab]) or (bucketdel(mz-labmz) == bucketdel(labdel[lab]) and it > labit[lab])):
                            labit[lab] = it
                            labdel[lab] = (mz-labmz)
                            labqual[lab] = resolution*(mz-labmz)/labmz
                data = {}
                for lab in ions:
                    data[lab] = (labit[lab],labdel[lab],labqual[lab])
                data['_total'] = sum(labit.values())
                if totalab > 0:
                    data['_frac'] = data['_total']/totalab
                else:
                    data['_frac'] = 0.0
                yield specid,data

    h.close()

def iterms2(infile):
    if infile.lower().endswith('.gz'):
        h = gzip.open(infile)
    else:
        h = open(infile,'rb')
        if h.read(2) == b'\x1f\x8b':
            h = gzip.open(infile)
        else:
            h.seek(0)
    for event, ele in ET.iterparse(h):
        if ele.tag == SPEC:
            specid = ele.attrib['id']
            # print(specid)
            cvparams = getcvparams(ele,specparams)
            msn = cvparams.get('MSn spectrum',False)
            msl = int(cvparams.get('ms level',0))
            if msn and msl == 2:
                mz = []; it = []
                for bda in ele.findall(BDA):
                    cvparams1 = getcvparams(bda,bdaparams)
                    assert cvparams1.get('64-bit float',False) or cvparams1.get('32-bit float',False)
                    ftype = ('f' if cvparams1.get('32-bit float',False) else 'd')
                    uncomp = ((lambda x: x) if not cvparams1.get('zlib compression',False) else zlib.decompress)
                    bintxt = bda.findtext(BINARY)
                    if bintxt:
                        values = array(ftype,uncomp(b64decode(bintxt)))
                        if sys.byteorder == 'big':
                            values.byteswap()
                        if cvparams1.get('m/z array',False):
                            mz  = values
                        if cvparams1.get('intensity array',False):
                            it = values
                assert(len(mz) == len(it))
                extraparams = dict()
                scan = ele.find(SCAN)
                cvparams = getcvparams(scan,scanparams)
                rt = cvparams['scan start time']
                prec = ele.find(PRE)
                precid = prec.attrib['spectrumRef']
                selected = prec.find(SEL)
                cvparams = getcvparams(selected,selionparams)
                precursormz = cvparams['selected ion m/z']
                try:
                    precursorz = int(cvparams.get('charge state'))
                except (TypeError,ValueError):
                    precursorz = ""
                precursorit = cvparams.get('peak intensity')
                activation = prec.find(ACT)
                cvparams = getcvparams(activation,activationparams)
                if cvparams.get('beam-type collision-induced dissociation',False):
                    extraparams['hcdenergy'] = cvparams['collision energy']
                yield dict(peaks=sorted(zip(mz,it)), 
                           precursor=dict(mz=precursormz, z=precursorz, it=precursorit, id=precid),
                           metadata=dict(rt=rt, id=specid, **extraparams))

    h.close()

def specid(idstr):
    scankvpairs = re.split(r' ([a-z]+)='," "+idstr)
    scandata = dict()
    for i in range(1,len(scankvpairs),2):
        scandata[scankvpairs[i]] = scankvpairs[i+1]
    return scandata

def write_mgf(infile,out=None):
    if not out:
        out = sys.stdout
    
    for s in iterms2(infile):
        peaks = s['peaks']
        if len(peaks) == 0:
            continue
        precursor = s['precursor']
        metadata = s['metadata']
        metadata['scan'] = specid(metadata['id'])['scan']
        metadata['precscan'] = specid(precursor['id'])['scan']
        print("BEGIN IONS",file=out)
        print("TITLE=Scan:%(scan)s RT:%(rt)s HCD:%(hcdenergy)s PrecursorScan:%(precscan)s"%metadata,file=out)
        print("CHARGE=%(z)s"%precursor,file=out)
        print("PEPMASS=%(mz)s"%precursor,file=out)
        for p in sorted(peaks):
            print("%f %f"%p,file=out)
        print("END IONS\n",file=out)

def xmlindent(current, parent=None, index=-1, depth=0):
    for i, node in enumerate(current):
        xmlindent(node, current, i, depth + 1)
    if parent is not None:
        if index == 0:
            parent.text = '\n' + ('  ' * depth)
        else:
            parent[index - 1].tail = '\n' + ('  ' * depth)
        if index == len(parent) - 1:
            current.tail = '\n' + ('  ' * (depth - 1))

def write_phosphors(mzmlfile,psmfile,modfile,out=None,tolerance=0.05,activation='HCD'):
    if not out:
        out = sys.stdout

    moddata = dict()
    headers = "ch,aas,matchdelta,name,delta,loss".split(',')
    # 1       S       79.966  Phospho 79.966331       97.976896
    # 1       T       79.966  Phospho 79.966331       97.976896
    for line in open(modfile):
        row = dict(zip(headers,line.split()))
        row['delta'] = float(row['delta'])
        row['loss'] = float(row['loss'])
        matchdelta = row['matchdelta']
        if float(matchdelta)>=0:
            matchdelta = "+"+matchdelta             
        else:
            matchdelta = "-"+matchdelta             
        row['matchdelta'] = matchdelta 
        if matchdelta in moddata:
            moddata[matchdelta]['aas'] += row['aas']
        else:
            moddata[matchdelta] = row

    psmdata = defaultdict(list)
    for row in csv.DictReader(open(psmfile),dialect='excel-tab'):
        if 'Peptide' in row:
            peptide = row['Peptide']
        else:
            peptide = row['PeptideSequence']
        splpep = re.split(r'([A-Z])',peptide)
        pepseq = "".join(splpep[1::2])
        mods = []
        for i,mi in enumerate(splpep[0::2]):
            if mi in moddata:
                if i == 0 and '[' in moddata[mi]['aas']:
                    mods.append(moddata[mi]['ch'])
                elif pepseq[i-1] in moddata[mi]['aas']:
                    mods.append(moddata[mi]['ch'])
                else:
                    raise LookupError(pepseq[i-1],mi)
            elif mi == "":
                mods.append('0')
            else:
                raise LookupError(mi)
        mods = mods[0]+"."+"".join(mods[1:])+".0"
        scannum = int(row['ScanNum'])
        psmdata[scannum].append(dict(scannum=scannum,pepseq=pepseq,mods=mods,peptide=peptide))

    root = ET.Element("phosphoRSInput")
    phosphorsxml = ET.ElementTree(root)
    ET.SubElement(root,"MassTolerance",Value=str(tolerance))
    ET.SubElement(root,"Phosphorylation",Symbol='1')
    spectra = ET.SubElement(root,"Spectra")
    
    for s in iterms2(mzmlfile):
        peaks = s['peaks']
        if len(peaks) == 0:
            continue
        precursor = s['precursor']
        metadata = s['metadata']
        metadata['scan'] = specid(metadata['id'])['scan']
        scan = int(metadata['scan'])
        if scan not in psmdata:
            continue
        spectrum = ET.SubElement(spectra,"Spectrum",ID=metadata['scan'],
                                                    PrecursorCharge=str(precursor['z']),
                                                    ActivationTypes=activation)
        peakselt = ET.SubElement(spectrum,"Peaks")
        peakstext = []
        for p in sorted(peaks):
            peakstext.append(":".join(map(str,p)))
        peakselt.text = ",".join(peakstext)
        peptides = ET.SubElement(spectrum,"IdentifiedPhosphorPeptides")
        for psm in psmdata[scan]:
            ET.SubElement(peptides,"Peptide",ID="_".join([metadata['scan'],psm['peptide']]),
                                             Sequence=psm['pepseq'],
                                             ModificationInfo=psm['mods'])

    modinfos = ET.SubElement(root,"ModificationInfos")
    for mod in moddata.values():
        if mod.get('loss',0.0) > 0:
            val = "%(ch)s:%(name)s:%(name)s:%(delta)s:PhosphoLoss:%(loss)s:%(aas)s"%mod
        else:
            val = "%(ch)s:%(name)s:%(name)s:%(delta)s:null:0:%(aas)s"%mod
        ET.SubElement(modinfos,"ModificationInfo",Symbol=mod['ch'],Value=val)

    xmlindent(root)
    phosphorsxml.write(out,encoding="unicode")

def write_reporters(infile,labelname,*args,**kw):
    out = sys.stdout
    labels,labelmd = get_reporter_ions(labelname)
    for sid,data in iterreporters(infile,labelname,*args,**kw):
        line = []
        scan = specid(sid)['scan']
        line.append(scan)
        if labelname.startswith('MS3-'):
            labelname = labelname.split('-',1)[1]
        if 'extra' not in labelmd:
            line.append(labelname)
        else:
            line.append("%s,+%d"%(labelname,len(labelmd['extra'])))
        for tag in labelmd['tags']:
            line.append(float("%.5e"%data[tag][0]))
        for tag in labelmd.get('extra',[]):
            line.append(float("%.5e"%data[tag][0]))
        for tag in labelmd['tags']:
            line.append("%.2f"%data[tag][2] if data[tag][2] != None else '?')
        for tag in labelmd.get('extra',[]):
            line.append("%.2f"%data[tag][2] if data[tag][2] != None else '?')
        line.append(float("%.5e"%data['_frac']) if max(data[t][0] for t in (labelmd['tags']+labelmd.get('extra',[]))) > 0 else "?")
        # line.append(float("%.5e"%data['_total']))
        print("\t".join(map(str,line)),file=out)

def add_spec_metadata(infile,psmfile):
    out = sys.stdout
    specmd = defaultdict(dict)
    for s in iterms2(infile):
        metadata = s['metadata']
        precursor = s['precursor']
        scan = int(specid(metadata['id'])['scan'])
        specmd[scan]['PrecursorScanNum'] = specid(precursor['id'])['scan']
        specmd[scan]['OriginalPrecursorMz'] = precursor['mz']
        specmd[scan]['OriginalCharge'] = precursor['z']
        if 'hcdenergy' in metadata:
            specmd[scan]['HCDEnergy'] = metadata['hcdenergy']
    reader = csv.DictReader(open(psmfile),dialect='excel-tab')
    writer = None
    for row in reader:
        if writer == None:
             fieldnames = reader.fieldnames
             # if 'PrecursorPurity' in fieldnames:
             #     fieldnames.remove('PrecursorPurity')
             # if 'FractionDecomposition' in fieldnames:
             #     fieldnames.remove('FractionDecomposition')
             writer = csv.DictWriter(out,fieldnames,extrasaction='ignore',dialect='excel-tab')
             writer.writeheader()
        scan = int(row['ScanNum'])
        if scan in specmd:
            row.update(specmd[scan])
        writer.writerow(row)

if __name__ == '__main__':

    import sys, os

    # for t in iterprecursor(sys.argv[1]):
    #     print('\t'.join(map(str,t)))

    cmd = sys.argv[1]
    args = sys.argv[2:]

    if cmd == "write_mgf":
        write_mgf(*args[:2])

    elif cmd == "write_reporters":
        write_reporters(*args[:2])
            
    elif cmd == "add_spec_metadata":
        add_spec_metadata(*args[:2])

    elif cmd == "write_phosphors":
        write_phosphors(*args[:3],tolerance=args[3],activation=args[4],)
