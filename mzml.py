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
    for chele in ele.findall(CVP):
        for p in paramlist:
            if chele.attrib['name'] == p[0] or \
               chele.attrib['accession'] == p[1]:
                values[chele.attrib['name']] = p[2](chele.attrib.get('value'))
    return values

ns = '{http://psi.hupo.org/ms/mzml}'
SPEC = ns+'spectrum'
CVP  = ns+'cvParam'
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
                if k not in ("tolerance","resolution","type","plex","tags","alias"):
                    reporters[sec][0][k] = float(v)
                elif k in ("tolerance","resolution"):
                    reporters[sec][1][k] = float(v)
                elif k in ("plex",):
                    reporters[sec][1][k] = int(v)
                elif k in ("type",):
                    reporters[sec][1][k] = v
                elif k in ("tags",):
                    reporters[sec][1][k] = [ vi.strip() for vi in v.split(',') ]
                elif k == "alias":
                    alias.append((v,sec))
            assert(len(reporters[sec][0]) == reporters[sec][1]['plex'])
            assert(len(reporters[sec][1]['tags']) == reporters[sec][1]['plex'])
        for a,sec in alias:
            reporters[a] = (dict(reporters[sec][0].items()),dict(reporters[sec][1].items()))
    return reporters[labelname]

def iterreporters(infile,labels,**kw):
    ions,ionmd = get_reporter_ions(labels)
    tolerance = kw.get('tolerance',ionmd['tolerance'])
    resolution = kw.get('resolution',ionmd['resolution'])

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
    for event, ele in ET.iterparse(h):
        if ele.tag == SPEC:
            specid = ele.attrib['id']
            cvparams = getcvparams(ele,specparams)
            msn = cvparams.get('MSn spectrum',False)
            if msn:
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
                        if mz >= labmz-tolerance and mz <= labmz+tolerance and it > labit[lab]:
                            labit[lab] = it
                            labdel[lab] = (mz-labmz)
                            labqual[lab] = resolution*(mz-labmz)/labmz
                data = {}
                for lab in ions:
                    data[lab] = (labit[lab],labdel[lab],labqual[lab])
                data['_total'] = sum(labit.values())
                data['_frac'] = data['_total']/totalab
                yield specid,data

    h.close()

def iterms2(infile):
    if infile.lower().endswith('.gz'):
        h = gzip.open(infile)
    else:
        h = open(infile,'rb')
    for event, ele in ET.iterparse(h):
        if ele.tag == SPEC:
            specid = ele.attrib['id']
            cvparams = getcvparams(ele,specparams)
            msn = cvparams.get('MSn spectrum',False)
            if msn:
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
                precursorz = int(cvparams.get('charge state'))
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

def write_reporters(infile,labelname,*args,**kw):
    out = sys.stdout
    labels,labelmd = get_reporter_ions(labelname)
    for sid,data in iterreporters(infile,labelname,*args,**kw):
        line = []
        scan = specid(sid)['scan']
        line.append(scan)
        line.append(labelname)
        for tag in labelmd['tags']:
            line.append(float("%.5e"%data[tag][0]))
        for tag in labelmd['tags']:
            line.append("%.2f"%data[tag][2] if data[tag][2] != None else '?')
        line.append(float("%.5e"%data['_frac']) if max(data[t][0] for t in labelmd['tags']) > 0 else "?")
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
