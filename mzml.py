#!/usr/bin/env python3

import sys, os, os.path, gzip
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
scanparams = [
    # convert to seconds...
    ('scan start time', 'MS:1000016', lambda s: 60.0*float(s)),
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
SCAN = ns+'scanList/' + ns + 'scan'

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

if __name__ == '__main__':

    import sys, os

    for t in iterprecursor(sys.argv[1]):
        print('\t'.join(map(str,t)))
