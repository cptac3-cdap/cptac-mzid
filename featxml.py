#!/usr/bin/env python3

import sys, os, os.path, gzip, zlib, math, re, csv
from base64 import b64decode
from collections import defaultdict
from array import array
from find_elementtree import ET

def iterfeat(infile):
    context = []
    for event, ele in ET.iterparse(open(infile),events=('start','end')):
        if event == "start":
            context.append(ele.tag)
            # print context
            continue
        if context == [ 'featureMap', 'featureList', 'feature' ]:
            # element is a base level feature.
            f = {}
            for posele in ele.findall('position'):
                if posele.attrib['dim'] == "0":
                    f['rt'] = posele.text
                elif posele.attrib['dim'] == "1":
                    f['mz'] = posele.text
            # f['area'] = ele.find('intensity').text
            f['z'] = int(ele.find('charge').text)
            for upele in ele.findall('UserParam'):
                if upele.attrib['name'] == "model_area":
                    f['area'] = upele.attrib['value']
                elif upele.attrib['name'] == "model_height":
                    f['intensity'] = upele.attrib['value']  
                elif upele.attrib['name'] == "model_status":
                    f['status'] = upele.attrib['value']
                elif upele.attrib['name'] == "model_FWHM":
                    f['fwhm'] =  upele.attrib['value']  
            for siele in ele.findall('PeptideIdentification'):
                f1 = dict(f.items())
                specref = siele.attrib['spectrum_reference']
                f1['scan'] = int(specid(specref)['scan'])
                phele = siele.find('PeptideHit')
                f1['pepseq'] = phele.attrib['sequence']
                for upele in phele.findall('UserParam'):
                    if upele.attrib['name'] == "calcMZ":
                        f1['calcmz'] = upele.attrib['value']
                    elif upele.attrib['name'] == "MS:1002054":
                        f1['qvalue'] = float(upele.attrib['value'])
                yield f1

            ele.clear()
        elif len(context) == 2 and context != ['featureMap', 'featureList']:
            ele.clear()
        context.pop(-1)        

def specid(idstr):
    scankvpairs = re.split(r' ([a-z]+)='," "+idstr)
    scandata = dict()
    for i in range(1,len(scankvpairs),2):
        scandata[scankvpairs[i]] = scankvpairs[i+1]
    return scandata

def write_features(infile,qvalue=0.01,out=None):
    if not out:
         out = sys.stdout
    else:
         out = open(out,'w')
    keys = ['scan','rt','mz','z','pepseq','calcmz','area','intensity','fwhm','status']
    headers = ['Scan','RT','MZ','Charge','Peptide','CalcMZ','Area','Intensity','FWHM','Status']
    print >>out, "\t".join(headers)
    for f in iterfeat(infile):
        if f['qvalue'] > qvalue:
            continue
        print >>out, "\t".join(map(lambda k: str(f.get(k,"")),keys))

    if out != sys.stdout:
        out.close()

if __name__ == '__main__':

    import sys, os

    args = sys.argv[1:]
    if len(args) > 1:
        args[1] = float(args[1])
    write_features(*args)
