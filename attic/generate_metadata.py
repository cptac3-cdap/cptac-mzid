#!/usr/bin/env python27
import re
from optparse import OptionParser
parser = OptionParser()

parser.add_option("--sheet",type=str,dest="sheet",default=None)
parser.add_option("--fractions",type=str,dest="fractions",default=None)
parser.add_option("--replicates",type=str,dest="replicates",default='1')
parser.add_option("--numberfmt",type=str,dest="numberfmt",default="%s")
parser.add_option("--pcc",type=str,dest="pcc",default=None)
parser.add_option("--pi",type=str,dest="pi",default=None)
parser.add_option("--asprot",type=str,dest="asprot",default=None)
parser.add_option("--lcprot",type=str,dest="lcprot",default=None)
parser.add_option("--msprot",type=str,dest="msprot",default=None)
parser.add_option("--output",type=str,dest="output",default=None)
opts,args = parser.parse_args()

assert opts.pcc
assert opts.pi
assert opts.asprot
assert opts.lcprot
assert opts.msprot
assert opts.sheet
assert opts.fractions
assert opts.replicates

fractions = []
for fs in opts.fractions.split(','):
    fs = fs.strip()
    m = re.search('^(\d+)-(\d+)$',fs)
    if m:
        f = map(lambda n: (n,opts.numberfmt%n),range(int(m.group(1)),int(m.group(2))+1))
    else:
        try:
            f = [ (fs,fs) ]
            f = [ (int(fs),opts.numberfmt%int(fs)) ]
        except ValueError:
            pass
    fractions.extend(f)

replicates = []
for fs in opts.replicates.split(','):
    fs = fs.strip()
    m = re.search('^(\d+)-(\d+)$',fs)
    if m:
        f = map(lambda n: (n,opts.numberfmt%n),range(int(m.group(1)),int(m.group(2))+1))
    else:
        try:
            f = [ (fs,fs) ]
            f = [ (int(fs),opts.numberfmt%int(fs)) ]
        except ValueError:
            pass
    replicates.extend(f)

folders = []
from dataset import XLSXFileTable
rows = XLSXFileTable(args[0],sheet=opts.sheet)
for r in rows:
    if not r.get("File Name"):
        continue
    if r.get("PCC","") != opts.pcc:
        continue
    # print r
    d = dict(Folder=r["Folder Name"])
    d["114-Biospecimen"] = r["114-Biospecimen"]
    d["115-Biospecimen"] = r["115-Biospecimen"]
    d["116-Biospecimen"] = r["116-Biospecimen"]
    d["117-Biospecimen"] = r["117-Biospecimen"]
    # JHU filename inconsistencies...
    # d["FilenameTemplate"] = re.sub(r'0?1\.raw$','%(fraction)s.raw',r["File Name"])
    # BI/PNNL filenames...
    # d["FilenameTemplate"] = r["File Name"].replace('_f01.raw','_f%s.raw')
    d["FilenameTemplate"] = r["File Name"].replace('_RUN1_','_RUN%(replicate)s_')
    date = str(r["Date"])
    d["Date"] = "/".join(map(str,[date[4:6],date[6:8],date[:4]]))
    folders.append(d)

headers = filter(None,map(str.strip,"""
PCC
Lab
114-Biospecimen
114-Aliquot
115-Biospecimen
115-Aliquot
116-Biospecimen
116-Aliquot
117-Biospecimen
117-Aliquot
Analytical Sample Protocol
Chromatography Protocol
Mass Spectrometry Protocol
Replicate
Fraction
Date
Instrument
Operator
Folder
Filename
""".splitlines()))
output = []
for fo in folders:
    for rep in replicates:
        for fr in fractions:
            d = dict(PCC=opts.pcc,Lab=opts.pi)
            d["Analytical Sample Protocol"] = opts.asprot
            d["Chromatography Protocol"] = opts.lcprot
            d["Mass Spectrometry Protocol"] = opts.msprot
            d["Fraction"] = fr[0]
            d["Replicate"] = rep[0]
            d["Date"] = fo["Date"]
            d["114-Biospecimen"] = fo["114-Biospecimen"]
            d["115-Biospecimen"] = fo["115-Biospecimen"]
            d["116-Biospecimen"] = fo["116-Biospecimen"]
            d["117-Biospecimen"] = fo["117-Biospecimen"]
            d["Folder"] = fo["Folder"]
            # print fo["FilenameTemplate"],fr
            d["Filename"] = fo["FilenameTemplate"]%dict(replicate=rep[1],fraction=fr[1])
            # JHU special
            d["Filename"] = re.sub(r'[fF]POOL\.raw','POOL.raw',d["Filename"])
            output.append(d)
            # print d

out = XLSXFileTable(opts.output,sheet="Metadata",headers=headers)
out.from_rows(output)
