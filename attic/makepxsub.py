#!/usr/bin/env python27

import os, os.path, sys, shutil, re, glob
from ConfigParser import RawConfigParser as ConfigParser
from collections import defaultdict

if len(sys.argv) < 2:
    sys.exit(1)
clean = False
if sys.argv[1] == '--clean':
    clean = True
    sys.argv.pop(1)
if len(sys.argv) < 5:
    sys.exit(1)
d = sys.argv[1].rstrip(os.sep)
study = sys.argv[2]
subproteome = sys.argv[3]
site = sys.argv[4]

iniFile = sys.argv[0].rsplit('.',1)[0]+'.ini'
config = ConfigParser()
config.read([iniFile])

from dataset import XLSXFileTable
prfile = glob.glob("%s/%s_Protocols.xlsx"%(d,d))[0]
mdfile = glob.glob("%s/%s_Metadata.xlsx"%(d,d))[0]
asp = XLSXFileTable(prfile,sheet='Analytical Sample Protocol')
chp = XLSXFileTable(prfile,sheet='Chromatography Protocol')
msp = XLSXFileTable(prfile,sheet='Mass Spectrometry Protocol')
metadatafile = XLSXFileTable(mdfile)

protocols = {}
for p,t in zip([asp,chp,msp],["ASP","CHP","MSP"]):
    headers = map(str.strip,p.headers())
    for r in p:
        field = str(r.get('Name')).strip()
        for h in headers[1:]:
            h = h.strip()
            val = r.get(h)
            if isinstance(val,basestring):
                val = val.strip()
            if not val:
                continue
            if (t,h) not in protocols:
                protocols[(t,h)] = dict()
            protocols[(t,h)][field] = str(val).strip()

metadata = {}
for r in metadatafile:
    f = r.get('LC/MS Acquisition Filename',r.get('Filename'))
    if not f:
        continue
    if f.rsplit('.')[-1].lower() in ('raw',):
        f = f.rsplit('.',1)[0]
    metadata[f] = dict(r.items())
    aspname = r['Analytical Sample Protocol']
    for k,v in protocols[("ASP",aspname)].items():
        if isinstance(v,basestring):
            v = v.strip()
        if v:
            metadata[f]["ASP:"+k] = v
    chpname = r['Chromatography Protocol']
    for k,v in protocols[("CHP",chpname)].items():
        if isinstance(v,basestring):
            v = v.strip()
        if v:
            metadata[f]["CHP:"+k] = v
    mspname = r['Mass Spectrometry Protocol']
    for k,v in protocols[("MSP",mspname)].items():
        if isinstance(v,basestring):
            v = v.strip()
        if v:
            metadata[f]["MSP:"+k] = v
    if 'iTRAQ' in metadata[f]["ASP:Type"]:
        bs = []
        for reporter in (113,114,115,116,117,118,119,121):
            if r.get('%d-Biospecimen'%reporter):
                if r.get('%d-Aliquot'%reporter):
                    bs.append("%d:%s-%s"%(reporter,r['%d-Biospecimen'%reporter],r['%d-Aliquot'%reporter]))
                else:
                    bs.append("%d:%s"%(reporter,r['%d-Biospecimen'%reporter]))
        metadata[f]["Biospecimen"] = ','.join(bs)
    else:
        if r.get('Aliquot'):
            metadata[f]["Biospecimen"] = "%s-%s"%(r['Biospecimen'],r['Aliquot'])
        else:
            metadata[f]["Biospecimen"] = "%s"%(r['Biospecimen'],)
    # JHU's POOL is not a fraction, per se...
    if metadata[f]['Fraction'] == 'POOL':
        metadata[f]['Fraction'] = metadata[f]['Fraction']
    elif metadata[f]['Fraction'] == 'N/A':
        metadata[f]['Fraction'] = "replicate %s"%(metadata[f]['Replicate'],)
    else:
        metadata[f]['Fraction'] = "fraction %s"%(metadata[f]['Fraction'],)

def getdict(config,section):
    return dict(config.items(section))

labhead = getdict(config,'Labhead:'+site)
project = getdict(config,'Project:'+study)
project['specieslist'] = re.findall(r'\[[^]]+\]',project['species'].strip())
project['tissuelist'] = re.findall(r'\[[^]]+\]',project['tissue'].strip())
project['diseaselist'] = re.findall(r'\[[^]]+\]',project['disease'].strip())
project['celltypelist'] = re.findall(r'\[[^]]+\]',project['celltype'].strip())
project['publications'] = filter(None,map(str.strip,project['publications'].split(',')))
project['links'] = filter(None,map(str.strip,project['links'].split(',')))

ontology = {}
for s in config.sections():
    if s.startswith('Ontology:'):
        key = s.split(':',1)[1]
        ont = config.get(s,"acc").split(':')[0]
        acc = config.get(s,"acc")
        name = config.get(s,"name")
        ontology[key] = (ont,acc,name)

# print >>wh, protocols

instruments = set()
quants = set()
protvals = defaultdict(lambda: defaultdict(dict))
protdesc = defaultdict(list)
for md in metadata.values():

    inst = md['MSP:Instrument']
    if inst in ontology:
        instruments.add(ontology[inst])

    lab = md['ASP:Labelling']
    if lab in ontology:
        quants.add(ontology[lab])

    for abbr,protype in [ ("ASP","Analytical Sample Protocol"),
                          ("CHP","Chromatography Protocol"),
                          ("MSP","Mass Spectrometry Protocol") ]:
        proname = md[protype]
        prokey = (abbr,proname)
        if proname not in protvals[abbr]:
            desc = []
            for k,v in sorted(protocols[prokey].items()):
                if k in ('Type','Updated','Document',):
                    continue
                if v in ('N/A','',None):
                    continue
                protvals[abbr][proname][k] = v.strip()
                desc.append("%s: %s"%(k,v))
            desc = protype + " - " + ", ".join(desc)
            if desc not in protdesc[abbr]:
                protdesc[abbr].append(desc)

protocolstr = []
for abbr in ["ASP","CHP","MSP"]:
    if len(protdesc[abbr]) == 1:
        protocolstr.append(protdesc[abbr][0])
        continue
    for i,di in enumerate(protdesc[abbr]):
        d1 = di.replace(" Protocol - "," Protocol %d - "%(i+1,))
        protocolstr.append(d1)
# print protocolstr
protocolstr = '; '.join(protocolstr)

project['instrument'] = []
for i in instruments:
    project['instrument'].append("[%s, %s, %s, ]"%i)
project['quantitation'] = []
for q in quants:
    project['quantitation'].append("[%s, %s, %s, ]"%q)
project['protocol'] = protocolstr

templatestr = """
MTD     submitter_name  Nathan Edwards
MTD     submitter_email nje5@georgetown.edu
MTD     submitter_affiliation   CPTAC Data Coordinating Center & Georgetown University
MTD     submitter_pride_login   nje5@georgetown.edu
MTD     lab_head_name   $labhead.name.strip()
MTD     lab_head_email  $labhead.email.strip()
MTD     lab_head_affiliation    $labhead.affiliation.strip()
MTD     project_title   $project.title.strip(): $subproteome, $site
MTD     project_description     $project.description.strip()
#if $tag
MTD     project_tag     $tag
#end if
MTD     keywords        CPTAC CPTAC-$study CPTAC-$site
MTD     sample_processing_protocol      $project.protocol. See associated protocols and methods documents at https://cptac-data-portal.georgetown.edu.
MTD     data_processing_protocol        Processed using the CPTAC Common Data Analysis Pipeline at NIST. Peptide identifications from MSGF+. Filtered at 1% FDR. See associated protocols and methods documents at https://cptac-data-portal.georgetown.edu.
#for l in $project.links
MTD     other_omics_link        $l
#end for
MTD     experiment_type [PRIDE, PRIDE:0000429, Shotgun proteomics, ]
MTD     submission_type COMPLETE
#for id in $project.publications
MTD     pubmed_id       $id
#end for
#for s in $project.specieslist
MTD     species $s
#end for
#for t in $project.tissuelist
MTD     tissue  $t
#end for
#for c in $project.celltypelist
MTD     cell_type       $c
#end for
#for d in $project.diseaselist
MTD     disease $d
#end for
#for i in $project.instrument
MTD     instrument      $i
#end for
#for q in $project.quantitation
MTD     quantification  $q
#end for
"""
import warnings
warnings.filterwarnings('ignore',category=UserWarning,module='Cheetah\.Compile')
from Cheetah.Template import Template

preamble = Template(templatestr.lstrip())
preamble.labhead = labhead
preamble.project = project
preamble.subproteome = subproteome
preamble.study = study
preamble.site = site
# preamble.tag = "Clinical Proteomic Tumor Analysis Consortium"
preamble.tag = "CPTAC Consortium"
# preamble.tag = None

wh = open(d+'.px','w')
print >>wh, str(preamble)

def findall(dir,ext=None,notext=None,regex=None,clean=False):
    if regex:
        regex = re.compile(regex)
    retval = []
    sfs = []
    for root, dirs, files in os.walk(dir, followlinks=True):
        dirs.sort()
        files.sort()
        if not ext:
            for d in dirs:
                if d.endswith('_mzIdentML') and d+'.cksum' not in files:
                    print >>sys.stderr, "Checksum file %s is missing"%(os.path.join(root,d+'.cksum'),)
                if d.endswith('_mzML') and d+'.cksum' not in files:
                    print >>sys.stderr, "Checksum file %s is missing"%(os.path.join(root,d+'.cksum'),)
                if d.endswith('_raw') and d+'.cksum' not in files:
                    print >>sys.stderr, "Checksum file %s is missing"%(os.path.join(root,d+'.cksum'),)
        nfiles = 0
        for name in files:
            if not ext and name.startswith('~'):
                if clean:
                    os.unlink(os.path.join(root,name))
                    print >>sys.stderr, "File %s removed"%(os.path.join(root,name),)
                else:
                    print >>sys.stderr, "File %s ignored"%(os.path.join(root,name),)
                continue
            if not ext and name.endswith('.bak'):
                if clean:
                    os.unlink(os.path.join(root,name))
                    print >>sys.stderr, "File %s removed"%(os.path.join(root,name),)
                else:
                    print >>sys.stderr, "File %s ignored"%(os.path.join(root,name),)
                continue
            if not ext and name.endswith('.log'):
                if clean:
                    os.unlink(os.path.join(root,name))
                    print >>sys.stderr, "File %s removed"%(os.path.join(root,name),)
                else:
                    print >>sys.stderr, "File %s ignored"%(os.path.join(root,name),)
                continue
            if not ext and name.endswith('_tsv.cksum'):
                if clean:
                    os.unlink(os.path.join(root,name))
                    print >>sys.stderr, "File %s removed"%(os.path.join(root,name),)
                else:
                    print >>sys.stderr, "File %s ignored"%(os.path.join(root,name),)
                continue
            if not ext and name.endswith('.psm'):
                if clean:
                    os.unlink(os.path.join(root,name))
                    print >>sys.stderr, "File %s removed"%(os.path.join(root,name),)
                else:
                    print >>sys.stderr, "File %s ignored"%(os.path.join(root,name),)
                continue
            if (ext == None or name.endswith(ext)) and \
               (regex == None or regex.search(os.path.join(root,name))) and \
               (notext == None or all(map(lambda ne: not name.endswith(ne), notext))):
                retval.append(os.path.join(root,name))
                nfiles += 1
        if not ext and nfiles == 0 and clean:
            print >>sys.stderr, "Empty directory %s removed"%(root,)
            os.rmdir(root)
    return retval

raw = findall(d,'.raw')
mzml = findall(d,'.mzML.gz')
# mzid = findall(d,'.mzid.gz',regex=r'_PSM\.CAP\.r1_mzIdentML/')
mzid = findall(d,'.mzid.gz')
assert len(raw) == len(mzml)
# assert len(raw) == len(mzid)
# nraw = len(raw)
rest = findall(d,notext=['.raw','.mzML.gz','.mzid.gz'],clean=clean)
# rest = sorted(set(rest) - set(raw) - set(mzml) - set(mzid))

def getbase(fn):
    if fn.endswith('.gz'):
        return os.path.split(fn)[1].rsplit('.',2)[0]
    return os.path.split(fn)[1].rsplit('.',1)[0]

bybase = defaultdict(dict)
file2ind = dict()
for fn in raw:
    bybase[getbase(fn)]['raw'] = fn
for fn in mzml:
    bybase[getbase(fn)]['mzml'] = fn

# def normpath(f):
#     return os.path.split(os.path.abspath(f))[1]
#
def normpath(f):
    return os.path.abspath(f)

print >>wh, '\t'.join("FMH      file_id file_type       file_path       file_mapping    ".split() + [])
ind = 0
for resultfile in sorted(mzid):
    base = getbase(resultfile)
    rawfile = bybase[base]['raw']
    peakfile = bybase[base]['mzml']
    needraw = False
    if rawfile not in file2ind:
        needraw = True
        file2ind[rawfile] = ind
        ind += 1
    needpeak = False
    if peakfile not in file2ind:
        needpeak = True
        file2ind[peakfile] = ind
        ind += 1
    file2ind[resultfile] = ind
    ind += 1
    if needraw:
        print >>wh, '\t'.join(map(str,["FME",file2ind[rawfile],"RAW",normpath(rawfile),""]))
    if needpeak:
        print >>wh, '\t'.join(map(str,["FME",file2ind[peakfile],"PEAK",normpath(peakfile),""]))
    print >>wh, '\t'.join(map(str,["FME",file2ind[resultfile],"RESULT",normpath(resultfile),"%d,%d"%(file2ind[rawfile],file2ind[peakfile])]))

for otheri in sorted(rest):
    print >>wh, '\t'.join(map(str,["FME",ind,"OTHER",normpath(otheri),""]))
    ind += 1
print >>wh, ""

metadata_seen = set()

print >>wh, '\t'.join("SMH      file_id species tissue  cell_type       disease modification    instrument      quantification  experimental_factor".split())
for resultfile in sorted(mzid):
    resibase = getbase(resultfile)
    md = metadata[resibase]
    metadata_seen.add(resibase)
    labelling = ""
    if md.get('ASP:Labelling') in ontology:
        labelling = "[%s, %s, %s, ]"%ontology[md['ASP:Labelling']]
    print >>wh, '\t'.join(map(str,["SME",
                             file2ind[resultfile],
                             project['species'],
                             project['tissue'],
                             project['celltype'],
                             project['disease'],
                             project['modification'],
                             "[%s, %s, %s, ]"%ontology[md['MSP:Instrument']],
                             labelling,
                             project['factor']%md,
                             ]))

wh.close()

for rb in sorted(metadata):
    if rb not in metadata_seen:
        print >>sys.stderr, "Metadata file %s missing"%rb
