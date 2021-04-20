#!/usr/bin/env python2
import hashlib, sys, re, os, os.path, gzip, glob, encodings

verbose = False
if sys.argv[1] == '-v':
    verbose = True
    sys.argv.pop(1)

def checksha1(f,sha1,st=None,ed=None):
    if ed == None:
        ed = 99999999999999999999
    if st == None:
        st = 0
    if not os.path.exists(f):
        return False
    if f.endswith('.gz'):
        h = gzip.open(f,'rb')
    else:
        h = open(f,'rb')
    h.seek(st)
    totalread = st
    fsha1=hashlib.sha1()
    while True:
        readlen = min(32*1024,ed-totalread)
        if readlen <= 0:
            break
        buf = h.read(readlen)
        if not buf:
            break
        totalread += len(buf)
        fsha1.update(buf)
    h.close()
    print("",sha1,fsha1.hexdigest(),(sha1 == fsha1.hexdigest()))
    return (sha1 == fsha1.hexdigest())

def getlastblock(f,bs):
    f.seek(0)
    buf0 = ""
    buf1 = ""
    buf2 = ""
    totalsize = 0
    while True:
        buf0 = buf1
        buf1 = buf2
        buf2 = f.read(32*1024)
        if not buf2:
            break
        totalsize += len(buf2)
    if len(buf1) >= bs:
        return buf1[-bs:],totalsize
    buf0 += buf1
    return buf0[-bs:],totalsize

allargs = []
for a in sys.argv[1:]:
    allargs.extend(glob.glob(a))

for f in allargs:
    if not f.lower().endswith('.mzml') and not f.lower().endswith('.mzml.gz'):
        continue
    base = os.path.split(f)[0]
    markedbad = False
    if f.endswith('.gz'):
        h = gzip.open(f,'rb')
    else:
        h = open(f,'rb')
    tocheck = []
    try:
        # Find initial block with sourceFile(s) in it...
        h.seek(0)
        buf = ""
        while '<spectrumList ' not in buf and '<chromatogramList ' not in buf:
            buf += h.read(1024)
        for m in re.finditer(r'<sourceFile .*?</sourceFile>',buf,re.DOTALL):
            sfblock = m.group(0)
            sfm = re.search(r'<sourceFile .* name="([^"]*)"',sfblock)
            sha1m = re.search(r' name="SHA-1" value="([0-9a-fA-F]+)"',sfblock)
            if sfm and sha1m:
                sf = sfm.group(1)
                sf = os.path.split(sf)[1]
                sf = os.path.join(base,sf)
                sha1 = sha1m.group(1).lower()
                tocheck.append((sf,sha1))
    except IOError:
        markedbad = True
        if verbose:
            print("Checking", f, "BAD (can't find sourceFile block.)", file=sys.stderr)
        else:
            print(f)
        raise
    if markedbad:
        continue
    bufsize = 200
    buf = ""
    if f.endswith('.gz'):
        try:
            buf,size = getlastblock(h,bufsize)
        except IOError:
            markedbad = True
            if verbose:
                print("Checking",f, "BAD", file=sys.stderr)
            else:
                print(f)
            raise
    else:
        h.seek(-bufsize,2)
        buf = h.read()
        size = os.path.getsize(f)
    if markedbad:
        continue

    m = re.search(r'<(fileChecksum|sha1)>([a-fA-F0-9]*)</(fileChecksum|sha1)>',buf,re.MULTILINE)
    if m:
        filehash = m.group(2).lower()
        ed = size-bufsize+m.start()+len(m.group(1))+2
        tocheck.append((f,filehash,0,ed))
    else:
        if verbose:
            print("Checking",f, "BAD (no checksum in XML)", file=sys.stderr)
        else:
            print(f)
    h.close()

    for f in tocheck:
        # print "\t".join(map(str,f))
        if verbose:
            print("Checking",f[0], end=' ', file=sys.stderr)
        if not checksha1(*f):
            if verbose:
                print("BAD", file=sys.stderr)
            else:
                print(f[0])
        elif verbose:
            print("GOOD", file=sys.stderr)
