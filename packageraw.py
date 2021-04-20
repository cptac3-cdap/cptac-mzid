#!/usr/bin/env python3

import sys, os, os.path, subprocess, shutil, encodings
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-r','--raw',type='string',default='wiff',dest='raw',
                  help='RAW file extension to convert from. Default: wiff.')
# parser.add_option('-x','--xml',type='string',default='mzML',dest='extn',
#                 help='XML format to convert to. Default: mzML.')
parser.add_option('-R','--remove',action='store_true',dest='remove',
                  help='Remove zipped files.')
parser.add_option('-v','--verbose',action='store_true',dest='verbose',
                  help='Verbose output from zip.')
opts,args = parser.parse_args()

opts.raw = opts.raw.lower()

zip_prog = 'zip'
zip_args = ['-9','-r']

if not opts.verbose:
    zip_args.append('-q')

curdir = os.path.abspath(os.getcwd())
for root, dirs, files in os.walk(curdir):
    for f in (files+dirs):
        try:
            b,e = f.rsplit('.',1)
        except ValueError:
            continue
        if e.lower() != opts.raw:
            continue
        allfs = []
        for f1 in (files+dirs):
            if f1.startswith(f) and not f1.endswith('.zip'):
                allfs.append(f1)
        if len(allfs) <= 1 and not os.path.isdir(os.path.join(root,f)):
            continue
        if os.path.isdir(os.path.join(root,f)):
            dirs.remove(f)
        cmd = [zip_prog]+zip_args+[f+'.zip']+allfs
        if opts.verbose:
            print(' '.join([s if ' ' not in s else '"'+s+'"' for s in cmd]))
            sys.stdout.flush()
        subprocess.call(cmd,cwd=root,shell=('win' in sys.platform))
        if opts.remove:
            for f1 in allfs:
                if os.path.isdir(os.path.join(root,f)):
                    shutil.rmtree(os.path.join(root,f),ignore_errors=True)
                else:
                    try:
                        os.unlink(os.path.join(root,f1))
                    except OSError:
                        pass
