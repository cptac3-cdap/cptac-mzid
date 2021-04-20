#!/usr/bin/env python3

import sys, os, os.path, subprocess, shlex, glob, encodings
from configfile import readconfig
from optparse import OptionParser

cfg = readconfig(sys.argv[0])

parser = OptionParser()
parser.add_option('-r','--raw',type='string',default='RAW',dest='raw',
                  help='RAW file extension to convert from. Default: RAW.')
parser.add_option('-v','--verbose',action='store_true',dest='verbose',
                  help='Verbose output from msconvert.')
opts,args = parser.parse_args()

opts.raw = opts.raw.lower()

iniKey = opts.raw

d = os.path.abspath(os.path.split(sys.argv[0])[0])
msconvert_path = os.path.join(d,cfg.get(iniKey,'Path'))
msconvert_prog = cfg.get(iniKey,'Binary')
msconvert_prog = os.path.join(msconvert_path,msconvert_prog)
msconvert_args = shlex.split(cfg.get(iniKey,'Arguements'))

if opts.verbose:
    msconvert_args.append('--verbose')

curdir = os.path.abspath(os.getcwd())
for root, dirs, files in os.walk(curdir):
    for f in (files + dirs):
        try:
            b,e = f.rsplit('.',1)
        except ValueError:
            continue
        if e.lower() != opts.raw:
            continue
        if glob.glob(os.path.join(root,b+'*'+cfg.get(iniKey,"Extn"))):
            continue
        # print [msconvert_prog,f]+msconvert_args
        if os.path.isdir(os.path.join(root,f)):
            dirs.remove(f)
        subprocess.call([msconvert_prog,f]+msconvert_args,cwd=root,shell=True)
