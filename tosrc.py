#!/usr/bin/env python3

import sys, os, shutil, stat, os.path

def rmminusrf(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    os.rmdir(top)

def clean(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            if name.endswith('.pyc'):
                os.remove(os.path.join(root, name))
            if name.endswith('~'):
                os.remove(os.path.join(root, name))
            if name.endswith('.so'):
                os.remove(os.path.join(root, name))
        for name in dirs:
            if name == '.svn':
                rmminusrf(os.path.join(root, name))
            if name == 'tests':
                rmminusrf(os.path.join(root, name))
            if name == '__pycache__':
                rmminusrf(os.path.join(root, name))

base = sys.argv[1]
base0 = os.path.split(base)[1].split('-')[0]

if os.path.exists(base):
    rmminusrf(base)

base1 = base
base = base+'/'+base0

try:
    os.makedirs(base)
except OSError:
    pass
for f in sys.argv[2:]:
    shutil.copyfile(f,base+'/'+f)
    os.chmod(os.path.join(base,f),stat.S_IEXEC|stat.S_IRWXU)
for pkg in ['peptidescan','parsimony.py','mzml.py','StringIndex.py','find_elementtree.py','psmformat.py','psmparse.py','version.py']:
    pdir,pname = os.path.split(pkg)
    if os.path.isdir(pkg):
        shutil.copytree(pkg,os.path.join(base,pname))
    else:
        shutil.copyfile(pkg,os.path.join(base,pname))

datadir = "%s.data"%(base0,)
if os.path.isdir(datadir):
    for f in os.listdir(datadir):
        if os.path.isdir(os.path.join(datadir,f)):
            shutil.copytree(os.path.join(datadir,f),os.path.join(base,f))
        else:
            shutil.copy2(os.path.join(datadir,f),os.path.join(base,f))

clean(base)

os.system("tar -czf %s.tgz -C %s %s"%(base1,base1,base0))
