#!/usr/bin/env python3

import sys, os, shutil, stat

def rmminusrf(top):
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))
    try:
        os.rmdir(top)
    except OSError:
        pass

base = sys.argv[1]
base0 = os.path.split(base)[1].split('-')[0]

if os.path.exists(base):
    rmminusrf(base)

base1 = base
base = base+'/'+base0

files = list(sys.argv[2:])

from cx_Freeze import setup, Executable

# Dependencies are automatically detected, but it might need fine tuning.
build_exe_options = {"packages": ["os","encodings.ascii","encodings.utf_8","encodings.utf_16_le","encodings.latin_1","encodings.string_escape","encodings.hex_codec","encodings","_strptime"], "excludes": ["tkinter"]}

sys.argv[1:] = ['build']
setup(executables = [Executable(f, base=None) for f in files])
shutil.copytree('build/exe.linux-x86_64-3.6',base)

for f in files:
    if os.path.exists(f[:-3] + '.ini'):
        shutil.copy(f[:-3]+'.ini',os.path.join(base,f[:-3]+'.ini'))

datadir = "%s.data"%(base0,)
if os.path.isdir(datadir):
    for f in os.listdir(datadir):
        if os.path.isdir(os.path.join(datadir,f)):
            if f != '.svn':
                shutil.copytree(os.path.join(datadir,f),os.path.join(base,f))
        else:
            shutil.copy2(os.path.join(datadir,f),base)
            # if f.endswith('.exe') or f.endswith('.sh'):
            #        os.chmod(os.path.join(base,f),stat.S_IEXEC|stat.S_IRWXU)

for root, dirs, files in os.walk(base):
    dirs1 = []
    for d in dirs:
        if d == '.svn':
            rmminusrf(os.path.join(root,d))
            # os.rmdir(os.path.join(root,d))
        else:
            dirs1.append(d)
    dirs = dirs1

os.system("tar -czf %s.tgz -C %s %s"%(base1,base1,base0))
