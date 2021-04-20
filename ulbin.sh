#!/bin/sh
VER=`python3 psmextract.py --version | tr -d -c '0-9.'`
OS=`uname`
AR=`uname -m`
XX="linux-$AR"
mkdir -p build dist
./tobin.py build/cptacmzid-${VER}.${XX} psmextract.py checkMzMLsha1.py msconvertall.py packageraw.py
mv build/cptacmzid-${VER}.${XX}.tgz dist

