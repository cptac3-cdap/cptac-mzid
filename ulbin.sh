#!/bin/sh
VER=`python3 version.py`
OS=`uname`
AR=`uname -m`
XX="linux-$AR"
mkdir -p build dist
./tobin.py build/cptacmzid-${VER}.${XX} psmextract.py checkMzMLsha1.py msconvertall.py packageraw.py mzml.py featxml.py version.py
mv build/cptacmzid-${VER}.${XX}.tgz dist

