#!/bin/sh
VER=`python3 version.py`
XX="src"
mkdir -p build dist
python3 ./tosrc.py build/cptacmzid-${VER}.${XX} psmextract.py checkMzMLsha1.py msconvertall.py packageraw.py
mv build/cptacmzid-${VER}.${XX}.tgz dist
