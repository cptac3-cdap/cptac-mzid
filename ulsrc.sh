#!/bin/sh
VER=`python3 psmextract.py --version | tr -d -c '0-9.'`
XX="src"
mkdir -p build dist
python3 ./tosrc.py build/cptac-mzid-${VER}.${XX} psmextract.py checkMzMLsha1.py msconvertall.py packageraw.py
mv build/cptac-mzid-${VER}.${XX}.tgz dist
