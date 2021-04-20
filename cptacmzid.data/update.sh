#!/bin/sh
VER="$1"
if [ "$VER" != "" ]; then
   VER="-$VER"
fi
DIR=`dirname "$0"`
cd $DIR
URL=http://cptac-cdap.georgetown.edu.s3-website-us-east-1.amazonaws.com
OS=`uname -o`
AR=`uname -m`
case $OS in
  Cygwin) MACH="win32"; EXT=".zip"; EXE=".exe";;
       *) MACH="linux-$AR"; EXT=".tgz"; EXE=".sh";;
esac
if [ -f ./psmextract.py ]; then
    MACH="src"; EXT=".tgz"; EXE=".py"
fi
ZIP=cptac-mzid$VER.$MACH$EXT
rm -f $ZIP
# Assume wget and unzip/tar are on the path...
wget --no-check-certificate -O $ZIP $URL/$ZIP
if [ "$EXT" = ".zip" ]; then
    unzip -o -d .. $ZIP
else
    tar -C .. -xvzf $ZIP
fi
./psmextract$EXE --version
rm -f $ZIP
