#!/bin/sh
BINPATH=`dirname $0`
cd $BINPATH
if [ "$1" = "" ]; then
  exec ./quickbuild.sh -q --incremental --without-mz5 --without-binary-msdata pwiz_tools/examples//textpsm2mzid pwiz_tools/commandline//idconvert pwiz_tools/commandline//msaccess pwiz_tools/commandline//msconvert pwiz_tools/examples//write_mzid_example_files
else
  exec ./quickbuild.sh -q --incremental --without-mz5 --without-binary-msdata $*
fi
