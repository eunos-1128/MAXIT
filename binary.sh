#!/bin/bash
#
export RCSBROOT="$(pwd)"
#
./bin/cif2bin
#
./bin/connect_main
if [ -e component.odb ]
then
    mv component.odb ./data/binary
fi
#
if [ -e variant.odb ]
then
    mv variant.odb ./data/binary
fi
#
if [ -e component_all.cif ]
then
    rm -f component_all.cif
fi
#
if [ -e variant_all.cif ]
then
    rm -f variant_all.cif
fi

