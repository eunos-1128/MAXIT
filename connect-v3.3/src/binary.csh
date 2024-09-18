#!/bin/tcsh -f
#
if ( ! ( -d data ) ) then
    ln -s data-stl data
endif
#
setenv RCSBROOT `pwd`
#
./utillib-v1.1/bin/cif2bin
#
./connect-v3.3/bin/connect_main
if ( -e component.odb ) then
    mv component.odb ./data/binary
endif
#
if ( -e variant.odb ) then
    mv variant.odb ./data/binary
endif
#
if ( -e component_all.cif ) then
    rm -f component_all.cif
endif
#
if ( -e variant_all.cif ) then
    rm -f variant_all.cif
endif

