#!/bin/bash
R="$1"
RBGEN="${2}"
if [[ ${R} == '' ]]; then
  R='R'
fi
if [[ ${RBGEN} == '' ]]; then
  RBGEN=build/R/rbgen
fi

${R} CMD INSTALL ${RBGEN}
${R} --vanilla << HERE_DOC
library( rbgen )

## Test we can load data
D = bgen.load( "example/example.16bits.bgen", ranges = data.frame( chromosome = '01', start = 0, end = 100000 ))
str( D )
head( D\$variants )
D\$data[1,1:10,1:3]

## Test we can load data on a subset of samples
samples = c( "sample_001", "sample_102", "sample_050", "sample_499" )
E = bgen.load( "example/example.16bits.bgen", ranges = data.frame( chromosome = '01', start = 0, end = 100000 ), samples = samples )
stopifnot( length( which( E\$data != D\$data[,samples,] ) ) == 0 )
stopifnot( length( which( D\$variants != E\$variants )) == 0 )

## Test we can load data for some specific rsid
rsids = c( "RSID_10", "RSID_20", "RSID_171", "RSID_9" )
F = bgen.load( "example/example.16bits.bgen", rsids = rsids )
str( F )
head( F\$variants )
F\$data

## Check we can load complex data
D = bgen.load(  "example/complex.bgen", data.frame( chromosome = '01', start = 0, end = 1000000 ), max_entries_per_sample = 50 )
str(D)

HERE_DOC
