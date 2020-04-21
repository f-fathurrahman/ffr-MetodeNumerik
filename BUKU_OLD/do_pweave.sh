#!/bin/bash
BASNAM=`basename $1 .texw`
pweave -f texminted $1 -o GEN_${BASNAM}.tex
